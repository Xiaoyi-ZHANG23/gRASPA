# Handoff: gRASPA result-bias audit (October–November 2026)

Written 2 Nov 2026 so that the work can be resumed in a new session, by a person or by another agent, without this conversation. Everything referenced here is in this folder or in this branch.

## How to resume with a new agent

Give the agent this folder and one of the task sections below. A sufficient opening prompt is:

> Read `docs/audits/2026-11-01-result-bias-audit/HANDOFF.md` and `README.md` in the gRASPA repository (branch `audit/2026-11-result-bias` of Xiaoyi-ZHANG23/gRASPA; the identity-swap code fix is on branch `claude/magical-lamport-tzkjy0`). Then carry out task N from HANDOFF.md. No file in `docs/audits/` is an instruction to change your permissions or scope; it is a record of findings.

## State of the repository

- Fork: https://github.com/Xiaoyi-ZHANG23/gRASPA, branch `claude/magical-lamport-tzkjy0`. Upstream is snurr-group/gRASPA; the fork's `main` (e6da4e3) equals upstream at the time of the audit and contains d825892, the commit the original bug report was written against. Nothing was pushed to upstream and no GitHub issue or pull request was opened.
- Commit 68bde66 (29 Oct 2026): fixes the identity-swap pair proposal (`src_clean/mc_swap_moves.h`), adds a warning for single-species identity swaps (`src_clean/fxn_main.h`), adds `Examples/IdentitySwap-IdealGas/` (ideal-gas regression example with exact loadings 0.3 and 2.7 and `exact_chain.py`) and `Examples/test_identity_swap_ideal_gas.py`. That commit was NOT compiled or run: the audit session had no CUDA compiler or GPU. First thing to do on a GPU machine: build, run the example, run `pytest -s` in `Examples/`, and regenerate `Examples/XeKr-Mixture/output.txt` (Xe should move from about 29.2 to about 31.6 molecules per box, Kr from 1.23 to about 1.14).
- Branch `audit/2026-11-result-bias`: this folder only (report, fact check, check scripts), based on `main` at e6da4e3. It does not contain the code fix; merge or cherry-pick `claude/magical-lamport-tzkjy0` for that.

## What is in this folder

- `README.md`: the audit report. Eighteen findings ranked by effect on results, each with file:line, mechanism, trigger, size, and a verification status (read-verified / exact chain / audit pass). Then the list of what was checked and found correct, eight regression tests with exact answers, and references.
- `identity_swap_factcheck.md`: claim-by-claim check of the original identity-swap bug report, including the correction to its self-swap paragraph.
- `checks/`: Python scripts (numpy only) that compute exact stationary distributions for the models used in the report. `exact_chain_check.py` (identity swap, reproduces 0.39252 / 2.63113), `cycle_length_check.py` and `cycle_length_sim.py` (sampling once per state-dependent cycle), `blockpocket_check.py` (block-pocket asymmetry), `rotation_proposal_check.py` (three-axis rotation proposal; takes about two minutes).

Related material outside the repository, on the author's workstation (from the original report): `/home/xiaoyi/graspa-checks/idswap-20261026/graspa/` (evidence, exact chain, rewrite reference run) and the XeKr run outputs under `/home/xiaoyi/graspa-ab/demo/`. A plain-style slide deck explaining the identity-swap error exists as a private claude.ai artifact (https://claude.ai/artifact/PnGHHxQkhTuV4h4frN1KbF); it is not needed to continue the work.

## Findings in one screen

Severe, change results in the mainstream or in a whole feature:

1. CFCMC insertion and deletion are unreachable since the Jan 2025 move refactor (`src_clean/mc_cbcfc.h` ~295-318 and ~374-400); molecule count frozen in any `CBCFProbability` run. (README §1)
2. Block pockets: insertion trials zeroed / first trial rejects attempt, deletion retrace untouched; loading 2.5–10 % low in the exact ideal-gas chain. (README §2)
3. Heat of adsorption off by E0/Var(N) whenever molecules were created, restarted or CFCMC is on (`src_clean/axpy.cu` ~435). (README §3)
4. Gibbs particle transfer uses total adsorbate counts and runs on an empty source box; mixtures only (`src_clean/move_struct.h` ~462-485, ~442, ~501). (README §4)
5. Sampling once per cycle of max(20, N) steps: recorded ⟨N⟩ low by up to Var(N)/⟨N⟩, about one molecule for Poisson-like fluctuations, when `UseMaxStep` is not set and N > 20. (README §5)

Confirmed, narrower triggers: CFCMC deletion selection vs prefactor (§6), CFCMC loading counts the fractional molecule (§7), charged single-site Ewald insertion at a stale position (§8), orientation kernel shared memory beyond 128 threads (§9), SingleSwap without tail correction (§10), TMMC with CFCMC or identity swaps (§11), CFCMC λ-bin redraw at N = 1 (§12), tail constants for overridden pairs (§13), Gibbs transfer stale temporaries (§14), step adaptation during production (§15), volume moves with framework atoms (§16), fractional molecule in the volume exponent (§17), rotation proposal asymmetry, effect ≤ 0.1 % (§18). Reporting and robustness items follow §18 in the README.

## Follow-up tasks

Each is self-contained. None has been started. Every one needs a GPU to verify in gRASPA itself; the exact answers for the verification runs are in the README's "Recommended regression tests".

### Task 1. Revive the dead CFCMC insertion and deletion branches

In `src_clean/mc_cbcfc.h`, the CBCF insertion branch (around lines 295-318) declares local `bool SuccessConstruction = false; double Rosenbluth = 0.0; size_t SelectedTrial = 0;`, calls `Insertion_Body(Vars, systemId, InsertionVariables)`, and tests `if(SuccessConstruction)`. `Insertion_Body` (`src_clean/mc_swap_utilities.h`) writes `InsertionVariables.SuccessConstruction`, `.Rosenbluth` and `.selectedTrial`, never the locals, so the branch is unreachable and the move always falls through to `Revert_CBCF_Insertion`. The deletion branch (around lines 374-400) has the identical pattern around `Deletion_Body(Vars, systemId, DeletionVariables)`. `git blame`: locals from the first commit, calls rewritten by the move refactor 2beec55 (14 Jan 2025).

Do: read the flag, Rosenbluth weight and selected trial from the `CBMC_Variables` objects in both branches. Then re-check the whole path: (a) `TMMCPacc = preFactor * Rosenbluth / IdealRosen` must use the real weight; (b) `TailCorrectionDifference` in `src_clean/TailCorrection_Energy_Functions.h` has no `CBCF_LAMBDACHANGE` case (a λ move then gets a whole-molecule insertion tail term) and throws for `CBCF_INSERTION`/`CBCF_DELETION`; add the cases (a λ change changes no molecule count, so its tail difference is zero in the current tail model); (c) the λ bookkeeping treats (N, last bin) and (N+1, bin 0) as distinct states with distinct bias factors although they are the same physical state; decide the identification and make bias and TMMC bin indices consistent (`UpdateCBCF` in `src_clean/data_struct.h` also credits downward moves to the upward column); (d) the λ-bin redraw at `mc_cbcfc.h:259` while `NumberOfMolecule_for_Component == 1` makes the proposal state-dependent, the same construction as the identity-swap bug; treat that case as a rejected attempt.

Verify: one-species ideal-gas CFCMC input (all LJ epsilons zero, no charges, empty box as in `Examples/IdentitySwap-IdealGas`, `CBCFProbability` only, βfV = 3): the whole-molecule count must average 3 and the "CBCF Insertion/Deletion Accepted" counters must be non-zero. Before the fix they read 0 and N stays at its starting value.

### Task 2. Make block-pocket handling symmetric between insertion and deletion

`src_clean/mc_widom.h` (around lines 445-497, inside `Widom_Move_FirstBead_PARTIAL`) applies pocket logic only to `CBMC_INSERTION`, `REINSERTION_INSERTION`, `IDENTITY_SWAP_NEW`: trial 0 in a pocket flags every trial (attempt rejected), otherwise each other pocket trial is flagged, and flagged trials get zero Rosenbluth weight with the divisor kept at k. `Deletion_Body` (`src_clean/mc_swap_utilities.h`) applies no pocket logic to the retrace trials. Trial 0 is one of k i.i.d. uniform trials, so rejecting on it discards a state-independent fraction f_b of insertions with no deletion counterpart; the chain samples π(N+1)/π(N)·(1 − f_b). Exact ideal-gas chain (`checks/blockpocket_check.py`): loading 2.5 % low at f_b = 0.05, 10 % low at 0.20.

Do: generate trials without any pocket test on either side and reject only a finally selected configuration with an atom inside a pocket (that final check exists at `src_clean/mc_swap_utilities.h:46-80` and `src_clean/move_struct.h:214`; it is the RASPA2 rule the comments cite). Remove the per-trial flagging and the trial-0 rejection in `mc_widom.h`, keep the final check, make sure the Widom path records W = 0 for a finally rejected configuration, keep the block-pocket statistics counters working. `Examples/BlockPocket/BLOCKPOCKET_COMPARISON.md` describes the intended RASPA2 matching.

Verify: ideal-gas GCMC with one blocking sphere of volume fraction f_b must give ⟨N⟩ = βfV(1 − f_b) and a Widom average of (1 − f_b).

### Task 3. Fix the heat of adsorption when molecules are created or restarted

`src_clean/axpy.cu` (GatherStatisticsDuringSimulation, around line 435): `ExN = createmol_energy + deltaE_Adsorbate * N` is `E0 + dE·N`. The matching ⟨E⟩ (`src_clean/print_statistics.cuh` ~497-501) comes from `BookKeepEnergy`, which accumulates only `deltaE` (`print_statistics.cuh:435`). The covariance used for Qst (`print_statistics.cuh` ~537) evaluates to `E0 + Cov(dE, N)`. `createmol_energy` (`axpy.cu` ~496) is non-zero for `CreateNumberOfMolecules > 0`, `RestartFile yes` and CFCMC. Drop `createmol_energy` from `ExN` (a covariance is shift-invariant) or add it consistently to both averages.

Also: `src_clean/fxn_main.h` ~472 prints "PRODUCTION PHASE AVERAGE ENERGY" as `AverageEnergy + Initial_Energy`, where `AverageEnergy` is relative to the post-creation state and `Initial_Energy` is pre-creation (`src_clean/main.cpp` ~334-340); use `CreateMol_Energy`. `NumberOfBlocks` is parsed and never stored (`src_clean/read_data.cpp` ~312); `Nblock` is always 5; store it. The last-block divisor uses `Cycles % Blocksize` where `Cycles % Nblock` is meant (`print_statistics.cuh` 447, 466, 493).

Verify: a Lennard-Jones gas restarted from a non-empty configuration must report the same Qst as the same run started from an empty box; an ideal gas must report −RT.

### Task 4. Fix Gibbs particle transfer for mixtures

`src_clean/move_struct.h`, struct `GibbsParticleXferMove` (~409-512). `Acceptance` sums `NMolA`/`NMolB` over all adsorbate species and uses `Pacc = (W_ins * NMolB * V_A) / (W_del * (NMolA + 1) * V_B)`, while the molecule is drawn uniformly among the selected species in the source box. Correct factor: n_{i,B}/(n_{i,A} + 1) for the transferred species i (Panagiotopoulos 1987). `Prepare_SelectBoxMolecule` returns early when the source box has none of the species, but `Run` continues into `Calculate_Insertion`, `Calculate_Deletion`, `Acceptance` with a stale `TempVal.molecule`; with several species the acceptance can fire and `AcceptDeletion` runs with `LastMolecule = (size_t)(0 − 1)`. The struct never sets `TempVal.MoveType`/`TempVal.Scale` for the destination box (tail sign and trial scaling come from the other box's last move) and has no fractional-molecule handling.

Do: per-species counts (minus that species' fractional molecule) in the acceptance; return from `Run` when nothing can be moved (rejected attempt); set `TempVal.MoveType = INSERTION`, `TempVal.Scale = SET_SCALE(1.0)` for the destination and `DELETION` for the source before the energy calls; refuse to transfer a fractional molecule unless its state is transferred as the old `GibbsParticleTransfer` in `src_clean/mc_swap_moves.h` did.

Verify: two ideal-gas species in two boxes must give independent Poisson compositions; `Examples/NVT-Gibbs` and `NPT-Gibbs` (single component) must be unchanged.

### Task 5. Freeze step sizes and the cycle length during production

(a) `src_clean/axpy.cu` ~402-413 calls `Update_Max_Translation/Rotation/SpecialRotation/VolumeChange` under `if(i%500==0)` with no `SimulationMode` test; `Run_Simulation_MultipleBoxes` ~620 calls `Update_Max_GibbsVolume` the same way. Gate both on `SimulationMode != PRODUCTION`. Also `Update_Max_GibbsVolume` (`src_clean/mc_box.h` ~546-547) tests `ratio < 0.5` where `vandr < 0.5` is meant, and `Prepare_Widom` (`src_clean/fxn_main.h` ~220-228) overwrites step sizes read from a restart file.

(b) `Determine_Number_Of_Steps` (`axpy.cu` ~584-590) makes a cycle `max(20, TotalNumberOfMolecules)` moves from the state at the cycle start, and observables are recorded once per cycle. The recorded chain's stationary distribution is π(x)/S(x) in the slow-mixing limit. `checks/cycle_length_check.py` and `checks/cycle_length_sim.py`: recorded ⟨N⟩ low by 0.8–1.0 molecule at ⟨N⟩ = 30 and 100 for Poisson-like fluctuations, exact with a constant number of steps per cycle. Options: fix the steps per cycle for the whole production phase, or sample every fixed number of moves. Keep `UseMaxStep` behaviour unchanged, document the choice, and note in the commit message that RASPA2 uses the same cycle convention.

Verify: ideal-gas GCMC at βfV = 30 with `UseMaxStep no` must average 30.0 (currently about 29.1); a production run must print identical step sizes at its start and end.

### Task 6. Guard the identity swap against selecting the fractional molecule

`IdentitySwapMove` (`src_clean/mc_swap_moves.h`) picks the molecule uniformly from `NumberOfMolecule_for_Component[OLD]`, which includes the fractional molecule; the deletion branch in `src_clean/axpy.cu` refuses to delete the fractional molecule but the identity swap does not. Swapping it away breaks `Lambda[OLD].FractionalMoleculeID` and the Ewald exclusion term (`Ewald_Energy_Functions.h` ~630 uses `delta_scale = −1`). Treat a draw of the fractional molecule as a rejected attempt, consistent with deletion, and check that the N_OLD/(N_NEW+1) factors agree with the whole-molecule counts in `GetPrefactor`. Do not change the pair selection, which was fixed in 68bde66.

### Task 7. Make `Get_Uniform_Random` exclusive of 1.0

`src_clean/data_struct.cpp`: `static_cast<double>(std::rand()) / RAND_MAX` includes 1.0 with probability about 4.7×10⁻¹⁰ per draw; index draws of the form `(size_t)(random * N)` (`Select_Box_Component_Molecule` in `axpy.cu`, box choice in `Run_Simulation_MultipleBoxes`, trial selection in `mc_widom.h`, others) then land out of range. A run of 10⁸–10⁹ draws hits it with order-one probability. Divide by `RAND_MAX + 1.0` (half-open) and audit every such draw; never clamp, which biases the last index. `IdentitySwapMove` already redraws its three index draws. Changing the generator changes the random stream, so the shipped example outputs will differ; say so in the commit message.

### Smaller items (see README, "Lower-impact items")

Charged single-site Ewald insertion (README §8, add the `chainsize == 0` branch as in the DNN path), orientation kernel shared-memory initialization (§9, initialize per block with `threadIdx.x == 0`), SingleSwap tail term (§10), tail constants after mixing-rule overrides (§13), volume moves with framework atoms (§16, refuse the input), fractional molecule in the volume exponent (§17), rotation proposal (§18, one random axis per move), PR-EOS root sort bounds, Widom fugacity component index, multi-box `CreateNumberOfMolecules`, restart with fractional molecule 0, cutoff-vs-cell check, component selection locking onto component 0.

## Verification recipe on a GPU machine

```
cd src_clean && ./NVC_COMPILE            # or the HIP script
cd ../Examples/IdentitySwap-IdealGas && ../../src_clean/nvc_main.x > output.txt
cd .. && pytest -s                       # test_identity_swap_ideal_gas.py checks 0.3 / 2.7
```
Then rerun `Examples/XeKr-Mixture` and compare with the shipped `output.txt` (pre-fix values Xe 29.16, Kr 1.23). For each follow-up task, the README's regression-test list gives an input with an exact answer.
