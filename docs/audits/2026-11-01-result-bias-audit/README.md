# gRASPA audit: bugs that shift simulation results

Date: 1 Nov 2026. Code: Xiaoyi-ZHANG23/gRASPA, branch `claude/magical-lamport-tzkjy0` at 68bde66 (= upstream main e6da4e3 plus the identity-swap fix). All paths are under `src_clean/`. Nothing in the repository was changed by this audit.

Method. The code was read in five areas (swap/CFCMC/Gibbs moves, CBMC machinery, single-particle and volume moves, thermodynamic inputs and statistics, energy functions). Each finding below carries a verification status: **read-verified** means I read the quoted lines myself and the argument is complete; **exact chain** means the size was computed with an exactly solvable Markov chain (scripts in `checks/` next to this file: `cycle_length_check.py`, `cycle_length_sim.py`, `blockpocket_check.py`, `rotation_proposal_check.py`); **audit pass** means a reviewer reported it with code quotes and I checked the quotes but did not re-trace every path. No GPU or CUDA compiler was available, so nothing was run.

Already known, not repeated: the identity-swap pair proposal (fixed in 68bde66), `Get_Uniform_Random()` returning exactly 1.0, the identity swap selecting a fractional molecule.

## Ranked summary

| # | Finding | Trigger | Effect on results | Status |
|---|---|---|---|---|
| 1 | CFCMC insertion and deletion can never be accepted (dead branches since the Jan 2025 move refactor) | `CBCFProbability > 0` | Molecule count frozen at its starting value for the whole run; any CFCMC isotherm point is wrong | read-verified |
| 2 | Block pockets: insertion trials in pockets are zeroed and an unlucky first trial rejects the whole insertion, deletion retrace has no pocket logic | `UseBlockPockets` | Detailed balance broken; ideal-gas chain: loading low by 2.5 % (5 % blocked volume) to 10 % (20 %); extra loading-dependent term from the trial-weight asymmetry | read-verified, exact chain |
| 3 | Heat of adsorption: creation energy added inside <E·N> but not <E> | `CreateNumberOfMolecules > 0`, `RestartFile yes`, or CFCMC | Qst shifted by E0/Var(N), typically tens to hundreds of kJ/mol more exothermic for restarted or pre-loaded runs | read-verified |
| 4 | Gibbs particle transfer uses total adsorbate counts, and runs on an empty source box | Gibbs with 2+ adsorbate species | O(1) composition bias toward pure boxes; memory corruption when the source box lacks the species | read-verified |
| 5 | Observables sampled once per cycle of max(20, N) steps | Default cycle definition (`UseMaxStep no`) with N > 20 | <N> low by up to Var(N)/<N>, about 1 molecule for Poisson-like fluctuations; exact chain gives −0.8 to −1.0 at <N> = 30 and 100 | read-verified, exact chain + simulation |
| 6 | CFCMC whole-molecule deletion picks among all N molecules but the prefactor uses N−1 | `CBCFProbability` + `SwapProbability` | Chain samples (n+1)·π(n): <N> high by about Var(N)/(<N>+1), up to +1 molecule | read-verified |
| 7 | CFCMC loading counts the fractional molecule as a whole molecule | `CBCFProbability > 0` | Reported loading high by 0.5 to 1 molecule per box and species | read-verified |
| 8 | Single-site charged adsorbate: CBMC insertion evaluates the Ewald Fourier term at a stale position | Monatomic ion with `ChargeMethod Ewald` | Wrong acceptance; visible as energy drift | read-verified |
| 9 | Trial-orientation kernel initializes shared first-bead data from global thread 0 only | `NumberOfTrialOrientations × (atoms−1) > 128` | Garbage trial orientations for large rigid molecules; no shipped example crosses the limit (p-xylene: 110 threads) | read-verified |
| 10 | SingleSwap acceptance omits the tail correction | `TURN_OFF_CBMC_SWAP yes` + tail corrections | Samples the untailed potential; loading low | read-verified |
| 11 | CFCMC + TMMC collection matrix and bias wrong in several ways; TMMC + identity swap bypasses the matrix | `RunTMMC yes` with CBCF or identity swaps | ln Π(N) wrong | audit pass, key lines read-verified |
| 12 | CFCMC λ-bin redraw when only the fractional molecule exists | CFCMC in the Henry regime | Same proposal-asymmetry class as the identity-swap bug; factor up to 2 on the N = 1 ↔ 2 flow (moot until #1 is fixed) | read-verified |
| 13 | Tail-correction constants not recomputed for pairs overridden in "mixing rules to overwrite" when the general tail rule is `yes` | That input combination | Tail uses Lorentz–Berthelot ε, σ while the pair uses the override; invisible to the drift check | audit pass |
| 14 | Gibbs transfer reuses the other box's stale `TempVal.MoveType`/`Scale` | Gibbs + tail corrections + swap or CFCMC moves per box | Wrong tail sign or fractional scaling on the inserted molecule | audit pass |
| 15 | Maximum step sizes (translation, rotation, special rotation, NPT and Gibbs volume) keep adapting every 500 cycles during production | Every run with production cycles ≥ 500 | Proposal depends on the chain's history; expected small for the ±5 % translation/rotation updates, larger for the ×0.5–1.5 volume updates; also drives `MaxRotation` to its 3.14 rad clamp in dilute systems | read-verified |
| 16 | Volume moves evaluate framework–guest energies from an unwritten half of the framework arrays | NPT or Gibbs volume moves in a box that contains framework atoms | Garbage host–guest energies; inert for bulk boxes, which is the supported use | audit pass |
| 17 | NPT/Gibbs volume acceptance exponent excludes the fractional molecule although its coordinates are scaled | CFCMC + volume moves | P(V) short by one factor of V, about 1/(N+1) relative on the volume | audit pass |
| 18 | Rotation proposal is Rz(c)·Ry(b)·Rx(a) with independent uniform angles, not inversion-symmetric on SO(3) | Every rotation move | Measured effect ≤ 0.1 % on an exactly solvable rotor observable at all step sizes; low practical impact, cheap to fix | read-verified, Metropolis test |

Lower-impact items (reporting, robustness, edge cases) are listed after the details.

## Details

### 1. CFCMC insertion and deletion are dead code

`mc_cbcfc.h:295-318` (insertion branch) declares `bool SuccessConstruction = false;` and `double Rosenbluth = 0.0;` as locals, calls `Insertion_Body(Vars, systemId, InsertionVariables)`, which writes `InsertionVariables.SuccessConstruction` and `InsertionVariables.Rosenbluth`, then tests `if(SuccessConstruction)` on the local that was set to false on line 307. The branch is unreachable; the move always takes `Revert_CBCF_Insertion`. The deletion branch at `mc_cbcfc.h:374-400` has the identical pattern around `Deletion_Body`. Even if the test were fixed, the acceptance uses the local `Rosenbluth = 0.0`, so the probability would be zero.

`git blame`: the locals date from the first commit (debd7e4, Aug 2024); the calls were rewritten by the move refactor 2beec55 (14 Jan 2025), which moved these outputs into `CBMC_Variables` and did not update the CFCMC callers. Every build since then is affected, upstream included.

Consequence: with CBCF moves only, N for that species stays at `CreateNumberOfMolecules + 1` for the entire run while λ random-walks in [0, 1]; the loading reported is the starting value. With CBCF plus whole-molecule swaps, N changes only through the swaps (subject to #6). The quickest runtime check is the `CBCF Insertion/Deletion Accepted` counters, which must read 0. No shipped example uses `CBCFProbability`, which is why it went unnoticed.

Related, same file: `CBCF_LambdaChange` calls `TailCorrectionDifference(..., CBCF_LAMBDACHANGE)`, which has no case for that move type, so with tail corrections on every λ move adds a whole-molecule insertion tail term and the insertion/deletion types throw. When the dead branches are revived, the λ bookkeeping across the N boundary also needs a decision: (N, last bin) and (N+1, bin 0) are treated as distinct states with distinct bias factors although they are the same physical state, so the cycle through them has no reverse edge.

### 2. Block pockets break detailed balance between insertion and deletion

`mc_widom.h:445-497`, executed for `CBMC_INSERTION`, `REINSERTION_INSERTION` and `IDENTITY_SWAP_NEW` only: after the k first-bead trials are drawn, if trial 0 lies in a blocked pocket every trial is flagged (the attempt is rejected); otherwise each other trial in a pocket is flagged. Flagged trials are dropped from the Rosenbluth sum while the divisor stays k (`mc_widom.h:47-52`, `369-370`). For `CBMC_DELETION` nothing of the kind exists: `Deletion_Body` (`mc_swap_utilities.h:135-225`) never calls `BlockedPocket`, so the k−1 random retrace trials that land in pockets keep their full Boltzmann weight.

Two asymmetries follow. (a) Trial 0 of an insertion is one of k i.i.d. uniform positions, so the insertion attempt is discarded with a state-independent probability f_b (blocked volume fraction) with no counterpart on the deletion side; the chain then samples π(N+1)/π(N) × (1 − f_b), as if the fugacity were lowered by (1 − f_b). (b) Insertion weights exclude pocket trials, deletion weights include them. The second check on the finally selected configuration (`mc_swap_utilities.h:46-80`) is the RASPA2 rule and is fine on its own; RASPA2 does not zero individual trials on either side. Reinsertion carries the same factor in both directions and is unaffected.

Exact ideal-gas chain (`blockpocket_check.py`, U = 0, k trials, exact answer a(1 − f_b)):

| a = βfV | f_b | gRASPA <N> | exact | relative |
|---|---|---|---|---|
| 5 | 0.05 | 4.604 | 4.750 | −3.1 % |
| 5 | 0.10 | 4.241 | 4.500 | −5.7 % |
| 5 | 0.20 | 3.599 | 4.000 | −10.0 % |
| 30 | 0.10 | 25.86 | 27.00 | −4.2 % |
| 30 | 0.20 | 22.56 | 24.00 | −6.0 % |

With interactions, asymmetry (b) adds a term whose sign depends on which direction saturates: in the Henry regime only (a) survives (loading low by about f_b(k−1)/k); near saturation, where deletions are the rarely accepted direction, pocket trials inflate W_old and suppress deletions, so the loading can come out high. The same zeroing lowers the Widom average by (1 − f_b(k−1)/k) and raises μ_ex accordingly. The shipped `Examples/BlockPocket` is affected. Fix: either apply no pocket logic to individual trials on either side and keep only the final-configuration rejection (RASPA2), or apply identical logic to the retrace trials.

### 3. Heat of adsorption carries the creation energy as a constant

`axpy.cu:434-436`:
```
double ExN = SystemComponents.createmol_energy + deltaE_Adsorbate * SystemComponents.NumberOfMolecule_for_Component[comp];
```
is `E0 + dE·N`, not `(E0 + dE)·N`. The matching `<E>` in `print_statistics.cuh:497-501` is built from `BookKeepEnergy`, which accumulates `deltaE` only (`print_statistics.cuh:435`). The covariance therefore evaluates to `E0 + Cov(dE, N)` and the printed Qst picks up E0 / Var(N) after unit conversion. `createmol_energy` (`axpy.cu:496`) is the adsorbate energy of the configuration after `CreateMolecule_InOneBox`, so E0 ≠ 0 for `CreateNumberOfMolecules > 0`, for `RestartFile yes` (restart molecules are in `CreateMol_Energy`) and for CFCMC (the fractional molecule is created before the run). Starting from an empty box with plain GCMC gives E0 = 0, which is why the examples look right. Size: 50 pre-loaded molecules at −20 kJ/mol with Var(N) ≈ 5 shifts Qst by about −200 kJ/mol; a restarted production run inherits the full adsorbate energy as E0. A covariance is shift-invariant, so the fix is to drop `createmol_energy` from the product (or add it to both).

### 4. Gibbs particle transfer with mixtures

`move_struct.h:462-485` (`GibbsParticleXferMove::Acceptance`): `NMolA` and `NMolB` are summed over all adsorbate species of the two boxes and
```
Pacc = (W_ins * NMolB * V_A) / (W_del * (NMolA + 1) * V_B);
```
while the molecule to move was drawn uniformly among the n_{i,B} molecules of the selected species i in the source box (`move_struct.h:444`). The correct factor is n_{i,B}/(n_{i,A} + 1). For one adsorbate species the two coincide. For mixtures the chain is reversible for π divided by the multinomial factors of each box's composition, which pushes each box toward a pure species; a binomial factor of order 10^20 for 100 molecules at 50/50 is an O(1) bias on coexistence compositions.

Second defect, same struct: `Prepare_SelectBoxMolecule` returns early when the source box holds no molecule of the species (`move_struct.h:442`), but `Run` (`:501-503`) continues into `Calculate_Insertion` and `Calculate_Deletion` with a stale molecule index. With one species `NMolB = 0` makes `Pacc = 0`; with several species `NMolB > 0`, so a deletion of a non-existent molecule can be accepted and `AcceptDeletion` runs with `LastMolecule = (size_t)(0 − 1)`. Both shipped Gibbs examples are single-component CO2 and unaffected. The struct also has no fractional-molecule handling, unlike the older `GibbsParticleTransfer` function that is no longer called.

### 5. Sampling once per cycle of max(20, N) steps

`axpy.cu:584-590` sets the number of moves in a cycle to `max(20, TotalNumberOfMolecules)` from the state at the start of the cycle (capped by `MaxStepPerCycle` when `UseMaxStep yes`), and `axpy.cu:633-638` records the observables once at the end of each cycle. Every individual move leaves π invariant, but the recorded chain Y_k = X_{T_k} with T_{k+1} = T_k + S(X_{T_k}) has kernel K^{S(x)}(x, ·), and π is not stationary for it when S depends on x. In the slow-mixing limit the recorded distribution is π(x)/S(x): states with more molecules are recorded less often because their cycles are longer. Two-state check: symmetric flip probability p, S(a) = 1, S(b) = 2 gives π′(a)/π′(b) = 2(1 − p).

Exact chain for an ideal gas with insertion/deletion plus N-preserving moves (`cycle_length_check.py`), and a direct simulation (`cycle_length_sim.py`):

| βfV | swap fraction | recorded <N> | bias | fixed 20-step cycles |
|---|---|---|---|---|
| 10 | any | 9.999 | −0.001 | 0 |
| 30 | 1 | 29.23 | −0.77 | 0 |
| 30 | 1/3 | 29.09 | −0.91 (simulation −0.86 ± 0.03) | 0 (simulation −0.03 ± 0.03) |
| 30 | 0.1 | 29.04 | −0.96 | 0 |
| 100 | 1/3 | 99.07 | −0.93 | 0 |

The bias is about −Var(N)/<N> when N decorrelates slowly compared with one cycle, so it is close to one molecule for Poisson-like fluctuations (low loading, large boxes) and small near saturation where Var(N) ≪ <N>; it vanishes for N ≤ 20 because S is then constant. Relative size: −3 % at <N> = 30 with ideal-gas-like fluctuations, −1 % at 100. It also biases every observable correlated with N. Almost all shipped examples set `UseMaxStep yes` with `MaxStepPerCycle 1`, which makes S constant and removes the effect; `TIP4PEW-MgMOF-LJ1264`, `NVT-Gibbs` and `NPT-Gibbs` use the default cycle and are exposed. The cycle convention is RASPA2's, so the same argument applies to any code that samples once per variable-length cycle; I did not verify RASPA2's sampling code in this audit. Remedies: a constant number of steps per cycle in production, or sampling every fixed number of steps.

### 6. CFCMC whole-molecule deletion: selection over N, prefactor with N−1

`axpy.cu:98` draws the molecule uniformly over `NumberOfMolecule_for_Component[comp]`, which includes the fractional molecule; `axpy.cu:247` turns the attempt into a no-op when the fractional molecule is drawn; `mc_utilities.h:322, 338` compute the deletion prefactor with the fractional molecule subtracted. With n whole molecules the per-molecule deletion proposal is 1/(n+1) while the acceptance pair is the standard one for a selection probability 1/n. The chain samples π′(n) ∝ (n+1)·π(n), so <N> is high by about Var(N)/(<N>+1), up to one molecule for Poisson-like fluctuations. Requires `CBCFProbability` together with `SwapProbability` (allowed unless TMMC is on). Given #1, this is currently the only path by which N changes in a CFCMC run. Fix: draw among whole molecules only, or keep the draw and use N_total in both prefactors.

### 7. CFCMC loading counts the fractional molecule

`axpy.cu:430` pushes `NumberOfMolecule_for_Component[comp]` into `MolAverage` every cycle; `data_struct.h:1150` documents that this count includes the fractional molecule; the same count feeds excess loading, density, N·N and E·N. The acceptance rules subtract the molecule (`mc_utilities.h:322`), the statistics never do. The reported loading is therefore 0.5 to 1 molecule per box and species above the CFCMC estimators (whole molecules, or whole molecules in the λ ≈ 0 / λ ≈ 1 bins). Matters most for small cells and low loadings.

### 8. Single-site charged adsorbate: stale Ewald position on insertion

`Ewald_Energy_Functions.h:480-486` sets `chainsize = Moleculesize − 1` for insertion and `Initialize_Copy_Positions_Together` (`:124-132`) copies only indices `ij < chainsize` into `Old`, assuming the first bead was written by the orientation kernel (`Mol.pos[0] = OldPos`, `mc_widom.h:240`). For `Moleculesize == 1` that kernel is never launched (`mc_swap_utilities.h:28`) and `chainsize = 0`, so `Old.pos[0]`, `Old.charge[0]` and `Old.scaleCoul[0]` hold whatever the previous move left. The Fourier difference of the insertion is then evaluated at a stale position while the stored molecule uses the real trial position and the deletion side reads from `d_a`. The DNN path has the missing `chainsize == 0` branch (`DNN_HostGuest_Energy_Functions.h:37-43`); the Ewald path does not. Triggers: a one-site species with a partial charge under Ewald (ions). Uncharged single sites (CH4, Ar, Xe, Kr) and all multi-site species are unaffected. The end-of-run drift check would expose it.

### 9. Trial-orientation kernel shared-memory initialization

`mc_widom.h:224-240` declares `__shared__ double3 OldPos; __shared__ double Chosenscale, ChosenscaleCoul;` and fills them under `if(i == 0)` with `i` the global thread index. The launch (`mc_widom.h:559`) uses `NumberWidomTrialsOrientations × chainsize` threads in blocks of 128 (`DEFAULTTHREAD`). With more than one block, blocks ≥ 1 never execute the initialization, and their trial orientations are built from uninitialized shared memory on both the insertion and the retrace side. Triggers: `NumberOfTrialOrientations × (atoms − 1) > 128`, e.g. 10 orientations with 14 or more atoms. The largest shipped molecule, p-xylene with 12 atoms at 10 orientations, needs 110 threads and is safe.

### 10. SingleSwap acceptance omits the tail correction

`mc_single_particle.h:236-237`: `Pacc = preFactor * exp(−β tot.total())` with no tail term anywhere in the single-body path, while the CBMC swap applies `TailCorrectionDifference` (`mc_swap_utilities.h:105-108`). Insertion and deletion are consistent with each other, so this is not a detailed-balance violation but a different Hamiltonian: the sampled ensemble is the truncated potential without tail correction, and the running energy lacks the term (drift-visible). Loading low by roughly exp(−βΔU_tail) per molecule in the swap acceptance. Triggers: `TURN_OFF_CBMC_SWAP yes` with tail corrections.

### 11. TMMC with CFCMC or identity swaps

`data_struct.h:209-224` (`UpdateCBCF`) adds `Pacc` to `CMatrix[bin].z` and `1 − Pacc` to `.y` regardless of the sign of `change`, so downward λ moves are credited to the upward column and the downward column stays empty; `AdjustTMBias` then builds ln Π from a one-sided matrix. `ApplyTMMCBias_UpdateCMatrix` (`data_struct.h:1227-1247`) subtracts the fractional molecule only in the swap branch, so "stay" counts from translations land in the wrong macrostate row. `ApplyTMBiasCBCF` (`data_struct.h:308`) is gated on `UseWLBias` instead of `UseTMBias`, so the flattening bias is never applied to CBCF moves. The bin index after a CBCF insertion is off by one across the N boundary. Separately, `IdentitySwapMove` contains no TMMC call and `TMMC::Update` has no `IDENTITY_SWAP` case, so with `RunTMMC yes` an identity swap changes N of two species outside the collection matrix and without the bias. Triggers: `RunTMMC yes` with `CBCFProbability` or `IdentityChangeProbability`.

### 12. CFCMC λ-bin redraw at N_total = 1

`mc_cbcfc.h:259`: `while(SelectednewBin == oldBin || (SelectednewBin < 0 && NumberOfMolecule_for_Component == 1))` redraws the bin change whenever the only molecule is the fractional one and the proposal would delete it. Redrawing on `Δ = 0` is state-independent and harmless; redrawing on `Δ < 0` conditioned on N = 1 is the same construction as the identity-swap bug: the surviving proposals from (N = 1, bin b) are rescaled by 1/(1 − P_del(b)) while their reverses are not. For the TMMC variant (±1 steps) the proposal from (1, bin 0) is +1 with probability 1 against 1/2 for its reverse. Over-populates N ≥ 2 relative to N = 1 in the Henry regime. Moot until #1 is fixed; treat the deletion proposal at N_total = 1 as a rejected attempt instead of redrawing.

### 13. Tail constants for overridden pairs (audit pass)

`main.cpp:99-102` runs `ForceField_Processing` (which fills `Mix_Epsilon/Mix_Sigma` by Lorentz–Berthelot and, if the general tail rule is `yes`, computes `Mix_Tail` from them, `read_data.cpp:1242`), then `OverWrite_Mixing_Rule`, which rewrites ε, σ and the shift for overridden pairs but not the tail, and `OverWriteTailCorrection`, which recomputes tails only for pairs listed in the tail block of `force_field.def`. A pair that is overridden and not re-listed uses the override inside the cutoff and the Lorentz–Berthelot tail outside. Running sum and final recomputation share the same constants, so the drift check cannot see it. The shipped examples use the general rule `no` with explicit per-pair tails and are consistent.

### 14. Gibbs transfer and stale per-box temporaries (audit pass)

`Insertion_Body` reads the tail sign from `TempVal.MoveType` (`mc_swap_utilities.h:106`) and the trial scaling from `TempVal.Scale`; `GibbsParticleXferMove` sets component and molecule for both boxes but never `MoveType` or `Scale`, and `Select_Box_Component_Molecule` resets `TempVal` only for the box running `RunMoves`. If the other box's last move was a CBMC deletion, the transfer's insertion half applies the tail difference with the deletion sign; if it was a CBCF attempt, the molecule is grown with a fractional scale. Pure NVT-Gibbs with translation, rotation, reinsertion and transfer is unaffected.

### 15. Step sizes are adapted during production

`axpy.cu:402-413`, inside `GatherStatisticsDuringSimulation`, which runs every cycle in every phase:
```
if(i%500==0)
{
  ... Update_Max_Translation(...); Update_Max_Rotation(...); Update_Max_SpecialRotation(...);
  if(SystemComponents.VolumeMoveAttempts > 0) Update_Max_VolumeChange(SystemComponents);
}
```
The neighbouring blocks (`axpy.cu:388`, `:417`) are gated on `SimulationMode == EQUILIBRATION`; this one is not, and `Update_Max_GibbsVolume` at `axpy.cu:620-622` is likewise ungated. The rules (`mc_utilities.h:612-711`, `mc_box.h:539-557`) multiply the step by 1.05 or 0.95 (translation, rotation) or by a factor between 0.5 and 1.5 (volume) from the acceptance over the last window, with clamps of 5 Å and 3.14 rad. A proposal whose width depends on the acceptance history makes the production chain non-stationary; RASPA2 freezes the steps for production. For the ±5 % translation and rotation updates the expected effect on averages is small; the volume update is coarser and its acceptance correlates with density. A practical side effect is that `MaxRotation` of small molecules in dilute systems climbs to the 3.14 rad clamp. Side bugs in the same code: `Update_Max_GibbsVolume` (`mc_box.h:546-547`) tests `ratio < 0.5` where it means the scaled factor, so the Gibbs volume step is halved whenever the acceptance is below 0.5; restart-file step sizes are overwritten by `Prepare_Widom` (`fxn_main.h:220-228`). Fix: gate the block on `SimulationMode != PRODUCTION`.

### 16. Volume moves with framework atoms (audit pass)

`mc_box.h:204` and `:388` leave the framework unscaled (`ScaleFirstComponentFramework = false`), but the new-state energy is evaluated with `UseOffset = true` for every component (`VDW_Coulomb.cu:1357-1367`, `Ewald_Energy_Functions.h:814-818`), which reads positions from the second half of each component's allocation. For the framework that half was never written (and `Allocate_size == size`, `read_data.cpp:1788`, so part of it is out of bounds). No input check refuses `NPTVolumeChangeProbability` or `GibbsVolumeChangeProbability` with a non-empty framework. For the supported bulk-fluid use (empty framework CIF) there are no host–guest pairs and the defect is inert.

### 17. Fractional molecule in the volume-move exponent (audit pass)

`Get_TotalNumberOfMolecule_In_Box` (`mc_box.h:183-195`) subtracts the fractional molecule, and the acceptance uses `exp((NMol+1) ln(V'/V) − β(ΔU + PΔV))` (`mc_box.h:284`, `:468`), while `ScalePositions` scales every adsorbate molecule including the fractional one. Each scaled molecule contributes one factor of V to the Jacobian, so the exponent should count it. Requires CFCMC together with volume moves.

### 18. Rotation proposal

`mc_utilities.h:517-525` builds the trial orientation as Rz(c)·Ry(b)·Rx(a) about the first atom with a, b, c independent and uniform in [−m, m]. For a Metropolis move without a Jacobian the density of the rotation Q must equal that of Q⁻¹ with respect to Haar measure. In ZYX Euler angles the Haar density is cos(b) da db dc, and Q⁻¹ = Rx(−a)Ry(−b)Rz(−c) has a different pitch when re-expressed in ZYX order, so the proposal density ratio is cos(b)/cos(b′) ≠ 1, and for a fraction of the cube (about 13 % at m = 30°) the reverse move lies outside the angle cube altogether. A single rotation about a random axis by a symmetric angle, or RASPA2's one-axis-per-move rule, is exactly symmetric.

Size: Metropolis test on a linear rotor in an orienting field, βU = −2 cos θ, exact ⟨cos θ⟩ = 0.53731 (`rotation_proposal_check.py`, 2×10⁴ chains × 2500 recorded steps per entry):

| proposal | m [rad] | ⟨cos θ⟩ | deviation |
|---|---|---|---|
| gRASPA, three axes | 0.524 | 0.53780 | +0.0005 |
| gRASPA, three axes | 1.0 | 0.53785 | +0.0005 |
| gRASPA, three axes | 1.5 | 0.53745 | +0.0001 |
| gRASPA, three axes | 3.1416 | 0.53734 | 0.0000 |
| control, random axis | 1.0 | 0.53733 | 0.0000 |
| control, random axis | 3.1416 | 0.53725 | −0.0001 |

The deviations are at most 1×10⁻³ relative, at the edge of the test's resolution. For non-linear rigid molecules in anisotropic pores the effect was not measured; given these numbers it is unlikely to matter at the percent level, but the fix is a few lines.

## Lower-impact items

- **Printed "PRODUCTION PHASE AVERAGE ENERGY"** (`fxn_main.h:472`) adds `Initial_Energy` (before molecule creation) to the mean of `deltaE` (relative to the state after creation), so the printed total omits the creation energy whenever molecules are created or CFCMC is on. Reporting only. (read-verified)
- **PR-EOS root handling** (`equations_of_state.h:241-261`) compares `Compressibility[0] < Compressibility[1]` without checking that a second root exists; with one real root (supercritical gases, the normal GCMC case) it reads past the vector. Undefined behavior that usually leaves the value intact; guard on `size()`. (read-verified)
- **`NumberOfBlocks` is parsed and never stored** (`read_data.cpp:312-316`); `Nblock` is always 5. Error bars only. The printed `ErrorBar` is 2× the population standard deviation of the five block means, about 2.2× wider than a 2σ standard error. (read-verified)
- **Last-block divisor** uses `Cycles % Blocksize` instead of `Cycles % Nblock` (`print_statistics.cuh:466, 493, 447`); wrong only when production has fewer than about 25 cycles. (audit pass)
- **Widom "Converted to Fugacity"** uses a hard-coded component index 1 and the final instantaneous N (`print_statistics.cuh:108-109`); the Henry coefficient and μ_ex lines are per component and fine. (audit pass)
- **Widom block average** divides by `Nblock` even when a block had no Widom attempts (`print_statistics.cuh:131-141, 177`). (audit pass)
- **Multi-box `CreateNumberOfMolecules`** index repair at `read_data.cpp:2690-2691` only fixes component 1, so Gibbs mixtures start with the wrong per-box counts; for NVT/NPT Gibbs the total N is an input, so coexistence shifts. (audit pass)
- **Restart with the fractional molecule as molecule 0** feeds a λ-scaled `scaleCoul` into the Ewald exclusion constants that every later insertion reuses (`ewald_preparation.h:267, 293, 352-367`). Drift-visible. (audit pass)
- **No check that the cutoff is below half the perpendicular cell widths**; `PBC` (`maths.cuh:427-450`) uses the single nearest fractional image. Consistent in both sums, so drift-blind; RASPA2 aborts in this case. (audit pass)
- **NPT or Gibbs volume moves with a charged framework** add the framework's full Fourier energy to the acceptance (`mc_box.h:264/448`); irrelevant for bulk boxes. (audit pass)
- **Orientation RNG debug write**: the first random triple of a run is overwritten with `{2.3, 4.5, 6.7}` (`data_struct.h:1280-1284`), one bad trial per run. (audit pass)
- **Component selection can lock onto component 0** (`axpy.cu:92-95`): `TempVal.Initialize()` sets `component = 0` and the redraw loop only runs while that component has no moves, so a framework component 0 given move probabilities (separated framework components) would be selected every step. (audit pass)
- **`SINGLE_DELETION` falls through into `SPECIAL_ROTATION`** in `get_new_position` (`mc_utilities.h:562-570`); harmless because deletion does not use the generated positions. (audit pass)
- **A 0.1 Å hard core** in the energy kernels that is not in the Hamiltonian (`VDW_Coulomb.cu:769` and others). Negligible. (audit pass)

## Checked and found correct

- Grand-canonical acceptance prefactors βyPφV/(N+1) and N/(βyPφV); β, pressure and volume units; mole-fraction normalization; component indexing of fugacity coefficient, mole fraction and ideal Rosenbluth weight.
- PR-EOS: mixture a and b with k_ij = 0 at the total pressure, φ_i expression, vapour/liquid root choice by Gibbs energy (apart from the index issue above).
- CBMC: retrace includes the old configuration as trial 0 with k−1 fresh trials for both first bead and orientations; weights normalized identically on both sides; trial selection exactly Boltzmann; overlaps zero-weight on both sides; orientations from Shoemake quaternions (uniform on SO(3)); first-bead positions uniform over triclinic cells (lower-triangular cell matrix plus fractional PBC); reinsertion uses a valid multiple-try scheme; Ewald and tail applied exactly once; Widom averages over all attempts including zero weights.
- Energy functions: one shared LJ kernel with consistent shift and soft-core form; tail formula and its per-move differences; Ewald real, Fourier, self and exclusion terms with one α and one k-set; per-move Fourier differences and their commit on acceptance; all accumulations in double; NIST SPC/E totals reproduced (those references cover only total-energy routines, not any move).
- Statistics accumulate only in production and are freshly sized at its start; Qst matrix algebra and unit conversion; mg/g, mol/kg, g/L conversions.
- Translation: one uniform displacement vector in [−Max, Max]³ applied to all atoms, symmetric, no wrapping; rotation about the first atom with an orthonormal Rodrigues matrix; special rotation about a fixed molecular axis with a symmetric angle; single-body acceptance uses the same energy that is committed.
- NPT volume move: proposal in ln V with the (N+1) ln(V′/V) exponent, molecules scaled rigidly by one reference point, cell, reciprocal vectors and tail correction updated before the energy; Gibbs volume move: uniform in ln(V₁/V₂) with the two-box exponent and exact conservation of the total volume; box choice 50/50 and independent of N.
- Move ladder and component selection (uniform over components with moves, independent of N); 50/50 insertion/deletion split; null deletions at N = 0 counted as stays in TMMC; TMMC collection matrix filled with the unbiased acceptance for ordinary swaps.

## Recommended regression tests

Each of the top findings has an exactly solvable control, in the style of `Examples/IdentitySwap-IdealGas` (U = 0, empty box, exact Poisson loadings):

1. CFCMC: one species, `CBCFProbability` only, βfV = 3; <N_whole> must equal 3 and the CBCF insertion/deletion counters must be non-zero.
2. Block pockets: ideal gas with one blocking sphere of volume fraction f_b; <N> must equal βfV(1 − f_b) and the Widom average (1 − f_b).
3. Cycle definition: ideal gas at βfV = 30 with `UseMaxStep no`; <N> must equal 30 within error (currently about 29.1).
4. Gibbs mixtures: two ideal-gas species in two boxes; compositions must be those of the independent Poisson distributions.
5. Qst: restart an ideal-gas run (or set `CreateNumberOfMolecules`); Qst must equal −RT (the ideal-gas value) instead of −RT + E0/Var(N).
6. A charged single-site species under Ewald; the drift check must stay at zero.
7. Step adaptation: a production run must print identical `MaxTranslation`/`MaxRotation` at its start and end once the gate is added.
8. Orientations: a rigid linear or non-linear molecule in a fixed orienting field (no framework, one tethered molecule) against the quadrature answer, with `MaxRotation` at 30° and at 3.14 rad.

## References

- N. Metropolis et al., J. Chem. Phys. 21, 1087 (1953); W. K. Hastings, Biometrika 57, 97 (1970).
- D. Frenkel and B. Smit, Understanding Molecular Simulation, 2nd ed. (2002), Secs. 3.1, 5.6, Ch. 13.
- V. I. Manousiouthakis and M. W. Deem, J. Chem. Phys. 110, 2753 (1999): balance versus detailed balance.
- W. Shi and E. J. Maginn, J. Chem. Theory Comput. 3, 1451 (2007): continuous fractional component MC.
- A. Poursaeidesfahani, A. Torres-Knoop, D. Dubbeldam, T. J. H. Vlugt, J. Chem. Theory Comput. 12, 1481 (2016): CFCMC estimators.
- G. C. A. M. Mooij, D. Frenkel, Mol. Phys. 74, 41 (1991); J. I. Siepmann, D. Frenkel, Mol. Phys. 75, 59 (1992): configurational-bias MC.
- A. Z. Panagiotopoulos, Mol. Phys. 61, 813 (1987): Gibbs ensemble and the n_B V_A / ((n_A + 1) V_B) transfer rule.
- D. Dubbeldam, S. Calero, D. E. Ellis, R. Q. Snurr, Mol. Simul. 42, 81 (2016): RASPA2.
- Z. Li, K. Shi, D. Dubbeldam, M. Dewing, C. Knight, Á. Vázquez-Mayagoitia, R. Q. Snurr, J. Chem. Theory Comput. 20, 10649 (2024): gRASPA.
