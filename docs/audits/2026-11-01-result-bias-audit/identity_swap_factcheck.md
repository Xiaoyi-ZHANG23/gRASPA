# Fact check of "IdentitySwapMove uses an asymmetric proposal ..." and what was done

Checked on 29 Oct 2026 against the local clone of Xiaoyi-ZHANG23/gRASPA (main = e6da4e3, which contains d825892 and e4edfc2) and graspa-rewrite (a448f28).

## Verdict

The central claim is correct and I reproduced its numbers independently. The fix recommended in the report, (a) draw once and reject on an empty source combined with (c) no self swaps, is implemented in commit 68bde66 on branch `claude/magical-lamport-tzkjy0` of Xiaoyi-ZHANG23/gRASPA, together with a U = 0 regression example. Nothing was pushed to snurr-group/gRASPA and no GitHub issue was filed.

## Claim-by-claim

| # | Claim in the report | Result |
|---|---|---|
| 1 | Quoted code at `mc_swap_moves.h:222-237` and `373-376` of d825892 | Verbatim match. Identity swap code is identical in e4edfc2, d825892 and current main apart from the `static_cast<int>` on the ExcludeList line. |
| 2 | `GetPrefactor` INSERTION = beta y P phi V/(1+N), DELETION = N/(beta y P phi V), N read from current counts (minus 1 if the species has a fractional molecule) | Correct (`mc_utilities.h:315-349`). |
| 3 | Pair probability 1/(n m) with the redraw; forward/reverse ratio 2 when a swap fills or empties a species; ratio 1 when N_A = 1 and N_B = 0 | Correct. The redraw loop depends only on `NOld`, so the kept pairs are uniform over the m n pairs with an occupied OLD. |
| 4 | Acceptance has no pair-selection factor; MH needs acc ratio = pi(y)/pi(x) x alpha(y->x)/alpha(x->y) = pi(y)/pi(x) x m_x/m_y | Correct. The sign of the bias follows: filling an empty species is over-proposed, so the often-empty species is over-populated. |
| 5 | Exact chain, B = 0.3 and 3.0, weights 0.5/0.5/4: symmetric 0.30000/3.00000, gRASPA rule 0.39252/2.63113, simpler redraw 0.39825/2.60978 | Reproduced to all five digits with an independent implementation (`checks/exact_chain_check.py`, dense 41 x 41 state space, stationary vector by linear solve). The correction rule (b), Pacc x m_x/m_y, also gives 0.30000/3.00000. Maximum relative detailed-balance violation of the identity swap alone: 0.333 for the gRASPA rule (a factor 2), 3e-16 for the symmetric rules. |
| 6 | +31 % / -12 % | 30.8 % and -12.3 %. |
| 7 | XeKr shift about -8 % Xe, +8 % Kr | Consistent with the table (29.1 vs 31.6, 1.22 vs 1.137). The runs themselves could not be repeated here (no GPU, no CUDA compiler). |
| 8 | graspa-rewrite uses the symmetric rule | Confirmed: `src/plugins/moves/identity_swap.cu` draws `from` over n and `to` over the other n-1 before reading any count, null attempt if `from` is empty. |
| 9 | `NumberOfMolecule_for_Component` includes the fractional molecule, so a species with only the fractional molecule counts as occupied | Correct. Also: the selected OLD molecule can itself be the fractional molecule; the deletion move refuses that case, the identity swap does not. Pre-existing, outside this fix (queued as a separate task). |
| 10 | Self swaps: product of the two prefactors is N/(N+1), correct factor 1 | Correct. |
| 11 | Self swaps: "detailed balance holds only for c < R < 1/c" | **Not correct.** With acc_f = min(1, cR), acc_r = min(1, c/R), the ratio is cR for R <= c, R^2 for c <= R <= 1/c, and R/c for R >= 1/c. It equals R only at R = 1. Detailed balance is violated for every R != 1, by a factor between N/(N+1) and (N+1)/N. The conclusion "error of order 1/N" stands; "only in the tails" does not. Suggested wording below. |
| 12 | "Single-adsorbate-component systems also have no pair choice to bias" | True for composition, but with one adsorbate species every identity-swap attempt was a self swap, i.e. a reinsertion with the biased factor of claim 10. Such inputs were sampling configurations with an O(1/N) error. After the fix the move does nothing in that case and gRASPA prints a warning. |
| 13 | Fix (a) correct; (c) proposal probability 1/(n(n-1)) | Correct. Both confirmed by the exact chain ("uniform" and "fixed" rules). |

One addition the report does not mention: `Get_Uniform_Random()` is `rand()/RAND_MAX`, which returns exactly 1.0 with probability 1/(RAND_MAX+1). The old `while` loop's range test also guarded against that. The new code keeps a range redraw that depends on the random number only, so the proposal stays state-independent.

## Suggested replacement for the "Self swaps" paragraph

> `GetPrefactor` reads the current count N in both calls, so for NEW == OLD the product is [beta f V/(1+N)] x [N/(beta f V)] = N/(N+1), whereas a reinsertion needs the factor 1. With c = N/(N+1) and R = W_new/W_old the acceptances are min(1, cR) forward and min(1, c/R) reverse. Their ratio is cR for R <= c, R^2 for c <= R <= 1/c and R/c for R >= 1/c, instead of R. Detailed balance therefore fails for every R != 1, by a factor between N/(N+1) and (N+1)/N. Excluding self swaps (fix (c)) removes the problem; a run with a single adsorbate species and IdentityChangeProbability > 0 was affected as well, because every attempt was a self swap.

## What the commit changes (68bde66)

- `src_clean/mc_swap_moves.h`: OLD uniform over adsorbate species, NEW uniform over the others, drawn once before any count is read; empty OLD = rejected attempt (counted in `IdentitySwap_Total_TO` etc.); no self swaps; single adsorbate species = return; range redraws for the 1.0 case on the three draws. The acceptance code is untouched; its `NEW == OLD` branch is now unreachable and was left in place to keep the diff small.
- `src_clean/fxn_main.h`: warning in `prepare_MixtureStats` when IdentityChangeProbability > 0 with one adsorbate species.
- `Examples/IdentitySwap-IdealGas/`: empty 30 A cubic box, two zero-epsilon species, beta P V = 3 with y = 0.1/0.9, so <N_A> = 0.3 and <N_B> = 2.7 exactly; ladder insertion 0.1, deletion 0.1, identity 0.8 per step; `exact_chain.py` prints old rule 0.38855/2.36899 and new rule 0.30000/2.70000; README with the derivation.
- `Examples/test_identity_swap_ideal_gas.py`: parses the last loading block of `output.txt`, requires |<N_c> - B_c| <= max(ErrorBar_c, 0.03 B_c); skipped until the example has been run; the parser is exercised on the shipped XeKr output in CI.

Not verified here: compilation (no nvcc or nvc++ in this container) and the example run. Please build, run `Examples/IdentitySwap-IdealGas` once, and run `pytest -s` in `Examples/`. The shipped `Examples/XeKr-Mixture/output.txt` predates the fix and should move from Xe 29.2 / Kr 1.23 to about 31.6 / 1.14 when regenerated.

## Statistical-mechanics background, briefly

- Grand-canonical weight of N non-interacting molecules in V: (beta f V)^N / N! (configurational integral V^N, de Broglie factors absorbed in f). Summing over N gives exp(beta f V), so N is Poisson with mean beta f V and the two species are independent. Any move set that satisfies detailed balance w.r.t. this weight, with insertion/deletion present for ergodicity in N, must return these means.
- Metropolis-Hastings: for a proposal alpha(x->y), acc(x->y) = min(1, pi(y) alpha(y->x) / (pi(x) alpha(x->y))). The identity-change rule N_OLD f_NEW/((N_NEW+1) f_OLD) x W_NEW/W_OLD already contains the molecule-selection factors 1/N_OLD and 1/(N_NEW+1) and the factorials of pi; what is missing in the old code is only the pair-selection factor m_x/m_y.
- Why insertion/deletion do not "repair" the bias: each move is reversible w.r.t. its own target. The old identity swap is reversible w.r.t. pi(x) m_x on each shell of fixed total N; insertion/deletion w.r.t. pi. The combined chain's stationary distribution is neither and depends on the move weights, which is why the exact-chain numbers depend on the 0.5/0.5/4 ladder.

## References (verified where marked)

1. N. Metropolis, A. W. Rosenbluth, M. N. Rosenbluth, A. H. Teller, E. Teller, J. Chem. Phys. 21, 1087 (1953).
2. W. K. Hastings, Biometrika 57, 97 (1970). Asymmetric proposals.
3. D. Frenkel, B. Smit, Understanding Molecular Simulation, 2nd ed., Academic Press (2002), Sec. 5.6 (grand-canonical MC), Ch. 13 (CBMC).
4. D. A. Kofke, E. D. Glandt, Mol. Phys. 64, 1105 (1988). Identity changes in the semigrand ensemble. (citation verified via Kofke's publication list)
5. M. G. Martin, J. I. Siepmann, J. Am. Chem. Soc. 119, 8921 (1997), doi:10.1021/ja964218q. CBMC identity exchange in mixtures. (citation verified; the search did not surface the paper's text)
6. D. Dubbeldam, S. Calero, D. E. Ellis, R. Q. Snurr, Mol. Simul. 42, 81 (2016), doi:10.1080/08927022.2015.1010082. RASPA2, whose identity-change rule gRASPA follows. (verified)
7. Z. Li, K. Shi, D. Dubbeldam, M. Dewing, C. Knight, A. Vazquez-Mayagoitia, R. Q. Snurr, J. Chem. Theory Comput. 20, 10649 (2024), doi:10.1021/acs.jctc.4c01058. gRASPA. (verified)
