# IdentitySwap-IdealGas: exact check of the identity-swap move

A two-component ideal gas in an empty 30 Å cubic box, sampled with insertions,
deletions and identity swaps. Every interaction is zero (both Lennard-Jones
epsilons are 0, no charges, no framework atoms), so the grand-canonical
answer is known exactly and any bias in a move shows up directly in the
loadings.

## Exact answer

For non-interacting molecules the configurational integral of N molecules in
a volume V is V^N, so the mu-V-T weight of a state with N_c molecules of each
species c is

    pi(N_A, N_B) = prod_c  B_c^{N_c} / N_c!   x const,      B_c = beta f_c V,

a product of Poisson distributions with means B_c (f_c is the fugacity,
here y_c P because the fugacity coefficient is 1). The input is chosen so that

    beta P V = 3.0  (P = 460216.3333 Pa, V = 27000 Å^3, T = 300 K, k_B = 1.380649e-23 J/K)
    y_A = 0.1, y_B = 0.9       ->      <N_A> = B_A = 0.3,   <N_B> = B_B = 2.7

(`MolFraction` values are normalised by gRASPA; 0.1 and 0.9 are already
normalised.) The exact P(N_A = 0) is exp(-0.3) = 0.7408, so species A is
empty most of the time, which is the regime in which the old proposal was
most biased.

## Move ladder

Per component: `SwapProbability 1.0` (insertion and deletion, half each) and
`IdentityChangeProbability 4.0`. A component is drawn uniformly and its
ladder is read, so per MC step: insertion 0.1, deletion 0.1, identity swap
0.8. `MaxStepPerCycle 1` makes one cycle one MC step; the loadings are
averaged over 2,000,000 production steps in 5 blocks.

The identity swap turns a molecule of OLD into a molecule of NEW at the same
place. Its acceptance rule in `src_clean/mc_swap_moves.h` is

    P_acc = min(1, [N_OLD / (N_NEW + 1)] x [f_NEW / f_OLD] x W_NEW / W_OLD),

the identity-change rule of RASPA2 and of Martin and Siepmann (W are
Rosenbluth weights, equal to 1 here). This rule assumes that the ordered
pair (OLD, NEW) is proposed with the same probability from every state, so
that the pair-selection probability cancels between a swap and its reverse.
The code now draws OLD uniformly over the adsorbate species and NEW
uniformly over the others, once, and counts a draw with an empty OLD as a
rejected attempt: every ordered pair has probability 1/(n(n-1)).

Before this change the pair was redrawn whenever OLD was empty. The pair
probability was then 1/(n m), with m the number of occupied species, so a
swap that fills an empty species was proposed twice as often (for n = 2) as
its reverse while the acceptance rule assumed equal proposal probabilities.
Detailed balance, pi(x) alpha(x->y) acc(x->y) = pi(y) alpha(y->x) acc(y->x),
then fails by a factor of 2 whenever a swap fills or empties a species.

## Expected result

`exact_chain.py` builds the exact Markov chain of this ladder for several
pair-selection rules and solves for its stationary distribution
(`python exact_chain.py`, needs numpy):

| pair-selection rule                                  | <N_A>   | <N_B>   | P(N_A = 0) |
|------------------------------------------------------|---------|---------|------------|
| exact (Poisson)                                      | 0.30000 | 2.70000 | 0.74082    |
| draw once, NEW != OLD, empty OLD rejected (current)  | 0.30000 | 2.70000 | 0.74082    |
| redraw while OLD is empty (gRASPA up to d825892)     | 0.38855 | 2.36899 | 0.65889    |

The old rule over-populates the minor species by about 30% and
under-populates the major one by about 12%. In `output.txt`, under
`BLOCK AVERAGES (LOADING: # MOLECULES)`, the `Overall: Average` lines of
`COMPONENT [1] (A)` and `COMPONENT [2] (B)` should therefore read about
0.300 and 2.700. With 2,000,000 steps the standard error of the two means is
of the order of 0.002 and 0.01 (rough estimates from the Poisson variances
and a correlation time of tens of steps), far below the bias of the old
rule (0.09 and 0.33).

## Running the check

    cd Examples/IdentitySwap-IdealGas
    ../../src_clean/nvc_main.x > output.txt      # or the binary you built
    cd ..
    pytest -s test_identity_swap_ideal_gas.py

The test reads the last loading block of `output.txt` and requires
|<N_c> - B_c| <= max(ErrorBar_c, 0.03 B_c); it is skipped when `output.txt`
is absent. (gRASPA's `ErrorBar` is twice the standard deviation of the
block averages.)

This example was added together with the fix of the pair selection and has
not yet been run on a GPU in this repository: the numbers above come from the
exact chain, not from a gRASPA run. The first run of the example should be
checked against them.

## References

- N. Metropolis, A. W. Rosenbluth, M. N. Rosenbluth, A. H. Teller, E. Teller,
  J. Chem. Phys. 21, 1087 (1953); W. K. Hastings, Biometrika 57, 97 (1970):
  the acceptance rule for an asymmetric proposal carries the ratio
  alpha(y->x)/alpha(x->y).
- D. Frenkel and B. Smit, Understanding Molecular Simulation, 2nd ed.,
  Academic Press (2002), Section 5.6 (grand-canonical MC) and Chapter 13
  (configurational-bias MC and the Rosenbluth weight).
- M. G. Martin and J. I. Siepmann, J. Am. Chem. Soc. 119, 8921 (1997),
  doi:10.1021/ja964218q: configurational-bias identity-exchange moves in
  multicomponent mixtures.
- D. Dubbeldam, S. Calero, D. E. Ellis, R. Q. Snurr, Mol. Simul. 42, 81
  (2016): RASPA2, whose identity-change move gRASPA follows.
- Z. Li, K. Shi, et al., J. Chem. Theory Comput. 20, 10649 (2024),
  doi:10.1021/acs.jctc.4c01058: gRASPA.
