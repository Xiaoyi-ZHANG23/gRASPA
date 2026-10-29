"""Exact stationary loadings of the IdentitySwap-IdealGas example for several
identity-swap proposal rules.

The example is a two-component ideal gas (U = 0) in the grand-canonical
ensemble, so the exact distribution of (N_A, N_B) is a product of Poisson
distributions with means B_c = beta * f_c * V (the configurational integral of
N non-interacting molecules in V is V^N, and the mu-V-T weight
(beta f V)^N / N! sums to exp(beta f V)).  Hence <N_c> = B_c, independent of
any move that satisfies detailed balance.

The script builds the exact Markov chain of the example's move ladder on the
truncated state space 0 <= N_c <= NMAX and solves pi T = pi.  Per attempt,
mirroring gRASPA's RunMoves for this input (a component is drawn uniformly,
then its ladder is read; both components here have the same ladder):

  insertion   (weight 0.5): component uniform, accept min(1, B_c/(N_c+1))
  deletion    (weight 0.5): component uniform, null if N_c == 0, else min(1, N_c/B_c)
  identity    (weight 4.0): draw (OLD, NEW) by RULE, pick a molecule of OLD uniformly,
                            accept min(1, (N_OLD/(N_NEW+1)) * B_NEW/B_OLD)
                            (the Rosenbluth ratio is exactly 1 at U = 0)

Rules for (OLD, NEW), with n = 2 components and m = number of occupied ones:

  fixed     : OLD uniform over n, NEW uniform over the n-1 others, drawn once;
              OLD empty -> rejected attempt.  This is what mc_swap_moves.h does now.
              Each ordered pair has probability 1/(n(n-1)) from every state.
  graspa    : redraw until OLD is occupied (gRASPA up to d825892): the pair is
              uniform over the m*n pairs with OLD occupied, so its probability
              1/(n m) depends on the state.  OLD == NEW is a self transition.
  uniform   : OLD, NEW uniform over n, drawn once; OLD empty -> rejected; OLD == NEW
              -> self transition.  Symmetric as well (probability 1/n^2).
  redraw    : OLD uniform over the m occupied components, NEW over the n-1 others.
  corrected : the graspa proposal with the acceptance multiplied by m_x/m_y,
              the Metropolis-Hastings correction for its asymmetry.

Only the state-independent rules (fixed, uniform, corrected) reproduce the
Poisson means.  Run:  python exact_chain.py
"""
import math
import numpy as np

B = np.array([0.3, 2.7])             # beta f_c V of components A and B in simulation.input
WEIGHTS = {"ins": 0.5, "del": 0.5, "id": 4.0}
NMAX = 40                            # P(Poisson(2.7) > 40) ~ 1e-27, truncation is harmless
n = 2


def pair_probs(rule, N):
    """{(old, new): probability} for drawing the ordered pair in state N.
    Probabilities may sum to less than one; the rest is a null attempt."""
    occ = [c for c in range(n) if N[c] > 0]
    m = len(occ)
    P = {}
    if rule == "fixed":
        for o in range(n):
            for w in range(n):
                if w != o:
                    P[(o, w)] = 1.0 / (n * (n - 1))
    elif rule in ("graspa", "corrected"):
        for o in occ:
            for w in range(n):
                P[(o, w)] = 1.0 / (m * n)
    elif rule == "uniform":
        for o in range(n):
            for w in range(n):
                P[(o, w)] = 1.0 / n ** 2
    elif rule == "redraw":
        for o in occ:
            for w in range(n):
                if w != o:
                    P[(o, w)] = 1.0 / (m * (n - 1))
    else:
        raise ValueError(rule)
    return P


def idx(i, j):
    return i * (NMAX + 1) + j


def transition_matrix(rule):
    S = (NMAX + 1) ** 2
    T = np.zeros((S, S))
    tot = sum(WEIGHTS.values())
    p_ins, p_del, p_id = (WEIGHTS[k] / tot for k in ("ins", "del", "id"))
    for i in range(NMAX + 1):
        for j in range(NMAX + 1):
            N = [i, j]
            s = idx(i, j)
            for c in range(n):                                   # insertion
                if N[c] < NMAX:
                    Nn = N.copy(); Nn[c] += 1
                    T[s, idx(*Nn)] += p_ins / n * min(1.0, B[c] / (N[c] + 1))
            for c in range(n):                                   # deletion
                if N[c] > 0:
                    Nn = N.copy(); Nn[c] -= 1
                    T[s, idx(*Nn)] += p_del / n * min(1.0, N[c] / B[c])
            m_x = sum(1 for c in range(n) if N[c] > 0)          # identity swap
            for (o, w), pp in pair_probs(rule, N).items():
                if N[o] == 0 or o == w or N[w] >= NMAX:
                    continue
                Nn = N.copy(); Nn[o] -= 1; Nn[w] += 1
                ratio = (N[o] / (N[w] + 1)) * B[w] / B[o]
                if rule == "corrected":
                    m_y = sum(1 for c in range(n) if Nn[c] > 0)
                    ratio *= m_x / m_y
                T[s, idx(*Nn)] += p_id * pp * min(1.0, ratio)
            T[s, s] += 1.0 - T[s].sum()
    return T


def stationary(T):
    S = T.shape[0]
    A = T.T - np.eye(S)
    A[-1, :] = 1.0                      # replace one redundant equation by normalisation
    b = np.zeros(S); b[-1] = 1.0
    pi = np.linalg.solve(A, b)
    assert pi.min() > -1e-12 and np.allclose(pi @ T, pi, atol=1e-12)
    return pi.reshape(NMAX + 1, NMAX + 1)


def main():
    N = np.arange(NMAX + 1)
    print(f"exact (Poisson): <N_A> = {B[0]:.5f}  <N_B> = {B[1]:.5f}  P(N_A = 0) = {math.exp(-B[0]):.5f}")
    print(f"{'rule':10s} {'<N_A>':>9s} {'<N_B>':>9s} {'P(N_A=0)':>9s} {'dA/A':>8s} {'dB/B':>8s}")
    for rule in ("fixed", "uniform", "corrected", "graspa", "redraw"):
        pi = stationary(transition_matrix(rule))
        mA, mB, p0 = pi.sum(1) @ N, pi.sum(0) @ N, pi[0, :].sum()
        print(f"{rule:10s} {mA:9.5f} {mB:9.5f} {p0:9.5f} {100 * (mA / B[0] - 1):+7.1f}% {100 * (mB / B[1] - 1):+7.1f}%")


if __name__ == "__main__":
    main()
