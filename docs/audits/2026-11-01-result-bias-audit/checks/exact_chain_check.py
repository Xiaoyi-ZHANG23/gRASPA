"""Independent exact Markov chain for a U = 0 two-component GCMC with identity swaps.

State (N0, N1), truncated at NMAX per component.  Reservoir activities
B_c = beta f_c V, so the exact grand-canonical answer is Poisson, <N_c> = B_c.

Moves per attempt (weights normalised):
  insertion  : pick component uniformly, accept min(1, B_c/(N_c+1))
  deletion   : pick component uniformly, null if N_c = 0, else min(1, N_c/B_c)
  identity   : pick (OLD, NEW) by RULE, pick a molecule of OLD, accept
               min(1, (N_OLD/(N_NEW+1)) * B_NEW/B_OLD)   (U = 0, Rosenbluth ratio 1)

RULES for drawing (OLD, NEW), n = 2 components, m = number occupied:
  uniform  : OLD, NEW each uniform over n (drawn once); OLD empty -> null attempt;
             OLD == NEW -> self transition.                          (symmetric)
  graspa   : redraw until OLD occupied, i.e. (OLD, NEW) uniform over the m*n pairs
             with OLD occupied; OLD == NEW -> self transition.       (gRASPA d825892)
  redraw   : OLD uniform over the m occupied components, NEW uniform over the n-1 others.
  fixed    : OLD uniform over n, NEW uniform over the n-1 others, OLD empty -> null.
             (proposed fix (a)+(c))
  corrected: graspa proposal with Pacc multiplied by m_x/m_y.        (proposed fix (b))
"""
import numpy as np

B = np.array([0.3, 3.0])
W = {"ins": 0.5, "del": 0.5, "id": 4.0}
NMAX = 40
n = 2

def idx(i, j):
    return i * (NMAX + 1) + j

def pair_probs(rule, N):
    """Return dict {(old,new): prob} for the given rule in state N; sum may be < 1,
    the remainder is a null attempt (no state change)."""
    occ = [c for c in range(n) if N[c] > 0]
    m = len(occ)
    P = {}
    if rule == "uniform":
        for o in range(n):
            for w in range(n):
                P[(o, w)] = 1.0 / n**2          # if N[o]==0 this is a null attempt
    elif rule in ("graspa", "corrected"):
        if m == 0:
            return P
        for o in occ:
            for w in range(n):
                P[(o, w)] = 1.0 / (m * n)
    elif rule == "redraw":
        if m == 0:
            return P
        for o in occ:
            for w in range(n):
                if w != o:
                    P[(o, w)] = 1.0 / (m * (n - 1))
    elif rule == "fixed":
        for o in range(n):
            for w in range(n):
                if w != o:
                    P[(o, w)] = 1.0 / (n * (n - 1))
    else:
        raise ValueError(rule)
    return P

def build(rule):
    S = (NMAX + 1) ** 2
    T = np.zeros((S, S))
    tot = sum(W.values())
    p_ins, p_del, p_id = W["ins"] / tot, W["del"] / tot, W["id"] / tot
    for i in range(NMAX + 1):
        for j in range(NMAX + 1):
            N = [i, j]
            s = idx(i, j)
            # insertion
            for c in range(n):
                if N[c] < NMAX:
                    acc = min(1.0, B[c] / (N[c] + 1))
                    Nn = N.copy(); Nn[c] += 1
                    T[s, idx(*Nn)] += p_ins / n * acc
            # deletion
            for c in range(n):
                if N[c] > 0:
                    acc = min(1.0, N[c] / B[c])
                    Nn = N.copy(); Nn[c] -= 1
                    T[s, idx(*Nn)] += p_del / n * acc
            # identity swap
            m_x = sum(1 for c in range(n) if N[c] > 0)
            for (o, w), pp in pair_probs(rule, N).items():
                if N[o] == 0 or o == w:
                    continue                      # null attempt / self transition
                if N[w] >= NMAX:
                    continue
                Nn = N.copy(); Nn[o] -= 1; Nn[w] += 1
                ratio = (N[o] / (N[w] + 1)) * B[w] / B[o]
                if rule == "corrected":
                    m_y = sum(1 for c in range(n) if Nn[c] > 0)
                    ratio *= m_x / m_y
                T[s, idx(*Nn)] += p_id * pp * min(1.0, ratio)
            # stay
            T[s, s] += 1.0 - T[s].sum()
    return T

def stationary(T):
    S = T.shape[0]
    A = T.T - np.eye(S)
    A[-1, :] = 1.0
    b = np.zeros(S); b[-1] = 1.0
    pi = np.linalg.solve(A, b)
    assert np.all(pi > -1e-12), pi.min()
    assert np.allclose(pi @ T, pi, atol=1e-12)
    return pi

def moments(pi):
    pi2 = pi.reshape(NMAX + 1, NMAX + 1)
    N0 = np.arange(NMAX + 1)
    return (pi2.sum(1) @ N0, pi2.sum(0) @ N0, pi2[0, :].sum())

print(f"exact Poisson: <N0> = {B[0]:.5f}  <N1> = {B[1]:.5f}  P(N0=0) = {np.exp(-B[0]):.5f}")
print(f"{'rule':10s} {'<N0>':>9s} {'<N1>':>9s} {'P(N0=0)':>9s} {'dN0/N0':>8s} {'dN1/N1':>8s}")
for rule in ["uniform", "graspa", "redraw", "fixed", "corrected"]:
    m0, m1, p0 = moments(stationary(build(rule)))
    print(f"{rule:10s} {m0:9.5f} {m1:9.5f} {p0:9.5f} {100*(m0/B[0]-1):+7.1f}% {100*(m1/B[1]-1):+7.1f}%")

# Detailed-balance check of the identity swap alone for the uniform/fixed rules:
# pi(x) P_id(x->y) == pi(y) P_id(y->x) with pi Poisson.
def db_violation(rule):
    worst = 0.0
    tot = sum(W.values())
    for i in range(1, NMAX):
        for j in range(0, NMAX - 1):
            x = [i, j]; y = [i - 1, j + 1]
            def flow(a, b_, o, w):
                pp = pair_probs(rule, a).get((o, w), 0.0)
                if a[o] == 0: return 0.0
                r = (a[o] / (a[w] + 1)) * B[w] / B[o]
                if rule == "corrected":
                    r *= sum(1 for c in range(n) if a[c] > 0) / sum(1 for c in range(n) if b_[c] > 0)
                return pp * min(1.0, r)
            pix = np.exp(-B.sum()) * B[0]**i * B[1]**j / (np.math.factorial(i) * np.math.factorial(j))
            piy = np.exp(-B.sum()) * B[0]**(i-1) * B[1]**(j+1) / (np.math.factorial(i-1) * np.math.factorial(j+1))
            f = pix * flow(x, y, 0, 1); g = piy * flow(y, x, 1, 0)
            if f + g > 0:
                worst = max(worst, abs(f - g) / (f + g))
    return worst
import math
np.math = math
for rule in ["uniform", "graspa", "redraw", "fixed", "corrected"]:
    print(f"max relative detailed-balance violation of the identity swap, rule {rule:10s}: {db_violation(rule):.3g}")
