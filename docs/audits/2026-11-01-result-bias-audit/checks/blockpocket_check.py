"""Exact chain for gRASPA's block-pocket treatment in an ideal gas.

Setting: U = 0, activity a = beta f V, a fraction f_b of the box is a blocked
pocket, k first-bead trials.  Exact answer (configurations in pockets have
zero weight): <N> = a (1 - f_b).

gRASPA insertion (mc_widom.h:445-497 + mc_swap_utilities.h): trial 0 blocked
-> whole attempt rejected (prob f_b); otherwise blocked trials get weight 0,
W_ins = (k - B)/k with B ~ Binomial(k-1, f_b); acc = min(1, a/(N+1) W_ins).
gRASPA deletion (Deletion_Body): no pocket logic, W_del = 1 at U = 0;
acc = min(1, N/a).
A symmetric treatment (RASPA2-style: no zeroing on either side, reject only
a selected configuration inside a pocket) gives the exact answer.
"""
import numpy as np
from math import comb

def stationary_mean(a, f_b, k, NMAX=None):
    if NMAX is None: NMAX = int(a + 12 * a**0.5 + 30)
    # E over B of min(1, c (k-B)/k)
    pB = [comb(k - 1, b) * f_b**b * (1 - f_b)**(k - 1 - b) for b in range(k)]
    def acc_ins(N):
        c = a / (N + 1)
        return (1 - f_b) * sum(pB[b] * min(1.0, c * (k - b) / k) for b in range(k))
    def acc_ins_sym(N):           # RASPA2-style: W_ins = 1, final check rejects P(selected in pocket) = f_b
        return (1 - f_b) * min(1.0, a / (N + 1))
    def solve(acc):
        T = np.zeros((NMAX + 1, NMAX + 1))
        for N in range(NMAX + 1):
            pi_ = 0.5 * acc(N) if N < NMAX else 0.0
            pd_ = 0.5 * min(1.0, N / a) if N > 0 else 0.0
            if N < NMAX: T[N, N + 1] = pi_
            if N > 0: T[N, N - 1] = pd_
            T[N, N] = 1 - pi_ - pd_
        A = T.T - np.eye(NMAX + 1); A[-1] = 1; b = np.zeros(NMAX + 1); b[-1] = 1
        pi = np.linalg.solve(A, b)
        return pi @ np.arange(NMAX + 1)
    return solve(acc_ins), solve(acc_ins_sym)

print(f"{'a':>5} {'f_b':>5} {'k':>3} {'exact a(1-f_b)':>15} {'gRASPA':>9} {'rel':>8} {'symmetric':>10}")
for a in (5.0, 30.0):
    for f_b in (0.05, 0.10, 0.20):
        for k in (8, 10):
            g, s = stationary_mean(a, f_b, k)
            ex = a * (1 - f_b)
            print(f"{a:5.0f} {f_b:5.2f} {k:3d} {ex:15.4f} {g:9.4f} {100*(g/ex-1):+7.1f}% {s:10.4f}")
