"""Does sampling once per 'cycle' of max(20, N) steps bias <N>?

Ideal-gas GCMC on N with activity LAM (exact answer: Poisson, <N> = LAM).
Per MC step: with probability p_swap attempt insertion (1/2) or deletion (1/2);
otherwise a move that leaves N unchanged (translation).  One-step kernel K is
exactly pi-invariant.  gRASPA/RASPA record observables once per cycle, and a
cycle is S(N) = max(20, N) steps with N read at the start of the cycle.  The
recorded chain has kernel Q(x, .) = K^{S(x)}(x, .); its stationary law is not
pi when S varies.  We solve it exactly on 0..NMAX.
"""
import numpy as np, math

def run(LAM, p_swap, NMAX=None, fixed_steps=None):
    if NMAX is None:
        NMAX = int(LAM + 12 * math.sqrt(LAM) + 30)
    n = np.arange(NMAX + 1)
    K = np.zeros((NMAX + 1, NMAX + 1))
    for N in range(NMAX + 1):
        p_ins = p_swap / 2 * min(1.0, LAM / (N + 1)) if N < NMAX else 0.0
        p_del = p_swap / 2 * min(1.0, N / LAM) if N > 0 else 0.0
        if N < NMAX: K[N, N + 1] = p_ins
        if N > 0:    K[N, N - 1] = p_del
        K[N, N] = 1.0 - p_ins - p_del
    S = np.maximum(20, n) if fixed_steps is None else np.full(NMAX + 1, fixed_steps)
    # Q[x, :] = (K^{S(x)})[x, :]
    powers = {}
    P = np.eye(NMAX + 1)
    for s in range(1, int(S.max()) + 1):
        P = P @ K
        powers[s] = P
    Q = np.vstack([powers[int(S[x])][x] for x in range(NMAX + 1)])
    A = Q.T - np.eye(NMAX + 1); A[-1, :] = 1.0
    b = np.zeros(NMAX + 1); b[-1] = 1.0
    pi_Y = np.linalg.solve(A, b)
    assert pi_Y.min() > -1e-12
    meanN = pi_Y @ n
    var_exact = LAM
    return meanN, pi_Y

print(f"{'LAM':>5} {'p_swap':>7} {'<N> recorded':>13} {'bias':>8} {'-Var/<N>':>9} {'fixed S=20':>11}")
for LAM in (10.0, 30.0, 100.0):
    for p_swap in (1.0, 1/3, 0.1):
        m, _ = run(LAM, p_swap)
        mfix, _ = run(LAM, p_swap, fixed_steps=20)
        print(f"{LAM:5.0f} {p_swap:7.3f} {m:13.4f} {m-LAM:+8.4f} {-LAM/LAM:9.3f} {mfix-LAM:+11.1e}")
