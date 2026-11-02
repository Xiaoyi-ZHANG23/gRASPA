"""Direct simulation cross-check of cycle_length_check.py (LAM=30, p_swap=1/3)."""
import random, math
random.seed(7)
LAM, p_swap = 30.0, 1/3
def run(cycles, fixed=None):
    N = 30; s = 0.0; s2 = 0.0; cnt = 0
    for c in range(cycles):
        S = max(20, N) if fixed is None else fixed
        for _ in range(S):
            u = random.random()
            if u < p_swap / 2:
                if random.random() < min(1.0, LAM / (N + 1)): N += 1
            elif u < p_swap:
                if N > 0 and random.random() < min(1.0, N / LAM): N -= 1
        s += N; s2 += N * N; cnt += 1
    m = s / cnt
    return m, math.sqrt((s2 / cnt - m * m) / cnt) * 3  # ~3 x naive SE (samples correlated)
for label, fixed in (("max(20,N) cycles", None), ("fixed 20-step cycles", 20)):
    m, se = run(400000, fixed)
    print(f"{label:22s} <N> = {m:.3f} (exact 30, bias {m-30:+.3f}; naive uncertainty ~{se:.3f})")
