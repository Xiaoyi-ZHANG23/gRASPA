"""Is gRASPA's rotation proposal symmetric?  Metropolis test on a linear rotor.

gRASPA (mc_utilities.h, get_new_position, case ROTATION) rotates the molecule
by Q = Rz(c) Ry(b) Rx(a) with a, b, c independent and uniform in [-m, m], and
accepts with min(1, exp(-beta dU)) and no Jacobian.  A proposal is symmetric
only if the density of Q equals that of Q^-1 with respect to Haar measure.
Test: linear rotor with U = -2 kT cos(theta); exact <cos theta> = coth(2) - 1/2
= 0.537315.  Many independent chains are advanced in parallel with numpy.
A single-axis rotation (random axis, uniform angle) is exactly symmetric and
serves as the control.
"""
import numpy as np, math
rng = np.random.default_rng(1)
EXACT = 1 / math.tanh(2.0) - 0.5

def rot(axis, ang):
    """Rotation matrices (n,3,3) about unit axes (n,3) by angles (n,)."""
    x, y, z = axis[:, 0], axis[:, 1], axis[:, 2]
    c, s = np.cos(ang), np.sin(ang); C = 1 - c
    R = np.empty((len(ang), 3, 3))
    R[:, 0, 0] = c + x*x*C;   R[:, 0, 1] = x*y*C - z*s; R[:, 0, 2] = x*z*C + y*s
    R[:, 1, 0] = y*x*C + z*s; R[:, 1, 1] = c + y*y*C;   R[:, 1, 2] = y*z*C - x*s
    R[:, 2, 0] = z*x*C - y*s; R[:, 2, 1] = z*y*C + x*s; R[:, 2, 2] = c + z*z*C
    return R

def run(m, nchains=20000, nsteps=3000, burn=500, mode="xyz"):
    # start from random directions on the sphere
    u = rng.normal(size=(nchains, 3)); u /= np.linalg.norm(u, axis=1)[:, None]
    ex = np.tile([1.0, 0, 0], (nchains, 1)); ey = np.tile([0, 1.0, 0], (nchains, 1)); ez = np.tile([0, 0, 1.0], (nchains, 1))
    acc_sum = 0.0; obs = 0.0; nobs = 0
    for step in range(nsteps):
        if mode == "xyz":
            a, b, c = (m * (2 * rng.random(nchains) - 1) for _ in range(3))
            Q = rot(ez, c) @ rot(ey, b) @ rot(ex, a)          # Rz(c) Ry(b) Rx(a), lab frame
        else:                                                 # control: random axis, uniform angle
            ax = rng.normal(size=(nchains, 3)); ax /= np.linalg.norm(ax, axis=1)[:, None]
            Q = rot(ax, m * (2 * rng.random(nchains) - 1))
        un = np.einsum("nij,nj->ni", Q, u)
        dU = -2.0 * (un[:, 2] - u[:, 2])                      # beta U = -2 cos theta
        acc = rng.random(nchains) < np.exp(-dU)
        u[acc] = un[acc]
        if step >= burn:
            obs += u[:, 2].mean(); nobs += 1; acc_sum += acc.mean()
    return obs / nobs, acc_sum / nobs

print(f"exact <cos theta> = {EXACT:.5f}; statistical error of each estimate about 2e-4")
print(f"{'proposal':34s} {'m [rad]':>8s} {'<cos theta>':>12s} {'dev':>9s} {'acc':>6s}")
for m in (0.5236, 1.0, 1.5, 3.1416):
    v, a = run(m, mode="xyz")
    print(f"{'gRASPA Rz(c)Ry(b)Rx(a), cube':34s} {m:8.4f} {v:12.5f} {v-EXACT:+9.5f} {a:6.3f}")
for m in (1.0, 3.1416):
    v, a = run(m, mode="axis")
    print(f"{'control: random axis, one angle':34s} {m:8.4f} {v:12.5f} {v-EXACT:+9.5f} {a:6.3f}")
