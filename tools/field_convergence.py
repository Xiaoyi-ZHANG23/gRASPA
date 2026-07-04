#!/usr/bin/env python3
"""field_convergence.py — decide whether a 3D density FIELD is converged.

Why not scalar uptake: a GCMC run can report a well-converged average loading
while the *spatial* field is still noisy — many voxels are visited too few times.
We therefore compare independent replicate runs of the SAME MOF/state-point
(different RandomSeed) voxel-by-voxel, and require spatial agreement.

Metrics (computed pairwise over all replicates, then aggregated worst-case):

  Continuous Tanimoto (a.k.a. Tanimoto/Jaccard coefficient for real vectors):
        T(a,b) = <a,b> / (<a,a> + <b,b> - <a,b>)
    For non-negative fields T in [0,1]; T=1 iff a==b. This is the DeepAPD-style
    criterion; we PASS a pair when T >= threshold (default 0.75).

  Jensen-Shannon divergence: normalise each field to a probability distribution
    p,q (sum=1, epsilon-smoothed), m=(p+q)/2, then
        JS = 0.5*KL(p||m) + 0.5*KL(q||m),  KL(p||q)=sum p*log(p/q)
    Symmetric, bounded in [0, ln2] nats; we also report it in bits ([0,1]).
    Reported for context; the pass/fail gate is Tanimoto by default.

  (optional) Binary Tanimoto: threshold each field at a fraction of its max
    occupancy (--binary-cutoff) and take Jaccard over the boolean "occupied"
    masks — a stricter geometric agreement of WHERE density sits.

Verdict = PASS iff the minimum pairwise (chosen) metric clears the threshold.

numpy only.

Examples:
  field_convergence.py rep0/.../grid.vtk rep1/.../grid.vtk
  field_convergence.py --run-dirs /data/pilot/rep0 /data/pilot/rep1 --threshold 0.75
  field_convergence.py --selftest
"""
from __future__ import annotations
import argparse
import itertools
import sys
import numpy as np

try:
    from gridio import load_grid, find_grid
except ImportError:
    sys.path.insert(0, __file__.rsplit("/", 1)[0])
    from gridio import load_grid, find_grid


def tanimoto(a, b):
    a = a.ravel(); b = b.ravel()
    ab = float(np.dot(a, b))
    denom = float(np.dot(a, a) + np.dot(b, b) - ab)
    if denom <= 0:
        return 1.0 if ab == 0 else 0.0     # both empty -> identical
    return ab / denom


def js_divergence(a, b, eps=1e-12):
    """Jensen-Shannon divergence in nats (0 = identical). Fields normalised to prob."""
    p = a.ravel().astype(np.float64); q = b.ravel().astype(np.float64)
    ps, qs = p.sum(), q.sum()
    if ps <= 0 or qs <= 0:
        return 0.0 if (ps <= 0 and qs <= 0) else float("inf")
    p = p / ps + eps; p /= p.sum()
    q = q / qs + eps; q /= q.sum()
    m = 0.5 * (p + q)
    kl = lambda x, y: float(np.sum(x * np.log(x / y)))
    return 0.5 * kl(p, m) + 0.5 * kl(q, m)


def binary_tanimoto(a, b, cutoff_frac):
    ma = a > (cutoff_frac * a.max()) if a.max() > 0 else np.zeros_like(a, bool)
    mb = b > (cutoff_frac * b.max()) if b.max() > 0 else np.zeros_like(b, bool)
    inter = np.logical_and(ma, mb).sum()
    union = np.logical_or(ma, mb).sum()
    return 1.0 if union == 0 else float(inter) / float(union)


def compare(paths, threshold=0.75, binary_cutoff=None, quiet=False):
    grids = []
    for p in paths:
        dens, voxvol, meta = load_grid(p)
        grids.append((p, dens, meta))
    shapes = {g[1].shape for g in grids}
    if len(shapes) != 1:
        raise ValueError(f"replicate grids differ in shape: "
                         f"{[(g[0], g[1].shape) for g in grids]}")
    # a field of ~0 total density is non-adsorbing, not "converged": two empty
    # replicates would score Tanimoto 1.0 and pass spuriously. Flag and fail them.
    totals = [float(g[1].sum()) for g in grids]
    empty = any(t <= 1e-9 for t in totals)
    pairs = []
    for (pa, a, _), (pb, b, _) in itertools.combinations(grids, 2):
        row = dict(a=pa, b=pb, tanimoto=tanimoto(a, b),
                   js_nats=js_divergence(a, b))
        row["js_bits"] = row["js_nats"] / np.log(2)
        if binary_cutoff is not None:
            row["binary_tanimoto"] = binary_tanimoto(a, b, binary_cutoff)
        pairs.append(row)
    # gate on continuous Tanimoto (DeepAPD criterion); binary is informational
    min_metric = min(r["tanimoto"] for r in pairs)
    ok = (not empty) and (min_metric >= threshold)
    result = dict(n_replicates=len(grids), gate="tanimoto", threshold=threshold,
                  min_metric=min_metric, ok=ok, empty=empty, pairs=pairs)
    if not quiet:
        _print(result)
    return result


def _print(r):
    print(f"replicates: {r['n_replicates']}   gate: {r['gate']} >= {r['threshold']}")
    for row in r["pairs"]:
        extra = f"  binT={row['binary_tanimoto']:.4f}" if "binary_tanimoto" in row else ""
        print(f"  {_short(row['a'])}  vs  {_short(row['b'])}: "
              f"Tanimoto={row['tanimoto']:.4f}  JS={row['js_bits']:.4f} bits{extra}")
    if r.get("empty"):
        verdict = "EMPTY / NON-ADSORBING (not a convergence result)"
    else:
        verdict = "CONVERGED (PASS)" if r["ok"] else "NOT CONVERGED (FAIL)"
    print(f"worst-pair {r['gate']} = {r['min_metric']:.4f}  ->  {verdict}")


def _short(p):
    # show the distinguishing run-dir plus the grid filename (the trailing
    # DensityGrid/System_0/<file> is identical across replicates, so hide it)
    if "DensityGrid" in p:
        run = p.split("/DensityGrid")[0].rstrip("/").split("/")[-1]
        return f"{run}/.../{p.split('/')[-1]}"
    parts = p.split("/")
    return "/".join(parts[-3:]) if len(parts) >= 3 else p


def _selftest():
    rng = np.random.default_rng(0)
    base = rng.random((16, 16, 16))
    # identical -> Tanimoto 1, JS 0
    r_id = compare_arrays([base, base.copy()])
    assert abs(r_id["min_metric"] - 1.0) < 1e-9, r_id["min_metric"]
    # near-identical (small noise) -> high Tanimoto, passes 0.75
    noisy = np.clip(base + 0.02 * rng.standard_normal(base.shape), 0, None)
    r_close = compare_arrays([base, noisy], threshold=0.75)
    assert r_close["ok"], r_close["min_metric"]
    # independent random -> should FAIL 0.75
    other = rng.random((16, 16, 16))
    r_far = compare_arrays([base, other], threshold=0.75)
    assert not r_far["ok"], r_far["min_metric"]
    print(f"field_convergence selftest: PASS "
          f"(identical={r_id['min_metric']:.3f}, close={r_close['min_metric']:.3f}, "
          f"independent={r_far['min_metric']:.3f})")


def compare_arrays(arrays, threshold=0.75, binary_cutoff=None):
    """In-memory variant for tests: arrays is a list of ndarrays."""
    pairs = []
    for a, b in itertools.combinations(arrays, 2):
        row = dict(tanimoto=tanimoto(a, b), js_nats=js_divergence(a, b))
        row["js_bits"] = row["js_nats"] / np.log(2)
        if binary_cutoff is not None:
            row["binary_tanimoto"] = binary_tanimoto(a, b, binary_cutoff)
        pairs.append(row)
    gate_key = "binary_tanimoto" if binary_cutoff is not None else "tanimoto"
    min_metric = min(r[gate_key] for r in pairs)
    return dict(n_replicates=len(arrays), gate=gate_key, threshold=threshold,
                min_metric=min_metric, ok=min_metric >= threshold, pairs=pairs)


def main(argv=None):
    ap = argparse.ArgumentParser(description="Replicate-cell field convergence (Tanimoto/JS).")
    ap.add_argument("grids", nargs="*", help="two or more .cube/.vtk grid files")
    ap.add_argument("--run-dirs", nargs="+", help="find one grid under each run dir instead")
    ap.add_argument("--threshold", type=float, default=0.75, help="pass if worst-pair gate metric >= this")
    ap.add_argument("--binary-cutoff", type=float, default=None,
                    help="also gate on binary Jaccard at this fraction-of-max occupancy (e.g. 0.05)")
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args(argv)
    if args.selftest:
        _selftest(); return 0
    paths = list(args.grids)
    if args.run_dirs:
        for d in args.run_dirs:
            g = find_grid(d)
            if not g:
                ap.error(f"no grid found under {d}")
            paths.append(g)
    if len(paths) < 2:
        ap.error("need >= 2 grids (or --run-dirs) to compare, or use --selftest")
    r = compare(paths, threshold=args.threshold, binary_cutoff=args.binary_cutoff)
    return 0 if r["ok"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
