#!/usr/bin/env python3
"""validate_grid.py — sanity-check a gRASPA 3D density grid.

This is the `density_grid_ab/validate_grid.py` the DensityGrid example README
refers to but that was never shipped. The core identity it checks:

    integral of the density field over the cell  ==  average adsorbate-ATOM count
    sum(voxel) * voxelVolume  ==  <N_atoms>       (molecules = atoms / sites)

For the reference case (CO2 in MFI, seed 0) this is 57.197 atoms = 19.06567 CO2.

It also reports basic non-degeneracy stats (min/max/mean, zero-voxel fraction) and
can optionally compare the integral to a known loading and save a slice PNG.

numpy only, unless --png is given (then matplotlib is imported lazily).

Examples:
  validate_grid.py DensityGrid/System_0/DensityGrid_1_CO2.vtk --sites 3 --expect 19.06567
  validate_grid.py --run-dir /data/pilot/run_042 --png slice.png
  validate_grid.py --selftest
"""
from __future__ import annotations
import argparse
import sys
import numpy as np

try:
    from gridio import load_grid, find_grid
except ImportError:                                   # allow running from elsewhere
    sys.path.insert(0, __file__.rsplit("/", 1)[0])
    from gridio import load_grid, find_grid


def validate(path, sites=3, expect=None, tol=1.0, png=None, quiet=False):
    dens, voxvol, meta = load_grid(path)
    if meta["format"] == "vtk" and not quiet:
        print("NOTE: integrating a .vtk — its SPACING assumes an orthorhombic cell; "
              "for a triclinic framework use the .cube (triclinic-correct voxel volume).")
    total_atoms = float(dens.sum() * voxvol)
    total_mols = total_atoms / sites
    stats = dict(
        path=path, fmt=meta["format"], dims=meta["dims"], voxvol=voxvol,
        vmin=float(dens.min()), vmax=float(dens.max()), vmean=float(dens.mean()),
        zero_frac=float((dens <= 0).mean()),
        integral_atoms=total_atoms, integral_molecules=total_mols,
    )
    ok = True
    reasons = []
    # non-degeneracy
    if stats["vmax"] <= 0:
        ok = False; reasons.append("grid is all-zero (degenerate)")
    if not np.isfinite(dens).all():
        ok = False; reasons.append("grid contains NaN/Inf")
    if stats["zero_frac"] >= 0.999:
        ok = False; reasons.append("virtually every voxel is zero")
    # optional loading comparison
    if expect is not None:
        err_pct = (total_mols - expect) / expect * 100.0 if expect else float("inf")
        stats["expect_molecules"] = expect
        stats["error_pct"] = err_pct
        if abs(err_pct) > tol:
            ok = False
            reasons.append(f"integral {total_mols:.5f} vs expected {expect:.5f} "
                           f"({err_pct:+.3f}% > {tol}% tol)")
    stats["ok"] = ok
    stats["reasons"] = reasons
    if not quiet:
        _print(stats)
    if png:
        _save_slices(dens, png, title=path.split("/")[-1])
        if not quiet:
            print(f"saved slice PNG -> {png}")
    return stats


def _print(s):
    print(f"grid      : {s['path']}  [{s['fmt']}]  dims={s['dims']}  voxvol={s['voxvol']:.5f} A^3")
    print(f"stats     : min={s['vmin']:.3e}  max={s['vmax']:.3e}  mean={s['vmean']:.3e}  "
          f"zero_frac={s['zero_frac']*100:.1f}%")
    print(f"integral  : {s['integral_atoms']:.4f} atoms -> {s['integral_molecules']:.5f} molecules")
    if "expect_molecules" in s:
        print(f"expected  : {s['expect_molecules']:.5f} molecules   error={s['error_pct']:+.3f}%")
    print(f"verdict   : {'PASS' if s['ok'] else 'FAIL'}" +
          ("" if s['ok'] else "  (" + "; ".join(s['reasons']) + ")"))


def _save_slices(dens, png, title=""):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    nx, ny, nz = dens.shape
    zproj = dens.sum(axis=2)          # (x,y)
    midz = dens[:, :, nz // 2]        # (x,y)
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.6))
    for a, img, t in ((ax[0], zproj, "z-integrated density"),
                      (ax[1], midz, f"mid-z plane (z={nz//2})")):
        im = a.imshow(img.T, origin="lower", cmap="magma", aspect="auto")
        a.set_title(t, fontsize=10); a.set_xlabel("x voxel"); a.set_ylabel("y voxel")
        fig.colorbar(im, ax=a, fraction=0.046)
    fig.suptitle(title, fontsize=11); fig.tight_layout()
    fig.savefig(png, dpi=110)


def _selftest():
    """Build a tiny known grid, write it as .vtk, and confirm the integral identity."""
    import tempfile, os
    nx = ny = nz = 4
    sx, sy, sz = 0.5, 0.5, 0.5                # A -> voxvol 0.125 A^3
    dens = np.zeros((nx, ny, nz))
    dens[1, 1, 1] = 8.0                       # one hot voxel: integral = 8*0.125 = 1.0 atom
    # write VTK (x-fastest). canonical dens[x,y,z] -> flat with idx=x+nx*(y+ny*z)
    flat = dens.transpose(2, 1, 0).reshape(-1)   # (z,y,x) then flatten = x-fastest
    d = tempfile.mkdtemp()
    p = os.path.join(d, "t.vtk")
    with open(p, "w") as fh:
        fh.write("# vtk DataFile Version 3.0\nselftest\nASCII\nDATASET STRUCTURED_POINTS\n")
        fh.write(f"DIMENSIONS {nx} {ny} {nz}\nORIGIN 0 0 0\nSPACING {sx} {sy} {sz}\n")
        fh.write(f"POINT_DATA {nx*ny*nz}\nSCALARS number_density double 1\nLOOKUP_TABLE default\n")
        fh.write("\n".join(f"{v:.8e}" for v in flat) + "\n")
    s = validate(p, sites=2, expect=0.5, tol=1e-6, quiet=True)   # 1 atom / 2 sites = 0.5 mol
    assert abs(s["integral_atoms"] - 1.0) < 1e-9, s["integral_atoms"]
    assert abs(s["integral_molecules"] - 0.5) < 1e-9, s["integral_molecules"]
    assert s["ok"], s["reasons"]
    print("validate_grid selftest: PASS "
          f"(integral={s['integral_atoms']:.6f} atoms, {s['integral_molecules']:.6f} mol)")


def main(argv=None):
    ap = argparse.ArgumentParser(description="Validate a gRASPA density grid (integral identity + stats).")
    ap.add_argument("grid", nargs="?", help=".cube or .vtk grid file")
    ap.add_argument("--run-dir", help="find the grid under a run directory instead")
    ap.add_argument("--sites", type=int, default=3, help="adsorbate pseudo-atoms per molecule (CO2=3)")
    ap.add_argument("--expect", type=float, default=None, help="expected loading (molecules) to compare against")
    ap.add_argument("--tol", type=float, default=1.0, help="allowed %% error vs --expect (default 1%%)")
    ap.add_argument("--png", help="save a 2-panel slice PNG here")
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args(argv)
    if args.selftest:
        _selftest(); return 0
    path = args.grid or (find_grid(args.run_dir) if args.run_dir else None)
    if not path:
        ap.error("give a grid file, or --run-dir, or --selftest")
    s = validate(path, sites=args.sites, expect=args.expect, tol=args.tol, png=args.png)
    return 0 if s["ok"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
