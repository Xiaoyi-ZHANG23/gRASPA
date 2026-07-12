#!/usr/bin/env python3
"""qmof_manifest.py — Phase-1 Stage D (assemble) + Stage P (PLD screen) for QMOF.

Reads the QMOF properties table (qmof.csv) and emits:
  qmof_manifest.csv  — id, name, formula, natoms, pld, lcd, density, volume,
                       crystal_system, angles_all_90 (bool), source, charge_source
  screen_pld.csv     — id, pld, retained, reason, xe_accessible, co2_accessible
  retained_mofs.txt  — retained working set (PLD >= CO2 floor), one id per line

Guest kinetic diameters: CO2 ~ 3.30 A, Xe ~ 3.96 A. Hard nonporous drop PLD < 2.4 A.
The .cube is always used for integration; angles_all_90 is metadata only (True only
for orthorhombic/tetragonal/cubic — ~12% of QMOF; the rest MUST use .cube).
"""
from __future__ import annotations
import argparse, csv, os

ORTHO_SYS = {"orthorhombic", "tetragonal", "cubic"}   # all cell angles == 90


def fnum(x):
    try: return float(x)
    except Exception: return None


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default="/data/qmof/qmof_database/qmof.csv")
    ap.add_argument("--out", default="/data/qmof")
    ap.add_argument("--nonporous", type=float, default=2.4, help="hard drop PLD < this")
    ap.add_argument("--co2-floor", type=float, default=3.30, help="CO2 kinetic diameter (retain floor)")
    ap.add_argument("--xe-floor", type=float, default=3.96, help="Xe kinetic diameter")
    a = ap.parse_args(argv)
    rows = list(csv.DictReader(open(a.csv)))
    os.makedirs(a.out, exist_ok=True)

    man = open(os.path.join(a.out, "qmof_manifest.csv"), "w", newline="")
    mw = csv.DictWriter(man, fieldnames=["qmof_id", "name", "formula", "natoms", "pld", "lcd",
                                         "density", "volume", "crystal_system", "angles_all_90",
                                         "source", "charge_source"])
    mw.writeheader()
    scr = open(os.path.join(a.out, "screen_pld.csv"), "w", newline="")
    sw = csv.DictWriter(scr, fieldnames=["qmof_id", "pld", "retained", "reason",
                                         "xe_accessible", "co2_accessible"])
    sw.writeheader()
    retained = []
    n_ret = n_xe = n_co2only = n_nonporous = n_below = n_badpld = 0
    for r in rows:
        qid = r["qmof_id"]; cs = r.get("info.symmetry.spacegroup_crystal", "")
        pld = fnum(r.get("info.pld"))
        mw.writerow({
            "qmof_id": qid, "name": r.get("name", ""), "formula": r.get("info.formula", ""),
            "natoms": r.get("info.natoms", ""), "pld": r.get("info.pld", ""),
            "lcd": r.get("info.lcd", ""), "density": r.get("info.density", ""),
            "volume": r.get("info.volume", ""), "crystal_system": cs,
            "angles_all_90": cs in ORTHO_SYS, "source": r.get("info.source", ""),
            "charge_source": "QMOF-DDEC(hse06_10hf/hle17)"})
        if pld is None:
            reason, ret = "no_pld", False; n_badpld += 1
        elif pld < a.nonporous:
            reason, ret = "nonporous", False; n_nonporous += 1
        elif pld < a.co2_floor:
            reason, ret = "below_co2_floor", False; n_below += 1
        else:
            reason, ret = "retained", True; n_ret += 1
        xe = ret and pld >= a.xe_floor
        co2 = ret
        if ret:
            retained.append(qid)
            if xe: n_xe += 1
            else: n_co2only += 1
        sw.writerow({"qmof_id": qid, "pld": r.get("info.pld", ""), "retained": ret,
                     "reason": reason, "xe_accessible": xe, "co2_accessible": co2})
    man.close(); scr.close()
    open(os.path.join(a.out, "retained_mofs.txt"), "w").write("\n".join(retained) + "\n")

    N = len(rows)
    print(f"Stage D: {N} QMOF structures -> qmof_manifest.csv")
    print(f"Stage P (PLD screen):")
    print(f"  dropped nonporous (PLD<{a.nonporous}): {n_nonporous} ({100*n_nonporous/N:.1f}%)")
    print(f"  dropped below CO2 floor ({a.nonporous}<=PLD<{a.co2_floor}): {n_below} ({100*n_below/N:.1f}%)")
    if n_badpld: print(f"  dropped no/invalid PLD: {n_badpld}")
    print(f"  RETAINED (PLD>={a.co2_floor}): {n_ret} ({100*n_ret/N:.1f}%)")
    print(f"    Xe+CO2 accessible (PLD>={a.xe_floor}): {n_xe}")
    print(f"    CO2-only ({a.co2_floor}<=PLD<{a.xe_floor}): {n_co2only}")
    print(f"-> retained_mofs.txt, screen_pld.csv in {a.out}")


if __name__ == "__main__":
    main()
