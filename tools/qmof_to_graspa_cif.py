#!/usr/bin/env python3
"""qmof_to_graspa_cif.py — convert QMOF relaxed CIFs to gRASPA-parseable P1 CIFs
with DDEC framework charges written to _atom_site_charge (needed for CO2 + Ewald).

QMOF CIFs are already P1 (all atoms explicit) but store charges under varying,
non-standard column names and functionals. This picks the best available charge
column per file (DDEC preferred), writes a minimal gRASPA CIF, and reports the
net framework charge (should be ~0 for a neutral framework).

Charge preference (best -> fallback):
  hse06_10hf_ddec > pbe_ddec > hle17_ddec > pbe_cm5 > pbe_bader

Usage:
  qmof_to_graspa_cif.py --list ids.txt --src /data/qmof/relaxed_structures --out /data/qmof/graspa_cifs
  qmof_to_graspa_cif.py --cif one.cif --out /tmp   # single file (prints net charge)
"""
from __future__ import annotations
import argparse, os, re

CHARGE_PREF = ["_atom_site_hse06_10hf_ddec_charge", "_atom_site_pbe_ddec_charge",
               "_atom_site_hle17_ddec_charge", "_atom_site_pbe_cm5_charge",
               "_atom_site_pbe_bader_charge"]


def _cellval(lines, key):
    for l in lines:
        s = l.split()
        if len(s) >= 2 and s[0].lower() == key:
            return re.sub(r"\(.*?\)", "", s[1])
    return None


def convert(cif_path):
    """Return (cif_text, n_atoms, net_charge, charge_source) or raise ValueError."""
    lines = open(cif_path, errors="replace").read().splitlines()
    cell = {k: _cellval(lines, "_cell_" + k) for k in
            ("length_a", "length_b", "length_c", "angle_alpha", "angle_beta", "angle_gamma")}
    if any(v is None for v in cell.values()):
        raise ValueError("missing cell params")
    # locate the atom_site loop: a loop_ whose header block contains _atom_site_fract_x
    hdr, data_start = None, None
    i = 0
    while i < len(lines):
        if lines[i].strip() == "loop_":
            j = i + 1; heads = []
            while j < len(lines) and lines[j].strip().startswith("_"):
                heads.append(lines[j].strip().lower()); j += 1
            if "_atom_site_fract_x" in heads:
                hdr, data_start = heads, j; break
            i = j
        else:
            i += 1
    if hdr is None:
        raise ValueError("no atom_site loop")
    idx = {h: k for k, h in enumerate(hdr)}
    def col(name): return idx.get(name)
    ti = col("_atom_site_type_symbol"); li = col("_atom_site_label")
    xi, yi, zi = col("_atom_site_fract_x"), col("_atom_site_fract_y"), col("_atom_site_fract_z")
    if None in (ti, xi, yi, zi):
        raise ValueError("missing type/fract columns")
    ci = next((col(c) for c in CHARGE_PREF if col(c) is not None), None)
    csrc = next((c for c in CHARGE_PREF if col(c) is not None), "none")
    atoms = []
    for l in lines[data_start:]:
        s = l.split()
        if not s or s[0].startswith("_") or s[0] == "loop_" or s[0].startswith("#"):
            if atoms:
                break
            continue
        if len(s) < len(hdr):
            break
        elem = re.sub(r"[0-9+\-].*$", "", s[ti]) or s[ti]   # strip any oxidation/charge suffix
        label = s[li] if li is not None else f"{elem}{len(atoms)}"
        q = float(s[ci]) if ci is not None else 0.0
        atoms.append((label, elem, s[xi], s[yi], s[zi], q))
    if not atoms:
        raise ValueError("no atoms parsed")
    net = sum(a[5] for a in atoms)
    out = []
    out.append(f"data_{os.path.splitext(os.path.basename(cif_path))[0]}")
    out.append("_symmetry_space_group_name_H-M   'P 1'")
    out.append("_symmetry_Int_Tables_number   1")
    out.append(f"_cell_length_a   {cell['length_a']}")
    out.append(f"_cell_length_b   {cell['length_b']}")
    out.append(f"_cell_length_c   {cell['length_c']}")
    out.append(f"_cell_angle_alpha   {cell['angle_alpha']}")
    out.append(f"_cell_angle_beta   {cell['angle_beta']}")
    out.append(f"_cell_angle_gamma   {cell['angle_gamma']}")
    out.append("loop_")
    out.append("_symmetry_equiv_pos_as_xyz")
    out.append("  'x, y, z'")
    out.append("loop_")
    out.append("_atom_site_label")
    out.append("_atom_site_type_symbol")
    out.append("_atom_site_fract_x")
    out.append("_atom_site_fract_y")
    out.append("_atom_site_fract_z")
    out.append("_atom_site_charge")
    for label, elem, x, y, z, q in atoms:
        out.append(f"{label:8s} {elem:4s} {x:>18s} {y:>18s} {z:>18s} {q:12.6f}")
    return "\n".join(out) + "\n", len(atoms), net, csrc


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--cif"); ap.add_argument("--list"); ap.add_argument("--src", default="/data/qmof/relaxed_structures")
    ap.add_argument("--out", required=True)
    ap.add_argument("--neutral-tol", type=float, default=0.10, help="warn if |net charge| > this")
    a = ap.parse_args(argv)
    os.makedirs(a.out, exist_ok=True)
    if a.cif:
        items = [a.cif]
    elif a.list:
        items = [os.path.join(a.src, l.strip() + ".cif") for l in open(a.list) if l.strip()]
    else:
        ap.error("give --cif or --list")
    ok = warn = fail = 0; srcs = {}
    log = open(os.path.join(a.out, "charge_report.csv"), "w")
    log.write("qmof_id,n_atoms,net_charge,charge_source,status\n")
    for cif in items:
        qid = os.path.splitext(os.path.basename(cif))[0]
        try:
            text, nat, net, csrc = convert(cif)
        except Exception as e:
            fail += 1; log.write(f"{qid},,,,FAIL:{e}\n"); continue
        open(os.path.join(a.out, qid + ".cif"), "w").write(text)
        srcs[csrc] = srcs.get(csrc, 0) + 1
        status = "ok" if abs(net) <= a.neutral_tol else "NONNEUTRAL"
        if status != "ok": warn += 1
        else: ok += 1
        log.write(f"{qid},{nat},{net:.4f},{csrc},{status}\n")
    log.close()
    print(f"converted: {ok} ok, {warn} non-neutral (|net|>{a.neutral_tol}), {fail} failed -> {a.out}")
    print(f"charge sources: {srcs}")


if __name__ == "__main__":
    main()
