#!/usr/bin/env python3
"""field_cost.py — log the real per-field GCMC cost of one gRASPA run.

This is deliverable "measure the real per-field cost on GPU". It parses one run's
inputs + log and appends a row to a CSV so the pilot produces a measured
cost-per-density-field number (the figure the multi-thousand-hour sweep is sized on).

What it reads:
  simulation.input  -> FrameworkName, adsorbate MoleculeName, Temperature, Pressure,
                       Init/Eq/Production cycles, DensityGridPoints
  run log (stderr)  -> "Work took N seconds"  (gRASPA's own compute time = the GPU cost)
                       "(nx x ny x nz, K snapshots)"                (grid sampling depth)
                       "LOADING: # MOLECULES" section -> adsorbate Overall Average/ErrorBar
                       optional GNU `/usr/bin/time -v` "Elapsed (wall clock)" line

GPU-hours = work_seconds * gpus / 3600  (default gpus=1; one sim per GPU).

Dependency-free (stdlib only).

Examples:
  field_cost.py --run-dir /data/pilot/smoke_densitygrid --log run_stderr.txt \
                --replicate 0 --gpu-id 0 --csv /data/pilot/pilot_costs.csv
  field_cost.py --selftest
"""
from __future__ import annotations
import argparse
import csv
import datetime as _dt
import os
import re
import sys

FIELDS = [
    "timestamp", "run_dir", "framework", "adsorbate", "temperature_K", "pressure_Pa",
    "replicate", "gpu_id", "gpus", "init_cycles", "eq_cycles", "prod_cycles",
    "grid_points", "snapshots", "work_seconds", "elapsed_seconds", "gpu_hours",
    "loading_molecules", "loading_error", "notes",
]


def _grep_input(text, key):
    m = re.search(rf"^\s*{re.escape(key)}\s+(\S+)", text, re.MULTILINE | re.IGNORECASE)
    return m.group(1) if m else None


def parse_input(path):
    d = {}
    if not path or not os.path.exists(path):
        return d
    t = open(path).read()
    d["framework"] = _grep_input(t, "FrameworkName")
    d["temperature_K"] = _grep_input(t, "Temperature")
    d["pressure_Pa"] = _grep_input(t, "Pressure")
    d["init_cycles"] = _grep_input(t, "NumberOfInitializationCycles")
    d["eq_cycles"] = _grep_input(t, "NumberOfEquilibrationCycles")
    d["prod_cycles"] = _grep_input(t, "NumberOfProductionCycles")
    m = re.search(r"MoleculeName\s+(\S+)", t)          # first adsorbate component
    d["adsorbate"] = m.group(1) if m else None
    g = re.search(r"DensityGridPoints\s+(\d+)\s+(\d+)\s+(\d+)", t)
    d["grid_points"] = "x".join(g.groups()) if g else None
    return d


def _elapsed_to_seconds(s):
    parts = s.split(":")
    try:
        parts = [float(p) for p in parts]
    except ValueError:
        return None
    sec = 0.0
    for p in parts:
        sec = sec * 60 + p
    return sec


def parse_log(path, adsorbate=None):
    d = {}
    if not path or not os.path.exists(path):
        return d
    t = open(path, errors="replace").read()
    m = re.search(r"Work took\s+([\d.]+)\s+seconds", t)
    d["work_seconds"] = float(m.group(1)) if m else None
    # GNU `/usr/bin/time -v` label has several colons ("(h:mm:ss or m:ss):"), so take
    # the substring after the final ": " on that line rather than a colon-greedy regex.
    eln = next((ln for ln in t.splitlines() if "Elapsed (wall clock)" in ln), None)
    d["elapsed_seconds"] = _elapsed_to_seconds(eln.rsplit(": ", 1)[-1].strip()) if eln else None
    m = re.search(r"\(\s*\d+\s*x\s*\d+\s*x\s*\d+\s*,\s*(\d+)\s+snapshots\)", t)
    d["snapshots"] = int(m.group(1)) if m else None
    d["loading_molecules"], d["loading_error"] = _parse_loading(t, adsorbate)
    return d


def _parse_loading(text, adsorbate=None):
    """Return (avg, err) for the adsorbate from the 'LOADING: # MOLECULES' block."""
    lines = text.splitlines()
    start = next((i for i, ln in enumerate(lines)
                  if "LOADING: # MOLECULES" in ln or "LOADING:# MOLECULES" in ln), None)
    if start is None:
        return None, None
    comp_re = re.compile(r"COMPONENT\s*\[(\d+)\]\s*\(([^)]+)\)")
    over_re = re.compile(r"Overall:\s*Average:\s*([-\d.eE]+),\s*ErrorBar:\s*([-\w.]+)")
    entries = []          # (index, name, avg, err)
    cur = None
    for ln in lines[start + 1:]:
        if "BLOCK AVERAGES" in ln and "# MOLECULES" not in ln:
            break         # next loading-unit section begins
        cm = comp_re.search(ln)
        if cm:
            cur = (int(cm.group(1)), cm.group(2))
            continue
        om = over_re.search(ln)
        if om and cur is not None:
            try:
                err = float(om.group(2))
            except ValueError:
                err = float("nan")           # '-nan' etc.
            entries.append((cur[0], cur[1], float(om.group(1)), err))
            cur = None
    if not entries:
        return None, None
    # prefer exact adsorbate-name match; else first non-framework (index>0); else max avg
    if adsorbate:
        for idx, name, avg, err in entries:
            if name.strip().lower() == adsorbate.strip().lower():
                return avg, err
    non_frame = [e for e in entries if e[0] >= 1]
    pick = max(non_frame or entries, key=lambda e: e[2])
    return pick[2], pick[3]


def collect(run_dir=None, log=None, sim_input=None, replicate=None, gpu_id=None,
            gpus=1, notes=""):
    if run_dir:
        sim_input = sim_input or os.path.join(run_dir, "simulation.input")
        if log and not os.path.isabs(log):
            log = os.path.join(run_dir, log)
        if not log:                      # auto-find a captured log
            for cand in ("run_stderr.txt", "run.log", "output.txt"):
                p = os.path.join(run_dir, cand)
                if os.path.exists(p):
                    log = p; break
    row = {k: "" for k in FIELDS}
    row["timestamp"] = _dt.datetime.now().isoformat(timespec="seconds")
    row["run_dir"] = run_dir or (os.path.dirname(log) if log else "")
    row["replicate"] = "" if replicate is None else replicate
    row["gpu_id"] = "" if gpu_id is None else gpu_id
    row["gpus"] = gpus
    row["notes"] = notes
    inp = parse_input(sim_input)
    row.update({k: v for k, v in inp.items() if v is not None})
    lg = parse_log(log, adsorbate=inp.get("adsorbate"))
    row.update({k: v for k, v in lg.items() if v is not None})
    ws = row.get("work_seconds")
    if isinstance(ws, (int, float)):
        row["gpu_hours"] = round(ws * gpus / 3600.0, 6)
    return row


def append_csv(row, csv_path):
    os.makedirs(os.path.dirname(os.path.abspath(csv_path)), exist_ok=True)
    new = (not os.path.exists(csv_path)) or os.path.getsize(csv_path) == 0
    with open(csv_path, "a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        if new:
            w.writeheader()
        w.writerow({k: row.get(k, "") for k in FIELDS})


def _selftest():
    import tempfile
    d = tempfile.mkdtemp()
    with open(os.path.join(d, "simulation.input"), "w") as fh:
        fh.write("FrameworkName MFI-2x2x2-P1\nTemperature 298\nPressure 1e4\n"
                 "NumberOfInitializationCycles 3000\nNumberOfEquilibrationCycles 0\n"
                 "NumberOfProductionCycles 3000\nDensityGridPoints 40 40 40\n"
                 "Component 0 MoleculeName CO2\n")
    with open(os.path.join(d, "run_stderr.txt"), "w") as fh:
        fh.write(
            "Work took 12.5 seconds\n"
            "Wrote 3D density grid for component 1 (CO2): a.cube / .vtk (40 x 40 x 40, 3000 snapshots)\n"
            "=================== BLOCK AVERAGES (LOADING: # MOLECULES)=============\n"
            "COMPONENT [0] (MFI-2x2x2-P1.cif)\nOverall: Average: 1.00000, ErrorBar: 0.00000\n"
            "COMPONENT [1] (CO2)\nOverall: Average: 19.06567, ErrorBar: 1.51800\n"
            "=====================BLOCK AVERAGES (LOADING: mg/g)=============\n"
            "COMPONENT [1] (CO2)\nOverall: Average: 1000.0, ErrorBar: -nan\n"
            "\tElapsed (wall clock) time (h:mm:ss or m:ss): 0:15.19\n")
    row = collect(run_dir=d, replicate=0, gpu_id=0, gpus=1, notes="selftest")
    assert row["framework"] == "MFI-2x2x2-P1", row
    assert row["adsorbate"] == "CO2", row
    assert row["prod_cycles"] == "3000", row
    assert row["snapshots"] == 3000, row
    assert abs(row["work_seconds"] - 12.5) < 1e-9, row
    assert abs(row["loading_molecules"] - 19.06567) < 1e-6, row          # NOT the framework 1.0 or mg/g 1000
    assert abs(row["loading_error"] - 1.518) < 1e-6, row
    assert abs(row["gpu_hours"] - round(12.5 / 3600, 6)) < 1e-9, row
    assert abs(row["elapsed_seconds"] - 15.19) < 1e-6, row
    csv_path = os.path.join(d, "costs.csv")
    append_csv(row, csv_path); append_csv(row, csv_path)
    with open(csv_path) as fh:
        lines = fh.read().strip().splitlines()
    assert len(lines) == 3, lines            # header + 2 rows
    print(f"field_cost selftest: PASS "
          f"(loading={row['loading_molecules']}, work={row['work_seconds']}s, "
          f"gpu_hours={row['gpu_hours']:.6f})")


def main(argv=None):
    ap = argparse.ArgumentParser(description="Log per-field gRASPA cost to a CSV.")
    ap.add_argument("--run-dir", help="gRASPA run directory (has simulation.input + a captured log)")
    ap.add_argument("--log", help="log file (captured stderr); default: auto-find in run-dir")
    ap.add_argument("--input", dest="sim_input", help="simulation.input path (default: in run-dir)")
    ap.add_argument("--replicate", type=int, default=None)
    ap.add_argument("--gpu-id", type=int, default=None)
    ap.add_argument("--gpus", type=int, default=1, help="GPUs used by this sim (for gpu_hours)")
    ap.add_argument("--notes", default="")
    ap.add_argument("--csv", default="/data/pilot/pilot_costs.csv")
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--dry-run", action="store_true", help="print the row, do not write CSV")
    args = ap.parse_args(argv)
    if args.selftest:
        _selftest(); return 0
    if not (args.run_dir or args.log):
        ap.error("need --run-dir or --log (or --selftest)")
    row = collect(run_dir=args.run_dir, log=args.log, sim_input=args.sim_input,
                  replicate=args.replicate, gpu_id=args.gpu_id, gpus=args.gpus, notes=args.notes)
    print("  ".join(f"{k}={row[k]}" for k in
                    ("framework", "adsorbate", "temperature_K", "pressure_Pa",
                     "prod_cycles", "snapshots", "work_seconds", "gpu_hours",
                     "loading_molecules", "loading_error") if row.get(k) != ""))
    if not args.dry_run:
        append_csv(row, args.csv)
        print(f"appended -> {args.csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
