#!/usr/bin/env python3
"""prescreen.py — Phase-1 Stage 0: cheap loading pre-screen of the QMOF library.

Per MOF: 2 grid-OFF Xe replicates (distinct NONZERO seeds), short cycles. Emits per
MOF -> screen.csv: loading_mean, loading_cv, adsorbing, n_atoms, orthorhombic,
unit_cells, and a pre-sized adaptive cap (cap_M) from the Week-0 loading->cycles fit.

Frozen config (Phase1_Handoff): Xe, 273 K, 25 bar, ChargeMethod None, grid OFF here.
Dispatches (MOF x seed) jobs across GPUs via a thread pool (one job per GPU at a time).

Usage:
  prescreen.py --mofs qmof-4836169,qmof-7e3c643 --gpus "0 1 2 3 4 5 6 7"   # validate a few
  prescreen.py --list /home/nvidia/Smoke_Test/target/mof_list.txt --gpus "0 1 2 3 4 5 6 7"
  prescreen.py --all --gpus "0 1 2 3 4 5 6 7"                              # full library
"""
from __future__ import annotations
import argparse, math, os, queue, re, shutil, subprocess, sys, threading, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import field_cost  # for parse_input / stdlib CSV schema reuse (not required)

NVHPC = os.environ.get("NVHPC", "/opt/nvidia/hpc_sdk/Linux_x86_64/24.5")
LIBS = f"{NVHPC}/compilers/lib:{NVHPC}/cuda/lib64:{NVHPC}/math_libs/lib64"
COLS = ["framework", "loading_mean", "loading_cv", "adsorbing", "n_atoms",
        "orthorhombic", "unit_cells", "cap_M", "seeds", "note"]


def parse_geom(path):
    """o.3812858 -> {name: {'uc':'nx ny nz', 'ortho':bool}} (4-line blocks)."""
    lines = [l.rstrip() for l in open(path) if l.strip()]
    g = {}; i = 0
    while i + 3 < len(lines) + 1:
        if i + 3 >= len(lines):
            break
        head = lines[i].split()
        if not head or not head[0].endswith(".cif"):
            i += 1; continue
        name = head[1] if len(head) > 1 else head[0][:-4]
        try:
            cp = list(map(float, lines[i + 1].split()))       # a al b be c ga
            uc = re.findall(r"\d+", lines[i + 3])
        except Exception:
            i += 1; continue
        if len(cp) >= 6 and len(uc) >= 3:
            angles = [cp[1], cp[3], cp[5]]
            ortho = all(abs(a - 90.0) < 0.5 for a in angles)
            g[name] = {"uc": f"{uc[0]} {uc[1]} {uc[2]}", "ortho": ortho}
        i += 4
    return g


def _cross(u, w): return (u[1]*w[2]-u[2]*w[1], u[2]*w[0]-u[0]*w[2], u[0]*w[1]-u[1]*w[0])
def _norm(u): return math.sqrt(u[0]**2 + u[1]**2 + u[2]**2)
def _dot(u, w): return u[0]*w[0] + u[1]*w[1] + u[2]*w[2]


def read_cif_cell(path):
    keys = ("_cell_length_a", "_cell_length_b", "_cell_length_c",
            "_cell_angle_alpha", "_cell_angle_beta", "_cell_angle_gamma")
    v = {}
    for line in open(path, errors="replace"):
        s = line.split()
        if len(s) >= 2 and s[0].lower() in keys:
            v[s[0].lower()] = float(re.sub(r"\(.*?\)", "", s[1]))
        if len(v) == 6:
            break
    if len(v) < 6:
        return None
    return tuple(v[k] for k in keys)


def cif_geom(cif_path, cutoff=12.8):
    """Unit cells (min-image over 2*cutoff) + orthorhombic flag from a CIF's cell.
    Validated against o.3812858 (15/15 exact). Returns {'uc','ortho'} or None."""
    cell = read_cif_cell(cif_path)
    if not cell:
        return None
    a, b, c, al, be, ga = cell
    ar, br, gr = math.radians(al), math.radians(be), math.radians(ga)
    ax = (a, 0.0, 0.0)
    bx = (b*math.cos(gr), b*math.sin(gr), 0.0)
    cx = c*math.cos(br)
    cy = c*(math.cos(ar) - math.cos(br)*math.cos(gr))/math.sin(gr)
    cz = math.sqrt(max(0.0, c*c - cx*cx - cy*cy))
    cv = (cx, cy, cz)
    V = abs(_dot(ax, _cross(bx, cv)))
    if V <= 0:
        return None
    wa = V/_norm(_cross(bx, cv)); wb = V/_norm(_cross(ax, cv)); wc = V/_norm(_cross(ax, bx))
    need = 2.0*cutoff
    uc = f"{max(1, math.ceil(need/wa))} {max(1, math.ceil(need/wb))} {max(1, math.ceil(need/wc))}"
    return {"uc": uc, "ortho": all(abs(x-90.0) < 0.5 for x in (al, be, ga))}


def cap_from_loading(L, tstop=0.85):
    """Pre-sized adaptive cap (millions of cycles): 2x expected converged cycles at tstop."""
    base = max(0.5, 0.012 * L - 0.49)          # cycles_M at T=0.75 (Week-0 fit)
    return round(2.0 * base * (0.25 / (1.0 - tstop)), 2)


def _write_csv(path, names, geom, results, adsorb_thresh, tstop, seeds_str):
    import csv as _csv
    nads = 0
    tmp = path + ".tmp"
    with open(tmp, "w", newline="") as fh:
        w = _csv.DictWriter(fh, fieldnames=COLS); w.writeheader()
        for n in names:
            r = results.get(n, {"loads": [], "nat": ""}); loads = r["loads"]
            if loads:
                mean = sum(loads) / len(loads)
                std = (sum((x - mean) ** 2 for x in loads) / len(loads)) ** 0.5
                cv = (std / mean) if mean > 1e-9 else 0.0
                ads = mean >= adsorb_thresh
            else:
                mean = cv = 0.0; ads = False
            nads += bool(ads)
            w.writerow({"framework": n, "loading_mean": round(mean, 4), "loading_cv": round(cv, 4),
                        "adsorbing": ads, "n_atoms": r["nat"], "orthorhombic": geom[n]["ortho"],
                        "unit_cells": geom[n]["uc"], "cap_M": cap_from_loading(mean, tstop),
                        "seeds": seeds_str, "note": "" if loads else "no_output"})
    os.replace(tmp, path)   # atomic
    return nads


def env(gpu):
    e = dict(os.environ); e["LD_LIBRARY_PATH"] = LIBS + ":" + e.get("LD_LIBRARY_PATH", "")
    e["CUDA_VISIBLE_DEVICES"] = str(gpu); return e


def run_one(cfg, name, uc, seed, gpu):
    """One grid-OFF Xe replicate. Returns (loading, n_atoms) or (None, None)."""
    R = os.path.join(cfg["out"], "runs", name, f"s{seed}")
    shutil.rmtree(R, ignore_errors=True); os.makedirs(R, exist_ok=True)
    for f in os.listdir(cfg["input_dir"]):
        if f.endswith(".def") or f == "simulation.input":
            shutil.copy(os.path.join(cfg["input_dir"], f), R)
    cif = os.path.join(cfg["cif_dir"], f"{name}.cif")
    if not os.path.exists(cif):
        return None, None
    shutil.copy(cif, R)
    p = os.path.join(R, "simulation.input"); txt = open(p).read()
    for pat, rep in {
        r"^FrameworkName .*": f"FrameworkName {name}",
        r"^UnitCells .*": f"UnitCells 0 {uc}",
        r"^RandomSeed .*": f"RandomSeed  {seed}",
        r"^NumberOfInitializationCycles .*": f"NumberOfInitializationCycles {cfg['init']}",
        r"^NumberOfProductionCycles .*": f"NumberOfProductionCycles     {cfg['prod']}",
        r"^ComputeDensityGrid .*": "ComputeDensityGrid       no",   # grid OFF for screen
        r"^UseChargesFromCIFFile .*": f"UseChargesFromCIFFile {cfg['charges']}",  # Xe:no, CO2:yes
    }.items():
        txt = re.sub(pat, rep, txt, flags=re.MULTILINE)
    open(p, "w").write(txt)
    try:
        subprocess.run([cfg["bin"]], cwd=R, env=env(gpu),
                       stdout=open(f"{R}/output.txt", "w"), stderr=open(f"{R}/err.txt", "w"),
                       timeout=cfg["timeout"], check=False)
    except subprocess.TimeoutExpired:
        return None, None
    out = open(f"{R}/output.txt", errors="replace").read()
    load = _parse_loading(out); nat = _parse_natoms(out)
    if not cfg["keep"]:
        shutil.rmtree(R, ignore_errors=True)
    return load, nat


def _parse_loading(text, adsorbate=None):
    """Adsorbate loading from the '# MOLECULES' block. Gas-agnostic: the framework is
    component [0]; the adsorbate is the highest-index component (or an exact name match)."""
    lines = text.splitlines()
    start = next((i for i, l in enumerate(lines) if "LOADING: # MOLECULES" in l), None)
    if start is None:
        return None
    entries = []; cur = None
    for l in lines[start + 1:]:
        if "BLOCK AVERAGES" in l and "# MOLECULES" not in l:
            break
        m = re.search(r"COMPONENT\s*\[(\d+)\]\s*\(([^)]+)\)", l)
        if m:
            cur = (int(m.group(1)), m.group(2)); continue
        o = re.search(r"Overall:\s*Average:\s*([-\d.eE]+)", l)
        if o and cur is not None:
            entries.append((cur[0], cur[1], float(o.group(1)))); cur = None
    if not entries:
        return None
    if adsorbate:
        for _idx, name, val in entries:
            if name.strip().lower() == adsorbate.strip().lower():
                return val
    non_frame = [e for e in entries if e[0] >= 1]      # framework is component 0
    return (non_frame or entries)[-1][2]


def _parse_natoms(text):
    m = re.search(r"Total atoms:\s*(\d+)", text)
    return int(m.group(1)) if m else ""


def worker(cfg, q, results, rlock, log):
    gpu = cfg["gpu_for_thread"][threading.current_thread().name]
    while True:
        try:
            name, uc, seed = q.get_nowait()
        except queue.Empty:
            return
        if cfg["abort"].is_set():
            q.task_done(); return
        load, nat = run_one(cfg, name, uc, seed, gpu)
        with rlock:
            results.setdefault(name, {"loads": [], "nat": ""})
            if load is not None:
                results[name]["loads"].append(load)
            if nat:
                results[name]["nat"] = nat
            cfg["counters"]["attempted"] += 1
            cfg["counters"]["parsed"] += int(load is not None)
            done = sum(len(v["loads"]) for v in results.values())
        if done % 50 == 0:
            log(f"  ...{done} replicate runs done")
        q.task_done()


def main(argv=None):
    ap = argparse.ArgumentParser(description="Phase-1 Stage 0 loading pre-screen (Xe, grid off).")
    ap.add_argument("--mofs", help="comma-separated MOF names (validation)")
    ap.add_argument("--list", help="file with one MOF name per line")
    ap.add_argument("--all", action="store_true", help="every MOF in the geometry table with a CIF")
    ap.add_argument("--limit", type=int, default=0, help="cap number of MOFs (0=all)")
    ap.add_argument("--seeds", default="7 999")
    ap.add_argument("--init", type=int, default=10000)
    ap.add_argument("--prod", type=int, default=10000)
    ap.add_argument("--charges", choices=["yes", "no"], default="no",
                    help="UseChargesFromCIFFile: no for Xe (chargeless), yes for CO2 (Ewald + framework charges)")
    ap.add_argument("--adsorb-thresh", type=float, default=0.5, help="loading_mean >= this => adsorbing")
    ap.add_argument("--tstop", type=float, default=0.85, help="production Tanimoto stop (for cap sizing)")
    ap.add_argument("--timeout", type=int, default=300)
    ap.add_argument("--gpus", default="0 1 2 3 4 5 6 7")
    ap.add_argument("--geom", default="/home/nvidia/Smoke_Test/target/o.3812858")
    ap.add_argument("--geom-mode", choices=["table", "cif"], default="cif",
                    help="cif: compute unit cells per CIF (QMOF, default); table: parse o.3812858")
    ap.add_argument("--cutoff", type=float, default=12.8)
    ap.add_argument("--input-dir", default="/home/nvidia/Smoke_Test/INPUT")
    ap.add_argument("--cif-dir", default="/home/nvidia/Smoke_Test/target")
    ap.add_argument("--out", default="/data/pilot/prescreen")
    ap.add_argument("--csv", default="/data/pilot/prescreen/screen.csv")
    ap.add_argument("--bin", default="/home/nvidia/gRASPA/src_clean/nvc_main.x")
    ap.add_argument("--keep", action="store_true", help="keep per-run dirs (default: delete)")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args(argv)

    if a.mofs:
        cand = [n.strip() for n in a.mofs.split(",") if n.strip()]
    elif a.list:
        cand = [l.split("#")[0].strip() for l in open(a.list) if l.split("#")[0].strip()]
    elif a.all:
        if a.geom_mode != "table":
            ap.error("--all needs --geom-mode table; for cif mode pass --list of MOF ids")
        cand = sorted(parse_geom(a.geom).keys())
    else:
        ap.error("give --mofs, --list, or --all")
    if a.geom_mode == "table":
        table = parse_geom(a.geom); geom = {n: table[n] for n in cand if n in table}
    else:  # compute unit cells + ortho from each CIF
        geom = {}
        for n in cand:
            cif = os.path.join(a.cif_dir, n + ".cif")
            if os.path.exists(cif):
                g = cif_geom(cif, a.cutoff)
                if g:
                    geom[n] = g
    names = [n for n in cand if n in geom and os.path.exists(os.path.join(a.cif_dir, n + ".cif"))]
    if len(names) < len(cand):
        print(f"note: {len(cand)-len(names)} of {len(cand)} candidates dropped (missing CIF or unparseable cell)")
    if a.limit:
        names = names[:a.limit]
    seeds = a.seeds.split(); gpus = a.gpus.split()
    print(f"pre-screen: {len(names)} MOFs x {len(seeds)} seeds = {len(names)*len(seeds)} runs, "
          f"{a.init}/{a.prod} cyc grid-off, {len(gpus)} GPUs, seeds {seeds}")
    if a.dry_run:
        for n in names[:20]:
            print(f"  {n}  uc={geom[n]['uc']}  ortho={geom[n]['ortho']}")
        return 0

    cfg = dict(out=a.out, input_dir=a.input_dir, cif_dir=a.cif_dir, bin=a.bin,
               init=a.init, prod=a.prod, timeout=a.timeout, keep=a.keep, charges=a.charges)
    cfg["counters"] = {"attempted": 0, "parsed": 0}; cfg["abort"] = threading.Event()
    os.makedirs(os.path.dirname(a.csv), exist_ok=True)
    q = queue.Queue()
    for n in names:
        for s in seeds:
            q.put((n, geom[n]["uc"], s))
    results = {}; rlock = threading.Lock(); log = lambda m: print(m, flush=True)
    cfg["gpu_for_thread"] = {}
    threads = []
    for i, gpu in enumerate(gpus):
        t = threading.Thread(target=worker, args=(cfg, q, results, rlock, log), name=f"gpu{gpu}")
        cfg["gpu_for_thread"][t.name] = gpu
        threads.append(t)
    t0 = time.time()
    for t in threads: t.start()
    stop = threading.Event()
    def _ckpt():                       # periodic atomic checkpoint (crash-safe) + sanity abort
        while not stop.wait(60):
            with rlock:
                _write_csv(a.csv, names, geom, results, a.adsorb_thresh, a.tstop, "|".join(seeds))
                att = cfg["counters"]["attempted"]; par = cfg["counters"]["parsed"]
            if att >= 40 and par == 0:
                print(f"SANITY ABORT: first {att} runs all produced no parseable loading "
                      f"(parser/config bug?). Stopping instead of burning the full sweep.", flush=True)
                cfg["abort"].set(); stop.set(); return
    ck = threading.Thread(target=_ckpt, daemon=True); ck.start()
    for t in threads: t.join()
    stop.set()
    dt = time.time() - t0

    with rlock:
        nads = _write_csv(a.csv, names, geom, results, a.adsorb_thresh, a.tstop, "|".join(seeds))
    print(f"\n=== pre-screen done: {len(names)} MOFs in {dt/60:.1f} min ===")
    print(f"adsorbing: {nads}/{len(names)}  ({100*nads/max(len(names),1):.0f}%)  -> {a.csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
