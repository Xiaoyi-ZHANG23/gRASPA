#!/usr/bin/env python3
"""run_adaptive.py — gRASPA density fields with DYNAMIC (adaptive) cycles.

Instead of a fixed NumberOfProductionCycles, each MOF runs a replicate pair with a
high cycle CAP while gRASPA flushes the *cumulative* density grid every --checkpoint
cycles (DensityGridWriteEvery). A monitor compares the two replicates' live .cube
grids (continuous Tanimoto, triclinic-correct) and STOPS both as soon as they agree
(>= --threshold for --stable consecutive checks). Each field therefore runs for
just-enough cycles — no fixed-count waste — and falls back to the cap if it never
converges. Non-adsorbing (empty-grid) MOFs are detected and stopped early.

Why this works: the replicate-Tanimoto convergence criterion is evaluated live off
the checkpoints gRASPA already writes; killing a process just means we keep the last
good cumulative grid, which IS the converged field.

Examples:
  run_adaptive.py --mof qmof-d475ee7 --uc "1 1 1"
  run_adaptive.py --mof-list /home/nvidia/Smoke_Test/pilot_mofs.txt --gpus "0 1 2 3 4 5 6 7"
  run_adaptive.py --mof qmof-b334eed --uc "3 3 3" --checkpoint 200000 --max-cycles 6000000
  run_adaptive.py --mof-list ... --dry-run
"""
from __future__ import annotations
import argparse, datetime, os, re, shutil, subprocess, sys, threading, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from gridio import load_grid, find_grid
import field_cost

NVHPC = os.environ.get("NVHPC", "/opt/nvidia/hpc_sdk/Linux_x86_64/24.5")
LIBS = f"{NVHPC}/compilers/lib:{NVHPC}/cuda/lib64:{NVHPC}/math_libs/lib64"


def sh_env(gpu):
    e = dict(os.environ)
    e["LD_LIBRARY_PATH"] = LIBS + ":" + e.get("LD_LIBRARY_PATH", "")
    e["CUDA_VISIBLE_DEVICES"] = str(gpu)
    return e


def stage(run_dir, input_dir, cif_dir, name, uc, seed, init, prod, sample, checkpoint):
    os.makedirs(run_dir, exist_ok=True)
    for f in os.listdir(input_dir):
        if f.endswith(".def") or f == "simulation.input":
            shutil.copy(os.path.join(input_dir, f), run_dir)
    shutil.copy(os.path.join(cif_dir, f"{name}.cif"), run_dir)
    p = os.path.join(run_dir, "simulation.input")
    txt = open(p).read()
    subs = {
        r"^FrameworkName .*": f"FrameworkName {name}",
        r"^UnitCells .*": f"UnitCells 0 {uc}",
        r"^RandomSeed .*": f"RandomSeed  {seed}",
        r"^NumberOfInitializationCycles .*": f"NumberOfInitializationCycles {init}",
        r"^NumberOfProductionCycles .*": f"NumberOfProductionCycles     {prod}",
        r"^DensityGridSampleEvery .*": f"DensityGridSampleEvery   {sample}",
        r"^DensityGridWriteEvery .*": f"DensityGridWriteEvery    {checkpoint}",
    }
    for pat, rep in subs.items():
        txt = re.sub(pat, rep, txt, flags=re.MULTILINE)
    open(p, "w").write(txt)


def cube_snapshots(path):
    try:
        with open(path) as fh:
            fh.readline(); line2 = fh.readline()
        m = re.search(r"(\d+)\s+snapshots", line2)
        return int(m.group(1)) if m else None
    except Exception:
        return None


def read_grid_safe(run_dir):
    """Return (density, snapshots) or (None, None) if the checkpoint isn't readable yet."""
    try:
        g = find_grid(run_dir)
        if not g or not os.path.exists(g):
            return None, None
        dens, _vv, _m = load_grid(g)          # raises if the flush is mid-write (bad value count)
        return dens, cube_snapshots(g)
    except Exception:
        return None, None


def tanimoto(a, b):
    a = a.ravel(); b = b.ravel(); ab = float(np.dot(a, b))
    d = float(np.dot(a, a) + np.dot(b, b) - ab)
    return (ab / d) if d > 0 else (1.0 if ab == 0 else 0.0)


def grid_loading(run_dir, sites):
    """Loading (molecules) from the converged grid integral — gRASPA's own summary is
    absent because we kill the run at convergence, so we read it off the density field."""
    try:
        dens, vv, _ = load_grid(find_grid(run_dir))
        return float(dens.sum() * vv) / sites
    except Exception:
        return None


def run_field(name, uc, gpus, cfg, log):
    """Run one MOF's replicate pair adaptively on two GPUs. Returns a result dict."""
    sa, sb = cfg["seeds"]
    base = os.path.join(cfg["out"], name)
    dirs = {sa: f"{base}/s{sa}", sb: f"{base}/s{sb}"}
    procs = {}
    for seed, gpu in ((sa, gpus[0]), (sb, gpus[1])):
        shutil.rmtree(dirs[seed], ignore_errors=True)
        stage(dirs[seed], cfg["input_dir"], cfg["cif_dir"], name, uc, seed,
              cfg["init"], cfg["max_cycles"], cfg["sample"], cfg["checkpoint"])
        procs[seed] = subprocess.Popen(
            [cfg["bin"]], cwd=dirs[seed], env=sh_env(gpu),
            stdout=open(f"{dirs[seed]}/output.txt", "w"),
            stderr=open(f"{dirs[seed]}/run_stderr.txt", "w"))
    log(f"[{name}] launched seeds {sa}/{sb} on GPU {gpus[0]}/{gpus[1]} (cap {cfg['max_cycles']:,})")
    t0 = time.time()

    stable = 0; converged = False; final_tan = None; conv_snap = None; empty = False
    while True:
        time.sleep(cfg["poll"])
        alive = {s: (p.poll() is None) for s, p in procs.items()}
        da, na = read_grid_safe(dirs[sa]); db, nb = read_grid_safe(dirs[sb])
        if da is not None and db is not None and da.shape == db.shape:
            if da.sum() <= 1e-9 and db.sum() <= 1e-9 and min(na or 0, nb or 0) >= cfg["min_snap"]:
                empty = True
                log(f"[{name}] EMPTY grid after {na}/{nb} snapshots -> non-adsorbing, stopping")
                break
            t = tanimoto(da, db)
            stable = stable + 1 if t >= cfg["threshold"] else 0
            log(f"[{name}] snaps {na}/{nb}  Tanimoto={t:.3f}  stable={stable}/{cfg['stable']}")
            if stable >= cfg["stable"]:
                converged = True; final_tan = t; conv_snap = min(na or 0, nb or 0); break
        if not any(alive.values()):        # both hit the cap without converging
            if da is not None and db is not None:
                final_tan = tanimoto(da, db); conv_snap = min(na or 0, nb or 0)
            log(f"[{name}] both runs reached cap; Tanimoto={final_tan}")
            break

    wall = time.time() - t0
    for s, p in procs.items():
        if p.poll() is None:
            p.terminate()
            try: p.wait(timeout=30)
            except subprocess.TimeoutExpired: p.kill()
    # Cost/loading: we killed the runs, so gRASPA never printed "Work took"/loading.
    # Record measured wall-time and read loading straight off the converged grid.
    cyc = (conv_snap or 0) * cfg["sample"]
    inp = field_cost.parse_input(os.path.join(dirs[sa], "simulation.input"))
    ts = datetime.datetime.now().isoformat(timespec="seconds")
    loading = None
    for gpu, seed in zip(gpus, (sa, sb)):
        load = grid_loading(dirs[seed], cfg["sites"])
        if loading is None: loading = load
        row = {k: "" for k in field_cost.FIELDS}
        row.update(timestamp=ts, run_dir=dirs[seed], framework=name,
                   adsorbate=inp.get("adsorbate"), temperature_K=inp.get("temperature_K"),
                   pressure_Pa=inp.get("pressure_Pa"), replicate=seed, gpu_id=gpu, gpus=1,
                   init_cycles=cfg["init"], prod_cycles=cyc, snapshots=conv_snap,
                   work_seconds=round(wall, 1), gpu_hours=round(wall / 3600.0, 6),
                   loading_molecules=(round(load, 4) if load is not None else ""),
                   notes=f"adaptive {name} s{seed} {'conv' if converged else ('empty' if empty else 'cap')}")
        field_cost.append_csv(row, cfg["cost_csv"])
    return dict(mof=name, converged=converged, empty=empty, tanimoto=final_tan,
                conv_cycles=cyc, snapshots=conv_snap, loading=loading,
                gpu_h_field=round(2 * wall / 3600.0, 4), wall=round(wall, 1))


def main(argv=None):
    ap = argparse.ArgumentParser(description="Adaptive-cycle gRASPA density-field runner.")
    ap.add_argument("--mof"); ap.add_argument("--uc", help='unit cells "nx ny nz"')
    ap.add_argument("--mof-list", help="file: 'name nx ny nz' per line (# comments ok)")
    ap.add_argument("--seeds", default="1 999")
    ap.add_argument("--init", type=int, default=50000)
    ap.add_argument("--sample", type=int, default=10, help="DensityGridSampleEvery")
    ap.add_argument("--checkpoint", type=int, default=200000, help="DensityGridWriteEvery (convergence-check granularity)")
    ap.add_argument("--max-cycles", type=int, default=8000000, help="production-cycle cap")
    ap.add_argument("--threshold", type=float, default=0.75)
    ap.add_argument("--sites", type=int, default=1, help="adsorbate pseudo-atoms/molecule (Xe=1, CO2=3)")
    ap.add_argument("--stable", type=int, default=2, help="consecutive passing checks required")
    ap.add_argument("--min-snap", type=int, default=200, help="min snapshots before declaring EMPTY")
    ap.add_argument("--poll", type=int, default=20, help="seconds between convergence checks")
    ap.add_argument("--gpus", default="0 1")
    ap.add_argument("--input-dir", default="/home/nvidia/Smoke_Test/INPUT")
    ap.add_argument("--cif-dir", default="/home/nvidia/Smoke_Test/target")
    ap.add_argument("--out", default="/data/pilot/xe_adaptive")
    ap.add_argument("--cost-csv", default="/data/pilot/pilot_costs.csv")
    ap.add_argument("--tools", default=os.path.dirname(os.path.abspath(__file__)))
    ap.add_argument("--bin", default="/home/nvidia/gRASPA/src_clean/nvc_main.x")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args(argv)

    mofs = []
    if a.mof_list:
        for ln in open(a.mof_list):
            ln = ln.split("#")[0].split()
            if len(ln) >= 4: mofs.append((ln[0], f"{ln[1]} {ln[2]} {ln[3]}"))
    elif a.mof and a.uc:
        mofs.append((a.mof, a.uc))
    else:
        ap.error("give --mof + --uc, or --mof-list")
    seeds = a.seeds.split(); gpus = a.gpus.split()
    cfg = dict(seeds=seeds, init=a.init, sample=a.sample, checkpoint=a.checkpoint,
               max_cycles=a.max_cycles, threshold=a.threshold, sites=a.sites, stable=a.stable,
               min_snap=a.min_snap, poll=a.poll, input_dir=a.input_dir, cif_dir=a.cif_dir,
               out=a.out, cost_csv=a.cost_csv, tools=a.tools, bin=a.bin)
    print(f"adaptive: {len(mofs)} MOFs, seeds {seeds}, checkpoint {a.checkpoint:,}, "
          f"cap {a.max_cycles:,}, threshold {a.threshold}, {len(gpus)} GPUs ({len(gpus)//2} MOFs concurrent)")
    if a.dry_run:
        for n, uc in mofs: print(f"  {n}  uc={uc}")
        return 0
    os.makedirs(a.out, exist_ok=True)

    # GPU pool: 2 GPUs per MOF, up to len(gpus)//2 concurrent
    free = list(gpus); lock = threading.Lock(); results = []; rlock = threading.Lock()
    plog = lambda m: print(m, flush=True)
    def worker(name, uc):
        while True:
            with lock:
                if len(free) >= 2: pair = [free.pop(0), free.pop(0)]; break
            time.sleep(2)
        try:
            r = run_field(name, uc, pair, cfg, plog)
            with rlock: results.append(r)
        finally:
            with lock: free.extend(pair)
    threads = [threading.Thread(target=worker, args=(n, uc)) for n, uc in mofs]
    for t in threads: t.start()
    for t in threads: t.join()

    print("\n=== ADAPTIVE RESULTS ===")
    conv_csv = os.path.join(a.out, "adaptive_convergence.csv")
    tot_gpuh = 0.0
    with open(conv_csv, "w") as fh:
        fh.write("framework,status,tanimoto,loading,conv_cycles,gpu_h_field\n")
        for r in sorted(results, key=lambda r: r["mof"]):
            tan = f"{r['tanimoto']:.3f}" if r["tanimoto"] is not None else "na"
            load = f"{r['loading']:.1f}" if r.get("loading") is not None else "na"
            status = "EMPTY" if r["empty"] else ("CONVERGED" if r["converged"] else "CAP-REACHED")
            tot_gpuh += r.get("gpu_h_field", 0) or 0
            print(f"  {r['mof']:14s} {status:11s} T={tan:>6s} load={load:>6s} "
                  f"@ {r['conv_cycles']:>9,} cyc  {r.get('gpu_h_field',0):.4f} GPU-h")
            fh.write(f"{r['mof']},{status},{tan},{load},{r['conv_cycles']},{r.get('gpu_h_field',0)}\n")
    print(f"total: {tot_gpuh:.3f} GPU-h for {len(results)} fields  ->  {conv_csv} ; costs -> {a.cost_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
