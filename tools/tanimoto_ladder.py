#!/usr/bin/env python3
"""tanimoto_ladder.py — label-quality ladder in ONE adaptive run per MOF.

For each MOF, run a replicate pair (distinct nonzero seeds) under a cycle cap while the
live cumulative density grids are compared (continuous Tanimoto). As T sweeps upward,
SNAPSHOT both replicates' .cube grids the moment each rung threshold is confirmed
(T >= rung at two polls separated by a new checkpoint), then continue to the next rung.
One run therefore yields labels at every rung for the price of the tightest one, and the
rungs differ ONLY in sampling noise — exactly what a label-quality ladder needs.

Outputs under <out>/<mof>/:
  T<rung>/s<seed>.cube  — snapshot grids per rung (both replicates)
  ladder.json           — per-rung cycles, achieved T, wall seconds; final status
Appends one line per MOF to <out>/ladder_results.jsonl.

Cap per MOF is pre-sized from the pre-screen loading via the Week-0 fit
  base_M = max(0.5, 0.0122*L - 0.49);  cycles(T) = base_M * 0.25/(1-T)
  cap_M  = clamp(2 * cycles(topT) * premium, 3, 40)   (premium: CO2 unknowns)
"""
from __future__ import annotations
import argparse, csv, json, os, shutil, subprocess, sys, threading, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from run_adaptive import stage, read_grid_safe, tanimoto, sh_env, grid_loading
from gridio import find_grid
import prescreen


def cap_from_loading(L, top_t, premium):
    base = max(0.5, 0.0122 * L - 0.49)
    return min(40.0, max(3.0, 2.0 * base * 0.25 / (1.0 - top_t) * premium))


def run_ladder(name, uc, gpus, cfg, log):
    sa, sb = cfg["seeds"]
    base = os.path.join(cfg["out"], name)
    dirs = {sa: f"{base}/s{sa}", sb: f"{base}/s{sb}"}
    cap_cycles = int(cfg["caps"][name] * 1e6)
    procs = {}
    for seed, gpu in ((sa, gpus[0]), (sb, gpus[1])):
        shutil.rmtree(dirs[seed], ignore_errors=True)
        stage(dirs[seed], cfg["input_dir"], cfg["cif_dir"], name, uc, seed,
              cfg["init"], cap_cycles, cfg["sample"], cfg["checkpoint"])
        procs[seed] = subprocess.Popen(
            [cfg["bin"]], cwd=dirs[seed], env=sh_env(gpu),
            stdout=open(f"{dirs[seed]}/output.txt", "w"),
            stderr=open(f"{dirs[seed]}/run_stderr.txt", "w"))
    log(f"[{cfg['gas']}:{name}] pair on GPU {gpus[0]}/{gpus[1]}, cap {cap_cycles/1e6:.1f}M")
    t0 = time.time()
    rungs = cfg["rungs"]; ri = 0
    pending_na = None                       # na at the first pass of the current rung
    rows = []; empty = False; last_t = None; last_na = 0
    while ri < len(rungs):
        time.sleep(cfg["poll"])
        alive = any(p.poll() is None for p in procs.values())
        da, na = read_grid_safe(dirs[sa]); db, nb = read_grid_safe(dirs[sb])
        if da is not None and db is not None and da.shape == db.shape:
            if da.sum() <= 1e-9 and db.sum() <= 1e-9 and min(na or 0, nb or 0) >= cfg["min_snap"]:
                empty = True; log(f"[{cfg['gas']}:{name}] EMPTY -> non-adsorbing"); break
            t = tanimoto(da, db); last_t = t; last_na = min(na or 0, nb or 0)
            rung = rungs[ri]
            if t >= rung:
                if pending_na is None:
                    pending_na = last_na               # first pass at this rung
                elif last_na > pending_na:             # confirmed on a NEW checkpoint
                    rd = f"{base}/T{rung:g}"; os.makedirs(rd, exist_ok=True)
                    for seed in (sa, sb):
                        g = find_grid(dirs[seed])
                        if not g: continue
                        for attempt in range(2):      # a copy can race the checkpoint flush:
                            try:                      # validate the grid parses before keeping it
                                from gridio import load_grid as _lg; _lg(g)
                                shutil.copy(g, f"{rd}/s{seed}.cube"); break
                            except Exception:
                                time.sleep(cfg["poll"])
                    rows.append(dict(rung=rung, cycles=last_na * cfg["sample"],
                                     tanimoto=round(t, 4), wall_s=round(time.time() - t0, 1)))
                    log(f"[{cfg['gas']}:{name}] rung {rung:g} @ {last_na*cfg['sample']/1e6:.2f}M cyc (T={t:.3f})")
                    ri += 1; pending_na = None
            else:
                pending_na = None                      # dipped below: reset confirmation
        if not alive:                                   # cap reached (or died)
            log(f"[{cfg['gas']}:{name}] cap reached at rung {rungs[ri]:g} (T={last_t})")
            break
    wall = time.time() - t0
    for p in procs.values():
        if p.poll() is None:
            p.terminate()
            try: p.wait(timeout=30)
            except subprocess.TimeoutExpired: p.kill()
    load = grid_loading(dirs[sa], 1)                    # molecular-center grid -> molecules
    rec = dict(mof=name, gas=cfg["gas"], empty=empty, rungs=rows,
               top_rung_reached=(len(rows) == len(rungs)), last_tanimoto=last_t,
               last_cycles=last_na * cfg["sample"], loading=load,
               gpu_h=round(2 * wall / 3600.0, 4), wall_s=round(wall, 1))
    json.dump(rec, open(f"{base}/ladder.json", "w"), indent=1)
    with cfg["jl_lock"]:
        with open(f"{cfg['out']}/ladder_results.jsonl", "a") as fh:
            fh.write(json.dumps(rec) + "\n")
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mofs", required=True, help="file with one MOF id per line")
    ap.add_argument("--gas", required=True, choices=["Xe", "CO2"])
    ap.add_argument("--screen", required=True, help="pre-screen csv (loading -> cap sizing)")
    ap.add_argument("--input-dir", required=True)
    ap.add_argument("--cif-dir", default="/data/qmof/graspa_cifs")
    ap.add_argument("--rungs", default="0.75,0.85,0.92,0.96")
    ap.add_argument("--premium", type=float, default=1.0, help="cap multiplier (CO2 unknowns)")
    ap.add_argument("--seeds", default="7 999")
    ap.add_argument("--init", type=int, default=50000)
    ap.add_argument("--sample", type=int, default=10)
    ap.add_argument("--checkpoint", type=int, default=200000)
    ap.add_argument("--poll", type=int, default=10)
    ap.add_argument("--min-snap", type=int, default=200)
    ap.add_argument("--gpus", required=True, help='e.g. "0 1 2 3"')
    ap.add_argument("--out", required=True)
    ap.add_argument("--bin", default="/home/nvidia/gRASPA/src_clean/nvc_main.x")
    ap.add_argument("--resume", action="store_true", help="skip MOFs that already have ladder.json")
    a = ap.parse_args()

    mofs = [l.strip() for l in open(a.mofs) if l.strip()]
    loads = {r["framework"]: float(r["loading_mean"]) for r in csv.DictReader(open(a.screen))}
    rungs = sorted(float(x) for x in a.rungs.split(","))
    caps = {m: cap_from_loading(loads.get(m, 100.0), rungs[-1], a.premium) for m in mofs}
    geom = {}
    for m in list(mofs):
        g = prescreen.cif_geom(os.path.join(a.cif_dir, m + ".cif"))
        if g: geom[m] = g["uc"]
        else: mofs.remove(m)
    os.makedirs(a.out, exist_ok=True)
    if a.resume:
        mofs = [m for m in mofs if not os.path.exists(f"{a.out}/{m}/ladder.json")]
    gpus = a.gpus.split()
    est = sum(caps[m] for m in mofs) / 2 * 0.042  # rough: expected ~half the cap, GPU-h/M-cyc pair
    print(f"{a.gas} ladder: {len(mofs)} MOFs, rungs {rungs}, {len(gpus)} GPUs "
          f"({len(gpus)//2} pairs concurrent); rough est ~{est:.0f} GPU-h", flush=True)

    cfg = dict(out=a.out, input_dir=a.input_dir, cif_dir=a.cif_dir, bin=a.bin, gas=a.gas,
               seeds=a.seeds.split(), init=a.init, sample=a.sample, checkpoint=a.checkpoint,
               poll=a.poll, min_snap=a.min_snap, rungs=rungs, caps=caps,
               jl_lock=threading.Lock())
    free = list(gpus); lock = threading.Lock()
    log = lambda m: print(m, flush=True)
    def worker(name):
        while True:
            with lock:
                if len(free) >= 2:
                    pair = [free.pop(0), free.pop(0)]; break
            time.sleep(3)
        try:
            run_ladder(name, geom[name], pair, cfg, log)
        except Exception as e:
            log(f"[{a.gas}:{name}] FAILED: {type(e).__name__}: {e}")
        finally:
            with lock: free.extend(pair)
    threads = [threading.Thread(target=worker, args=(m,)) for m in mofs]
    for t in threads: t.start()
    for t in threads: t.join()
    done = sum(1 for m in mofs if os.path.exists(f"{a.out}/{m}/ladder.json"))
    print(f"{a.gas} ladder complete: {done}/{len(mofs)} MOFs -> {a.out}", flush=True)


if __name__ == "__main__":
    main()
