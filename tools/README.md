# tools/ — density-field convergence, cost, and adaptive cycling

Post-processing and orchestration for the 3D adsorbate **density grid**
(`ComputeDensityGrid`). These decide when a spatial density field ρ(r) is
*converged* (not just when scalar uptake is), measure the real GPU cost per field,
and run each MOF for **just enough cycles** to converge — no fixed-count waste.

Dependency-light: `numpy` only (matplotlib is imported lazily, and only if you ask
for a slice PNG). Every script has `--help` and a `--selftest`.

## The tools

| script | what it does |
|---|---|
| `gridio.py` | Read a `.cube`/`.vtk` density grid into a canonical `[x,y,z]` array. **Defaults to `.cube`** (see gotcha below). |
| `validate_grid.py` | Integrate a grid and check `∫ρ·dV == avg atom count` (÷ sites → molecules). Non-degeneracy stats; optional slice PNG. |
| `field_convergence.py` | Compare ≥2 **replicate** grids: continuous Tanimoto (DeepAPD criterion) + Jensen–Shannon. PASS iff worst-pair Tanimoto ≥ threshold (default 0.75). Flags empty/non-adsorbing grids. |
| `field_cost.py` | Parse a run's log → wall-time, GPU-hours, cycles, loading; append a row to a costs CSV. |
| `run_pilot.sh` | Fixed-cycle pilot: dispatch MOF×seed runs across GPUs (GPU-slot semaphore, no scheduler), log cost, gate convergence. |
| `cycle_map.sh` | Map Tanimoto vs production cycles (build a cost-vs-loading model). |
| `run_adaptive.py` | **Adaptive (dynamic) cycling — run each field to just-enough convergence.** |

## Adaptive cycling — `run_adaptive.py`

Instead of a fixed `NumberOfProductionCycles`, each MOF runs a replicate pair with a
high cycle **cap** while gRASPA flushes the *cumulative* density grid every
`--checkpoint` cycles (`DensityGridWriteEvery`). A monitor compares the two
replicates' live `.cube` grids and **stops both the moment they hold Tanimoto ≥
`--threshold` for `--stable` consecutive checkpoints**. This rides on the periodic
grid writes gRASPA already does — no gRASPA source change, no restart.

Because the runs are killed at convergence, gRASPA never prints its end-of-run
summary, so the runner **measures wall-time itself and reads loading from the
converged grid integral** (`--sites`: Xe=1, CO2=3). Non-adsorbing (empty-grid) MOFs
are detected and stopped early.

```bash
# one MOF
tools/run_adaptive.py --mof qmof-d475ee7 --uc "1 1 1" --gpus "0 1"

# a set (2 GPUs per MOF, len(gpus)//2 MOFs concurrent)
tools/run_adaptive.py --mof-list mofs.txt --gpus "0 1 2 3 4 5 6 7" \
    --checkpoint 200000 --max-cycles 8000000 --threshold 0.75 --sites 1
```
Writes `adaptive_convergence.csv` (status, Tanimoto, loading, converged-cycles,
GPU-h/field) and appends per-replicate cost rows to the costs CSV.

`mofs.txt` lines: `framework_name  nx ny nz` (`# comments` ok).

## Gotchas (learned the hard way — heed these)

- **Integrate the `.cube`, not the `.vtk`, for non-orthorhombic frameworks.** The
  `.vtk` STRUCTURED_POINTS geometry assumes an orthorhombic cell; its declared
  SPACING overstates the voxel volume by `|a||b||c| / V_cell`, so re-integrating a
  **triclinic** MOF's `.vtk` overcounts (measured +20% on a triclinic QMOF). The
  `.cube` voxel vectors are triclinic-correct. Density *values* are identical in both.
- **Independent replicates need distinct *nonzero* seeds — never pair `0` with `1`.**
  glibc `srand(0)` and `srand(1)` produce the identical stream, so `RandomSeed 0`
  and `1` are byte-identical runs.
- **Converged uptake ≠ converged field.** Scalar loading converges far sooner than
  ρ(r); always gate on the replicate field metric, never on uptake alone. Cycles to
  field-convergence rise with loading (empirically ≈ 0.012·loading million cycles).
