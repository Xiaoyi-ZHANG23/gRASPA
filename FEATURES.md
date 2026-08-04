# What this fork adds (vs. [snurr-group/gRASPA](https://github.com/snurr-group/gRASPA))

All development happens on the single feature branch **`feature/energy-grid-ewald`**.
Every feature is opt-in from `simulation.input` — with all keywords absent, the code
behaves exactly like upstream gRASPA.

| Feature | Enable with | Output |
|---|---|---|
| 3D adsorbate density grid | `ComputeDensityGrid yes` | `DensityGrid/System_<i>/*.cube` + `.vtk` |
| Framework energy grid | `ComputeEnergyGrid yes` | `EnergyGrid/System_<i>/*.grid` + `.cube` |
| Widom insertion-energy grid | `ComputeWidomEnergyGrid yes` | `WidomEnergyGrid/System_<i>/*.cube` + `.dat` |
| Field-workflow toolkit | — | `tools/` (CLI scripts) |

## 1. 3D adsorbate density grid

Accumulates the time-averaged adsorbate number density ρ(r) on a voxel grid during
GCMC production (RASPA2 `ComputeDensityProfile3DVTKGrid` parity).

```
ComputeDensityGrid       yes
DensityGridPoints        40 40 40
DensityGridSampleEvery   1        # sample every N cycles
DensityGridWriteEvery    100000   # periodic flush of the cumulative grid
DensityGridMolecularCenter yes    # optional: bin one point per molecule (its center)
                                  # instead of one point per atom site
```

Per adsorbate component it writes a Gaussian `.cube` (true triclinic voxel vectors,
geometry in Bohr) and a `.vtk` STRUCTURED_POINTS file. The `.cube` integrates back to
the simulation's reported average loading (validated to 0.000%). **For triclinic
frameworks always integrate the `.cube`** — the `.vtk` format assumes an orthorhombic
cell and overstates the voxel volume by `|a||b||c|/V_cell`.

Regression case: `Examples/DensityGrid/` (CO2 in MFI 2×2×2, deterministic seed, ships a
vetted reference `output.txt`; grid integral = 19.06567 CO2 = block-average loading).

## 2. Framework energy grid

Precomputes the rigid-framework potential over one unit cell, once, before any MC
cycle: a VdW grid per probe pseudo-atom type plus an Ewald Coulomb grid (real **and**
reciprocal space) for a unit +1e probe. Bit-for-bit RASPA2 `SimulationType MakeGrid`
parity — a RASPA2 MakeGrid input can be reused directly:

```
ComputeEnergyGrid  yes
SpacingVDWGrid     1.0            # Å
SpacingCoulombGrid 1.0            # Å
GridTypes          C_co2 O_co2    # probe pseudo-atom types
EnergyGridOverlapCriteria 1.0e7   # optional energy cap (RASPA2 parity)
```

Outputs ASCII `.grid` files (RASPA2 column layout, `?` markers at capped/overlap
points) and `.cube` files (values in K) for VMD/OVITO. Full specification — formats,
units, grid-geometry math, implementation map, and validation record — in
[`ENERGY_GRIDS.md`](ENERGY_GRIDS.md). Example: `Examples/EnergyGrid/`.

## 3. Widom insertion-energy grid

A spatial map of Widom test-insertion energies and Rosenbluth weights, binned by the
sample's position — the 3D generalization of RASPA2's 1D `ComputeFreeEnergyProfile`.
Includes the full Ewald contribution. Accumulated during production, written at the end:

```
ComputeWidomEnergyGrid   yes
WidomEnergyGridPoints    50 50 50
WidomEnergyGridRawDump   yes      # optional: one line per successful sample
```

Writes `*_EnergyW.cube`, `*_BetaF.cube` (local free-energy map), `*_EwaldW.cube`, and
optionally the raw per-sample dump. Details in [`ENERGY_GRIDS.md`](ENERGY_GRIDS.md).

## 4. Field-workflow toolkit (`tools/`)

CLI scripts for running density-field campaigns (each has `--help`; most have `--selftest`):

- `gridio.py` / `validate_grid.py` — parse `.cube`/`.vtk` grids; gate a run on
  *grid integral == reported loading*.
- `field_convergence.py` — replicate-field convergence between two independent runs
  (continuous/binary Tanimoto, symmetric KL). A converged **uptake** does *not* imply a
  converged **field**; this tool gates on the field itself.
- `field_cost.py` — per-field GPU-hour cost logger (CSV).
- `run_adaptive.py` — adaptive cycling: runs a replicate pair with a high cycle cap,
  monitors the periodically flushed grids, and stops both runs as soon as the field
  convergence gate holds (large compute savings vs. any fixed cycle count).
- `run_pilot.sh` / `cycle_map.sh` — multi-GPU batch orchestration (plain GPU-slot
  semaphore, no scheduler required).
- `qmof_manifest.py` / `qmof_to_graspa_cif.py` / `prescreen.py` — QMOF high-throughput
  pipeline: manifest assembly, CIF conversion with DDEC partial-charge injection, and a
  cheap loading pre-screen to filter non-adsorbing structures.

## 5. Documentation and fixes

- [`ENERGY_GRIDS.md`](ENERGY_GRIDS.md) — full spec of both energy-grid features.
- [`KNOWN_ISSUES.md`](KNOWN_ISSUES.md) — documents a pre-existing **upstream** Ewald
  staging bug affecting single-bead charged species (not introduced or fixed here).
- `src_clean/NVC_COMPILE` — link flags updated to match the NVIDIA HPC SDK 24.5
  toolchain (the CUDA library path previously pointed at a 22.5 install).
