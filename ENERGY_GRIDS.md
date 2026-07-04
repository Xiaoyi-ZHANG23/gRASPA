# Energy-Grid Features (branch `feature/energy-grid-ewald`)

This branch adds two grid features on top of `feature/density-grid-3d`:

| Feature | Keyword | What it computes | When | Output dir |
|---|---|---|---|---|
| **Framework energy grid** | `ComputeEnergyGrid yes` | Precomputed framework potential over one unit cell: VDW per probe type + Ewald Coulomb (real **and** reciprocal) for a unit +1e probe. RASPA2 `SimulationType MakeGrid` parity. | Once, before any MC cycle | `EnergyGrid/System_<i>/` |
| **Widom energy grid** | `ComputeWidomEnergyGrid yes` | Spatial map of Widom-insertion energies (incl. Ewald) and Rosenbluth weights, binned by sample position. 3D generalization of RASPA2's 1D `ComputeFreeEnergyProfile`. | Accumulated during production, written at end | `WidomEnergyGrid/System_<i>/` |

(The base branch's density grid — `ComputeDensityGrid`, RASPA2 `ComputeDensityProfile3DVTKGrid` parity — is documented in `Examples/DensityGrid/README.md`.)

Both features default **OFF** and are provably surgical: with the keywords absent, output
is byte-identical to the base branch (verified with `sweep_compare.sh` and an empty
expected-diff manifest across CO2-MFI, Henrys_coefficient, Tail-Correction).

**Units.** gRASPA's internal energy unit is 10 J/mol. Multiply by `1.2027242847`
(`Units::energy_to_kelvin`, `src_clean/data_struct.h`) to get Kelvin. All energies in the
files described below are **already converted to Kelvin**. The Ewald prefactor
`Boxsize::Prefactor = 138935.483496` equals RASPA2's `COULOMBIC_CONVERSION_FACTOR`.

---

## 1. Framework energy grid (`ComputeEnergyGrid`)

### Input keywords (RASPA2 names — a RASPA2 MakeGrid input can be reused directly)

| keyword | default | meaning |
|---|---|---|
| `ComputeEnergyGrid yes` | off | enable |
| `SpacingVDWGrid <x>` | 0.15 Å | requested VDW grid spacing |
| `SpacingCoulombGrid <x>` | 0.15 Å | requested Coulomb grid spacing |
| `GridTypes <name> [name...]` | none | pseudo-atom probes, one VDW grid per name. The number of names defines the number of grids (RASPA2's `NumberOfGrids` count keyword is unnecessary and ignored). |
| `EnergyGridOverlapCriteria <x>` | 1.0e7 | cap, **internal** energy units (RASPA2 `EnergyOverlapCriteria` parity) |

The Coulomb grid is made automatically when charges are on (`ChargeMethod Ewald`,
charged atoms present); skipped under `noCharges`.

### Grid geometry — RASPA2 `grids.c` parity, bit for bit

- Enclosing rectangular box of **one unit cell**: `Size_axis = |a_axis|+|b_axis|+|c_axis|`
  (component-wise absolute values); negative components accumulate into `Shift`.
  Unit-cell vectors are recovered as `Box.Cell rows / NumberofUnitCells`.
- `N = (int)(Size/Spacing)`, bumped to the next **odd** number; `Delta = Size/N`.
- `N+1` inclusive planes per axis (plane N is the periodic image of plane 0).
- Point `(i,j,k)` at `pos = Shift + (i·Δx, j·Δy, k·Δz)`; loop order i(x) outer → j(y) → k(z).
- Energies evaluated with minimum-image over the **full supercell** (like RASPA2, which
  grids one unit cell but sums over the replicated box).

### Output files

**`EnergyGrid_VDW_<name>.grid`** — ASCII, RASPA2 `MakeASCIGrid` parity, no header:
```
x  y  z  E  Fx  Fy  Fz          (Å, K, K/Å; F = -gradient, analytic)
x  y  z  ?  ?   ?   ?           (where E >= EnergyGridOverlapCriteria)
```
`awk '{print $4}'` yields the single energy column with `?` overlap markers — the shape
consumed by the PSED-CNN ASCII pipeline.

**`EnergyGrid_Coulomb_Ewald.grid`** — ASCII, 6 columns, unit +1e probe:
```
x  y  z  E_total  E_real  E_recip     (K)
x  y  z  ?        ?       ?           (E_real > cap  OR  nearest framework atom < 1 Å)
```
- `E_real` = `Σ Prefactor·q_j·erfc(α·r)/r` within `CutOffCoulomb` — **exactly** the
  quantity RASPA2 `MakeGrid` tabulates in its Coulomb grid
  (`CalculateDerivativesAtPositionReal`, unit probe hardcoded).
- `E_recip` = framework↔probe reciprocal-space cross term, computed from the framework
  structure factors (`FrameworkEik`) with the same k-vector packing, gating
  (`ksqr>1e-10 && ksqr<ReciprocalCutOff`, LAMMPS-Ewald override) and prefactors as
  `TotalFourierEwald_CalculateEnergy`. RASPA2 only adds this at run time — it is never in
  RASPA2's grid file.
- `E_total = E_real + E_recip` is the complete Ewald electrostatic energy of the probe
  and is **independent of the Ewald splitting parameter α** — this is the built-in
  correctness check (see Validation).

**`*.cube`** — Gaussian cube files of each grid (values in K) for visualization
(VMD/OVITO). Capped points hold `2 × EnergyGridOverlapCriteria` (RASPA2 MakeGrid storage
parity). Voxel vectors are diagonal (`Delta`), origin at `Shift`, geometry in Bohr.

### Implementation map

| piece | file | symbol |
|---|---|---|
| struct members | `src_clean/data_struct.h` (Components) | `ComputeEnergyGrid`, `EnergyGridSpacingVDW/Coulomb`, `EnergyGridOverlapCriteria`, `EnergyGridTypes` |
| keyword parsing + banner | `src_clean/read_data.cpp` | after the DensityGrid keyword block |
| GPU kernels | `src_clean/VDW_Coulomb.cu` (end of file) | `EnergyGrid_Framework_VDW`, `EnergyGrid_Framework_CoulReal`, `EnergyGrid_Framework_Fourier` (block-per-grid-point, threads stride framework atoms / k-vectors, shared-memory reduction) |
| host driver + writers | `src_clean/VDW_Coulomb.cu` | `Compute_EnergyGrid` (decl. in `VDW_Coulomb.cuh`), `EnergyGrid_Dimension`, `EnergyGrid_WriteCube` |
| hook | `src_clean/main.cpp` | right after `Check_Simulation_Energy(..., INITIAL, ...)` — framework atoms and Ewald structure factors are on the device and the box holds no adsorbates yet (pure framework field, same state RASPA2 MakeGrid uses) |
| example | `Examples/EnergyGrid/` | MFI + CO2 probes |

The kernels reuse the production energy functions verbatim: `VDW()` / `CoulombReal()` /
`PBC()` from `maths.cuh` (which is why the code lives in `VDW_Coulomb.cu` — the only
translation unit that includes `maths.cuh`; `PBC` is a non-static symbol, so including
`maths.cuh` elsewhere would violate the one-definition rule). VDW gradients are analytic
(`DF = (1/r)dU/dr`, RASPA2 `potentials.c` LJ formula; 12-6-4 path handled).

### Validation record (Quest A100, MFI-2x2x2-P1 + CO2 probes, spacing 1.0 Å)

Artifacts: `/gpfs/projects/p32082/gRASPA_Debug/energy_grid_ab/` (RASPA2 reference run,
`compare_grids.py` binary-.grid parser, job scripts).

- vs **RASPA2 MakeGrid binary grids** (same cif, FF, cutoffs, EwaldPrecision): VDW grids
  (C_co2, O_co2) max relative deviation ≈ 5e-6 = RASPA2's float32 storage precision;
  **zero** cap-set mismatches. Coulomb `E_real` max relative deviation 4.8e-6.
- **α-independence** of `E_total`: EwaldPrecision 1e-6 → 1e-9 moves `E_real` by up to
  8165 K per point while `E_total` changes ≤ 2.4 K (k-vectors 4140 → 17575).
- Surgical when OFF (see above).
- RASPA2 gotcha for re-running the reference: `simulate` writes grids under
  `$RASPA_DIR/share/raspa/grids/...` — export `RASPA_DIR=$PWD` and pre-create
  `share/raspa/`, or it segfaults on a NULL `fopen`.

---

## 2. Widom energy grid (`ComputeWidomEnergyGrid`)

### What is recorded

Every **production** Widom sample (the `WIDOM` branch in `RunMoves`, `src_clean/axpy.cu`):
the selected trial's **first-bead position** bins the sample into a voxel grid over the
simulation cell (fractional coordinates, density-grid convention), accumulating

- `W` — the CBMC Rosenbluth weight, **including** the Ewald reciprocal and tail-correction
  factors (`Insertion_Body` multiplies `W` by `exp(-β·ΔE_Ewald)` and `exp(-β·ΔE_tail)`),
- `E_total` — the full insertion energy `MoveEnergy` (HG/GG VDW + real Coulomb +
  `GPU_EwaldDifference_General` Ewald reciprocal+exclusion + tail),
- `E_ewald` — the Ewald part alone (`HGEwaldE + GGEwaldE`).

Failed CBMC constructions count toward the global totals with `W = 0` (exactly matching
`RecordRosen`, so the feature's `<W>` is the same statistic gRASPA prints) but cannot be
binned spatially — a failed multi-trial CBMC sample has no single defined position. Their
count is reported separately.

### Input keywords

| keyword | default | meaning |
|---|---|---|
| `ComputeWidomEnergyGrid yes` | off | enable |
| `WidomEnergyGridPoints nx ny nz` | 50 50 50 | voxels spanning the simulation cell |
| `WidomEnergyGridRawDump yes` | off | also write every sample as one line |

### Output files (written at end of production)

**`WidomEnergyGrid_<comp>_<name>.grid`** — one line per voxel (voxel centers, Cartesian):
```
x  y  z  count  meanW  betaF  E_W  E_ewald_W  minE
```
- `meanW = Σ W / count` — local mean Rosenbluth weight
- `betaF = -ln(meanW)` — free energy in kT units; this is RASPA2's `FreeEnergyProfile`
  quantity (RASPA2 only has it in 1D; there is **no** 3D Widom grid in RASPA2 — audited)
- `E_W = Σ(E_total·W) / Σ W` — Rosenbluth-weighted mean insertion energy (K); the local
  version of gRASPA's printed "AVERAGE WIDOM ADSORPTION ENERGY"
- `E_ewald_W` — same weighting, Ewald part only (K); exactly 0 for chargeless probes
- unsampled voxels: `count = 0`, remaining columns `?`

**`*_EnergyW.cube`, `*_BetaF.cube`, `*_EwaldW.cube`** — the three maps for visualization.
Unsampled voxels: `2.4054e7` (K maps) / `50.0` (betaF map).

**`WidomSamples_<comp>_<name>.dat`** (RawDump) — one line per successful sample:
```
x  y  z  W  E_total  E_vdw  E_real  E_ewald  E_tail      (Å, -, K...)
```
— literally the Widom insertion energy at each sampled point.

A summary line in the main output gives the global cross-check numbers (see below).

### Statistics: how this relates to gRASPA's printed Widom numbers

- The samples are the **same** ones entering the printed statistics. The feature's global
  `<W>` (ΣW over all attempts / N attempts) equals the printed
  "Averaged Rosenbluth Weight" (verified to 3e-10).
- The printed "AVERAGE WIDOM ADSORPTION ENERGY" is a **block-averaged** estimator:
  mean over blocks of `Σ_b(E·W) / N_b / <W>_b` (`print_statistics.cuh`, `OverallAvgWr`
  is the per-block `<W>`). The feature's `E_W` is the **pooled** estimator `Σ(E·W)/ΣW`.
  They agree in the limit of small block-to-block scatter; with strong scatter they can
  differ by ~0.1% (verified: rebuilding the block estimator from the raw dump reproduces
  the printed value to 0.003 K). Neither is wrong — they are different estimators.
- Henry coefficient per voxel: `KH(r) ∝ meanW(r) / (R·T·ρ_framework·W_idealgas)` — the
  same conversion RASPA2/gRASPA apply to the global `<W>` (`ComputeFugacity`,
  `print_statistics.cuh`).

### Implementation map

| piece | file | symbol |
|---|---|---|
| struct members | `src_clean/data_struct.h` (Components) | `ComputeWidomEnergyGrid`, `WidomEnergyGridPoints`, `WidomEnergyGridRawDump`, `WidomGrid*` accumulators |
| keyword parsing + banner | `src_clean/read_data.cpp` | after the EnergyGrid keyword block |
| sample capture | `src_clean/axpy.cu` (WIDOM branch of `RunMoves`) | after `RecordRosen`; one 24-byte `cudaMemcpy` per sample: `Sims.New.pos[selectedTrial]` (single bead) or `Sims.Old.pos[0]` (chains) — the same device layout `Insertion_Body`'s block-pocket check reads |
| accumulation + writers | `src_clean/write_data.h` | `Accumulate_WidomEnergyGrid`, `WriteWidomEnergyGrid` |
| flush | `src_clean/axpy.cu` (`MCEndOfPhaseSummary`) | next to the density-grid flush |
| example | `Examples/WidomEnergyGrid/` | CO2 in MFI |

### Validation record (Quest A100, CO2 in MgMOF74_charged 5x3x3, 50k samples)

Artifacts: `/gpfs/projects/p32082/gRASPA_Debug/energy_grid_ab/run_WidomGrid/`, `wg_job.sh`.

- global `<W>` vs printed "Averaged Rosenbluth Weight": relative difference 2.8e-10.
- pooled `E_W` = −2826.21 K; block estimator rebuilt from the raw dump reproduces the
  printed −2828.16 K to 0.003 K (difference = estimator shape, see above).
- Ewald column nonzero in 100% of samples (charged probe).
- voxel counts sum exactly to the successful-sample count; `betaF ≡ −ln(meanW)` per voxel.
- Surgical when OFF (same 3-case sweep).

---

## 3. Notes for AI agents (and future maintainers)

1. **Surgical-OFF invariant.** Both features must remain byte-identical when their
   keyword is absent. Banners print **only when a feature is ON**. Any change here must
   re-run the pre-merge gate: build base and feature binaries, run
   `gRASPA-Debug/sweep_compare.sh BASE_RUNS FEAT_RUNS empty_manifest.txt` over standard
   Examples (see `/gpfs/projects/p32082/gRASPA_Debug/energy_grid_ab/eg_surgical_job.sh`
   for a template). Verdict must be CONSTANT everywhere.
2. **Parser gotcha.** `read_data.cpp` matches keywords by *substring* (`str.find`), then
   guards with an exact `termsScannedLined[0] ==` compare. When adding keywords, check
   both directions for substring collisions (e.g. `ComputeEnergyGrid` vs
   `ComputeWidomEnergyGrid` — safe only because "Widom" separates the substrings and the
   exact-token guard exists). Keep the guard.
3. **Banner braces.** A brace-less `if` around a banner once swallowed the next statement
   (density-grid bug, fixed in bdcb953). Keep braces around all banner prints.
4. **Units discipline.** Everything internal is 10 J/mol. Convert to K only at the
   output boundary (×1.2027242847). `FF.CutOffVDW` / `FF.CutOffCoul` are stored
   **squared**. `Beta = 1/(0.831446·T)`.
5. **Where device data lives.** Framework atoms: `Sims.d_a[comp<NComponents.y]`.
   Ewald: `Sims.Box.{Alpha,Prefactor,kmax,ReciprocalCutOff,Volume}` (host-readable
   values), `Sims.Box.FrameworkEik` (device pointer, populated at the INITIAL energy
   check when `!noCharges`). Host mirrors: `Vars.Box[i]` (host `Cell`/`InverseCell`).
6. **Widom trial-position layout** (needed if you extend the Widom grid): single bead →
   `Sims.New.pos[CBMC.selectedTrial]`; chain → first bead `Sims.Old.pos[0]`, atom i ≥ 1
   at `Sims.New.pos[selectedTrialOrientation·(molsize−1)+(i−1)]`.
7. **Validation recipes that worked** (reuse them):
   - RASPA2 binary `.grid` layout parser: `energy_grid_ab/compare_grids.py` (header:
     double, 3×int32, 4×3 doubles, 3×int32 [+ double EwaldPrecision for Coulomb], then
     8·(Nx+1)(Ny+1)(Nz+1) float32, slot-0-first; slot 0 = energy, internal units).
   - α-independence: run twice with different `EwaldPrecision`; `E_total` must be
     invariant while `E_real`/`E_recip` shift.
   - Widom cross-check: the feature's summary line vs the printed Rosenbluth average
     (must match ~1e-10); pooled-vs-block estimator reconciliation from the raw dump.
8. **Known issues** are tracked in `KNOWN_ISSUES.md` — read it before touching the Ewald
   difference code or single-bead charged species.

### Caveats (by design, documented not fixed)

- Serial multi-box mode (`SingleSimulation yes` with several systems run one-by-one):
  `MCEndOfPhaseSummary` runs after each box and loops all sims, so earlier boxes' grids
  are re-written (byte-identical) and their summary line prints again. Cosmetic; same
  behavior as the density grid.
- `UseDNNforHostGuest` runs: `DNN_Replace_Energy()` folds host–guest terms into `DNN_E`,
  so the Widom raw dump's `E_vdw/E_real/E_ewald` columns no longer sum to `E_total`
  (which stays correct) and `E_ewald` loses its HG part. Non-DNN runs are unaffected.
- Widom-grid voxel statistics are over successfully constructed samples only (failed
  samples have W=0 and no defined position); the global `<W>` includes them, matching
  the printed statistic.
