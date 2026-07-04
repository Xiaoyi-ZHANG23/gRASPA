# WidomEnergyGrid example — spatial map of Widom insertion energies (incl. Ewald)

Pure Widom sampling (`WidomProbability 1.0`, no molecules created) of CO2 in MFI.
Every PRODUCTION Widom sample is recorded: the selected trial's **first-bead position**
bins the sample's full insertion energy (VDW + real Coulomb + **Ewald reciprocal** +
tail, exactly the `Insertion_Body` MoveEnergy that enters the Rosenbluth weight) into a
voxel grid over the simulation cell. This is a 3D generalization of RASPA2's 1D
`ComputeFreeEnergyProfile` (whose stored quantity is betaF(q) = -ln<W(q)>).

## Keywords

| keyword                        | default    | meaning                                             |
|--------------------------------|------------|-----------------------------------------------------|
| `ComputeWidomEnergyGrid yes`   | off        | enable the Widom insertion-energy grid              |
| `WidomEnergyGridPoints nx ny nz` | 50 50 50 | voxels spanning the simulation cell (fractional)    |
| `WidomEnergyGridRawDump yes`   | off        | also write every sample as one line (see below)     |

## Output (under `WidomEnergyGrid/System_<id>/`, written at the end of production)

- `WidomEnergyGrid_<comp>_<name>.grid` — ASCII, one line per voxel (voxel centers):
  `x y z count meanW betaF E_W[K] E_ewald_W[K] minE[K]`
  - `meanW`  = sum(W)/count — local mean Rosenbluth weight
  - `betaF`  = -ln(meanW) — RASPA2 FreeEnergyProfile quantity, in kT units
  - `E_W`    = sum(E_total*W)/sum(W) — Rosenbluth-weighted mean insertion energy
  - `E_ewald_W` = same weighting for the Ewald (reciprocal + exclusion) part alone
  - unsampled voxels: `count = 0`, remaining columns `?`
- `..._EnergyW.cube`, `..._BetaF.cube`, `..._EwaldW.cube` — the three maps for
  visualization (unsampled voxels: 2.4e7 K in the energy cubes, 50 in the betaF cube).
- `WidomSamples_<comp>_<name>.dat` (RawDump) — one line per successful sample:
  `x y z W E_total[K] E_vdw[K] E_real[K] E_ewald[K] E_tail[K]`
  — literally the Widom insertion energy at each sampled point.

## Consistency guarantees (validated on Quest A100, CO2 in MgMOF74_charged, 50k samples)

- The samples are the SAME ones entering gRASPA's printed Widom statistics: the global
  `<W>` in the feature's summary line matches the printed "Averaged Rosenbluth Weight"
  to 3e-10.
- The feature's global `sum(E*W)/sum(W)` is the POOLED version of the printed "AVERAGE
  WIDOM ADSORPTION ENERGY". The printed value is a mean over blocks of per-block
  `sum(E*W)/N_b/<W>_b`; the two estimators agree in the limit of small block scatter
  (rebuilding the block estimator from the raw dump reproduced the printed value to
  0.003 K).
- Failed CBMC constructions count toward `<W>` with W = 0 (matching `RecordRosen`) but
  cannot be binned spatially: a failed multi-trial CBMC sample has no single position.
  Their count is reported in the `.grid` header and the summary line.

## Notes

- Binning uses the first bead (RASPA2's FreeEnergyProfile pins the first bead too, so
  the betaF map is the closest 3D analogue; note RASPA2 PINS the bead per bin while
  Widom sampling scatters it uniformly — voxel statistics are therefore sample-count
  weighted, see the `count` column).
- For chargeless probes (e.g. CH4), `E_ewald_W` is exactly 0 and the Ewald path is
  skipped upstream (`noCharges`/`hasPartialCharge` gate in `Insertion_Body`).
