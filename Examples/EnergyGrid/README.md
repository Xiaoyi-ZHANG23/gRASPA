# EnergyGrid example (RASPA2 `MakeGrid` parity, Ewald Coulomb)

Computes the precomputed framework energy grid for MFI (2x2x2 P1 supercell treated as
one unit cell) once, before any MC cycle, then runs a tiny 10+10-cycle GCMC of CO2 as a
smoke test that the feature does not disturb the simulation.

## Keywords (RASPA2 names, so a RASPA2 MakeGrid input can be reused directly)

| keyword                     | default | meaning                                                        |
|-----------------------------|---------|----------------------------------------------------------------|
| `ComputeEnergyGrid yes`     | off     | enable the one-shot framework energy grid                      |
| `SpacingVDWGrid x`          | 0.15 A  | requested VDW grid spacing (RASPA2 dimensioning: N odd, N+1 planes) |
| `SpacingCoulombGrid x`      | 0.15 A  | requested Coulomb grid spacing                                 |
| `GridTypes name1 name2 ...` | (none)  | pseudo-atom probes, one VDW grid per name (RASPA2 `NumberOfGrids` count is not needed) |
| `EnergyGridOverlapCriteria` | 1.0e7   | cap, internal energy units (RASPA2 `EnergyOverlapCriteria` parity) |

The Coulomb grid (unit +1e probe) is made automatically when charges are on
(`ChargeMethod Ewald` and charged atoms present); it is skipped for `noCharges` runs.

## Output (under `EnergyGrid/System_0/`)

- `EnergyGrid_VDW_<name>.grid` — ASCII, RASPA2 `MakeASCIGrid` parity: 7 columns
  `x y z E Fx Fy Fz` (E in Kelvin, F = -gradient in K/A), `?` in the four value columns
  where E >= EnergyOverlapCriteria. Loop order i(x) outer, j(y), k(z) inner; N+1 planes
  per axis (last plane = periodic image of the first).
- `EnergyGrid_Coulomb_Ewald.grid` — ASCII, 6 columns `x y z E_total E_real E_recip`
  (Kelvin, for a unit +1e probe charge).
  - `E_real` is exactly what RASPA2 `MakeGrid` tabulates in its Coulomb grid
    (`CalculateDerivativesAtPositionReal`: sum of Prefactor*q*erfc(alpha*r)/r).
  - `E_recip` is the framework<->probe reciprocal-space cross term (computed from the
    framework structure factors) that RASPA2 only adds at run time; `E_total =
    E_real + E_recip` is the complete, alpha-independent Ewald electrostatic energy.
  - cap rule is RASPA2's: `?` where E_real > criteria OR the nearest framework atom is
    closer than 1 A.
- `*.cube` — Gaussian cube files (values in Kelvin) of each grid for visualization;
  capped points hold 2*EnergyOverlapCriteria (RASPA2 MakeGrid storage parity).

## Pipeline note

`awk '{print $4}' EnergyGrid_VDW_<name>.grid` gives the single energy column (with `?`
markers) in the same shape the PSED-CNN ASCII pipeline consumes; for electrostatics,
column 4 of `EnergyGrid_Coulomb_Ewald.grid` is the total Ewald energy of a +1e probe.
