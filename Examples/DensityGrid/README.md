# DensityGrid example — 3D spatial adsorbate-density grid

Regression/demonstration case for the `ComputeDensityProfile3DVTKGrid`-parity feature
(`ComputeDensityGrid`). CO2 GCMC in rigid silicalite (MFI, 2x2x2 P1) at 298 K, 1e4 Pa —
a real rigid framework, where the spatial density "inside the MOF" is meaningful.

## What it exercises
- `ComputeDensityGrid yes` turns the feature on (default is `no`).
- `DensityGridPoints 40 40 40` — 40x40x40 voxels over the unit cell.
- `DensityGridSampleEvery 1` — bin every production cycle (so the integrated grid
  is directly comparable to the run's reported average loading).
- `DensityGridWriteEvery 100000` — large, so the grid is written once at end of production.
- `RandomSeed 0` — deterministic, so the run is bit-comparable across builds.

## Output
At the end of the production phase gRASPA writes, per adsorbate component, under
`DensityGrid/System_0/`:
- `DensityGrid_1_CO2.cube` — Gaussian cube (triclinic-correct voxel vectors, geometry in Bohr).
- `DensityGrid_1_CO2.vtk`  — VTK `STRUCTURED_POINTS` (ASCII; assumes an orthorhombic cell).

The scalar field is the **average number density** of adsorbate atoms,
`count / (Nsnapshots * voxelVolume)` in molecules/Angstrom^3. Both files open in
ParaView/VisIt.

## Validation
Integrating the density over the cell recovers the average number of adsorbate atoms:
`sum(voxel) * voxelVolume ~= Nsnapshots-averaged atom count`, and dividing by the 3 sites
of CO2 reproduces the "# MOLECULES" block average printed in the output. See
`density_grid_ab/validate_grid.py`.

A vetted reference `output.txt` (build noted at the top of the file) is shipped for the
A/B `score.py` gate.

## Reference build

The shipped `output.txt` was produced by the vanilla (non-ML) build at commit `bdcb953`
(branch `feature/density-grid-3d`, off `main` `3fc256d`) on a Quest A100. Validated:
integral over the cell = 57.197 atoms = 19.06567 CO2 molecules, matching the reported
`# MOLECULES` block average to 0.00%. The first line of `output.txt` is the absolute
exe-path banner (machine-specific) — drop it when byte-matching the shipped format.
