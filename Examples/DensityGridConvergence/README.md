# DensityGridConvergence — converged density-field calculation (optional, default OFF)

CO2 GCMC in rigid silicalite (MFI, 2x2x2 P1) at 298 K / 1e4 Pa, the same physical case as
`Examples/DensityGrid/`, but the run **stops production automatically once the 3D density
field itself is converged** instead of running a fixed number of cycles.

## Why

A converged *uptake* is not a converged *density field*: the scalar loading (a sum over the
grid) settles long before every occupied voxel has good statistics. If the density field is
the quantity you care about, cycle counts sized on uptake systematically under-sample it —
and the right number of cycles varies strongly between systems (roughly linearly in the
loading). This feature lets each run decide its own length.

## Keywords (all optional; **defaults reproduce the original behavior byte-for-byte**)

| keyword | default | meaning |
|---|---|---|
| `DensityGridTanimotoStop yes` | `no` | end PRODUCTION early once the density field is converged |
| `DensityGridTanimotoThreshold 0.85` | `0.85` | convergence gate on the split-half Tanimoto similarity |
| `DensityGridMolecularCenter yes` | `no` | bin one point per molecule (its geometric center = the C atom for CO2) so the grid integrates to the **molecule** count |

The check runs every `DensityGridWriteEvery` cycles (here 5000), together with the periodic
grid flush.

## How the convergence test works

The accumulated snapshots are split into two halves by alternating `DensityGridWriteEvery`-
sized blocks (block-interleaving, so each half spans the whole production phase). With
**a** and **b** the two half-grids, the continuous Tanimoto similarity

    T(a,b) = <a,b> / ( <a,a> + <b,b> - <a,b> )

is 1 only when the halves agree voxel-for-voxel; 1 - T is (to first order) the fractional
noise power of the field, and it decays as 1/N with sampling. The run prints T at every
check and stops production once `T >= threshold` for **two consecutive checks**, e.g.:

    Density-grid split-half Tanimoto = 0.8612 (threshold 0.850, consecutive pass 2/2, 35000 snapshots)
    DENSITY GRID CONVERGED (split-half Tanimoto >= 0.850 twice) -- ending production early.

The final grid (`DensityGrid/System_0/*.cube` / `*.vtk`) contains every sample drawn up to
the stop and is normalized by the actual snapshot count, so it is the converged product.

## Two caveats (read before relying on it)

1. **Split-half is an internal, slightly optimistic estimate.** The two halves share one
   trajectory, so residual autocorrelation makes them agree a little sooner than two fully
   independent runs would. For publication-grade convergence, run two replicates with
   different `RandomSeed` values and compare their grids with the same Tanimoto formula;
   use this in-run stop as the cheap, single-process version of that test.
2. **With an early stop, read the loading from the grid, not the block averages.** The
   BLOCK AVERAGES printed at the end divide by the full planned cycle count and will
   under-report if the run stopped early. The density-grid integral
   (`sum(voxel) * voxelVolume`, from the `.cube`) is the correct average loading over the
   sampled portion — with `DensityGridMolecularCenter yes` it equals the average molecule
   count directly.

## Files

Identical inputs to `Examples/DensityGrid/` (same CIF, force field, CO2 model); only
`simulation.input` differs: a large production cap (100k cycles), `DensityGridWriteEvery
5000` as the check cadence, and the three keywords above. With `RandomSeed 0` the run is
deterministic and typically converges well before the cap.
