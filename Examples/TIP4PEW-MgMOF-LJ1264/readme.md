# TIP4PEW Water in Mg-MOF with 12-6-4 LJ Potential (GENERIC2_HC)

## Overview

GCMC simulation of TIP4PEW water adsorption in qmof-aa922ec (Mg-MOF) at 298 K / 2000 Pa,
using the 12-6-4 Lennard-Jones potential from Du et al. (JCTC 2020) for Mg-Ow interactions.

The 12-6-4 potential adds a charge-induced dipole term (r^-4) to the standard 12-6 LJ,
capturing polarization effects at open metal sites:

```
U(r) = A/r^12 - B/r^6 - C/r^4
```

## Key Input Settings

- `UseLJ1264 yes` in simulation.input enables polynomial VDW mode
- `force_field.def` contains GENERIC2_HC override for Mg-OwH2O_TIP4PEW pair
- All other pairs use standard Lorentz-Berthelot LJ mixing

## GENERIC2_HC Parameters (from RASPA2 reference run)

```
Mg  OwH2O_TIP4PEW  GENERIC2_HC  0  0  5458  2234  0  -14427061
```

Mapping: `p0*exp(-p1*r) - p2/r^4 - p3/r^6 - p4/r^8 - p5/r^12`
- p2 = 5458 [K A^4] (C_4 charge-induced dipole)
- p3 = 2234 [K A^6] (C_6 dispersion)
- p5 = -14427061 [K A^12] (C_12 repulsion, note: negative sign)

## Reference

- Du, H.; Rodriguez, A.; Lin, L.-C.; Chen, J. *J. Chem. Theory Comput.* 2020, 16, 6060-6072
- RASPA2 implementation: https://github.com/haoyuanchen/RASPA-tools/tree/master/LJ1264Potential

## Comparison: gRASPA vs RASPA2

Both codes run 100k production cycles, 298 K, 2000 Pa, CutOff VDW 12.8 A, truncated, no tail corrections.

| Property | gRASPA | RASPA2 v2.0.47 | Note |
|----------|--------|-----------------|------|
| Loading (molecules) | 473.9 +/- 13.9 | 487.4 +/- 10.4 | Within error bars |
| Loading (mol/kg) | 40.7 +/- 1.2 | 41.8 +/- 0.9 | Within error bars |
| GG VDW (K) | +524,906 | +521,432 | 0.7% diff |
| GG Coulomb (K) | -2,924,905 | -2,932,894 | 0.3% diff |
| Energy drift (K) | 0.0 | ~5.6e-9 | Both excellent |

- Formula verified to machine precision (<10^-13 K) against haoyuanchen potentials.c
- GG VDW close match confirms standard LJ mixing is correct
- Loading overlap within error bars confirms physical consistency
- HG energies differ due to different Ewald implementations and stochastic sampling

## Files

- `simulation.input` - Main configuration (UseLJ1264 yes)
- `force_field.def` - GENERIC2_HC pair override
- `force_field_mixing_rules.def` - Standard LJ parameters (truncated, no tail)
- `pseudo_atoms.def` - Atom types (charges from CIF)
- `qmof-aa922ec.cif` - Framework structure with DDEC charges
- `TIP4PEW.def` - TIP4P/Ew water molecule
- `output.txt` - gRASPA output (100k production cycles)
- `RASPA2-reference-output.txt` - RASPA2 v2.0.47 output (same conditions)
