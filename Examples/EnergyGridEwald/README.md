# EnergyGridEwald — Ewald electrostatic energy grid on a DDEC-charged MOF

Companion to `Examples/EnergyGrid/` (which documents the full feature on silicalite).
This example exercises the part that matters for charged adsorbates: the **complete Ewald
electrostatic grid** of a real, DFT-charged MOF — and shows how to verify it.

- **Framework:** `qmof-d475ee7` from the QMOF database (Rosen et al., *Matter* 2021,
  DOI 10.1016/j.matt.2021.02.015; data CC-BY 4.0), DFT-DDEC partial charges written in
  `_atom_site_charge`. Single unit cell (~36.4 x 36.8 x 36.8 A), 244 atoms — large enough
  that `UnitCells 1 1 1` satisfies the 12.8 A cutoff with no replication.
- **Grids computed once before MC** (`ComputeEnergyGrid yes`, spacing 1.0 A here):
  two VDW probe grids (`GridTypes C_co2 O_co2`) and, because `ChargeMethod Ewald` with a
  charged framework, the Coulomb grid `EnergyGrid_Coulomb_Ewald.grid` for a unit +1e probe:
  6 columns `x y z E_total E_real E_recip` (Kelvin).
- Then a tiny 10+10-cycle CO2 GCMC runs as a smoke test that the feature leaves the
  simulation untouched.

## What the three Coulomb columns mean

Ewald splits the electrostatic energy into a short-range real-space part and a smooth
reciprocal-space part, weighted by the splitting parameter alpha:

- `E_real`  — sum over framework atoms of `Prefactor * q * erfc(alpha r)/r` (this is the
  quantity RASPA2's `MakeGrid` tabulates in its Coulomb grid).
- `E_recip` — the framework<->probe reciprocal-space cross term, computed from the
  framework structure factors (RASPA2 adds this only at run time; here it is stored).
- `E_total = E_real + E_recip` — the complete electrostatic energy of the probe. Only this
  sum is physical: the real/reciprocal split depends on alpha, the total must not.

## The validation this example is for: alpha-independence

`EwaldPrecision` sets alpha, so changing it redistributes energy between the two parts
while the sum stays fixed. Run the example twice (`EwaldPrecision 1e-6`, then `1e-8`) and
compare the Coulomb grid point by point. Measured on this system (53,971 uncapped points):

| quantity | median shift | max shift |
|---|---|---|
| `E_real`  | 7.5 K | 1251.1 K |
| `E_recip` | 7.5 K | 1251.2 K |
| **`E_total`** | **0.0095 K** | **0.100 K** |

The split moves by up to ~1250 K per point; the total moves by less than 0.1 K over a
~39,000 K energy range — i.e. the stored reciprocal complement is correct. (Grid
computation itself takes well under a second here.)

Quick check after your own two runs:

    paste run1/EnergyGrid/System_0/EnergyGrid_Coulomb_Ewald.grid \
          run2/EnergyGrid/System_0/EnergyGrid_Coulomb_Ewald.grid \
      | awk '$4!="?" && $10!="?" {d=$4-$10; if(d<0)d=-d; if(d>m)m=d} END{print "max |dE_total| =", m, "K"}'

## Why this matters for the density-field workflow

For a charged adsorbate (CO2 here), the full electrostatic energy of a probe at every
voxel — not just the real-space part — is what makes an energy grid usable as a
model input or Boltzmann-prior source: `E_total` is alpha-independent and therefore a
well-defined physical field, co-registered with the same framework's density grids.

Force-field files are the general-purpose set (all-element UFF-style mixing rules +
TraPPE CO2 with C +0.6512 / O -0.3256), so any `_atom_site_charge`-bearing P1 CIF can be
dropped in by changing `FrameworkName`.
