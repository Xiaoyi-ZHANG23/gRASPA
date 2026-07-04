# Known issues

## 1. UPSTREAM BUG (pre-existing, NOT fixed on this branch): stale position/charge in the Ewald difference for single-bead charged species in CBMC insertion / Widom

**Status:** documented only. Discovered 2026-07-04 during adversarial review of the
Widom energy-grid feature; verified by direct code reading (two independent passes).
It exists on upstream `main` and predates every commit on this branch. Planned to be
fixed separately.

### Mechanism

`GPU_EwaldDifference_General` (`src_clean/Ewald_Energy_Functions.h`) stages the inserted
molecule's atoms into `Sims.Old` before building the wave vectors:

- `Initialize_Copy_Positions_Together` (same file, `MoveType == INSERTION` branch) copies
  trial-orientation atoms with `if(ij < chainsize) Old.pos[ij+1] = New.pos[Location*chainsize+ij]`
  — for a **single-bead** molecule `chainsize = molsize - 1 = 0`, so **nothing is copied**.
- The first bead slot `Old.pos[0]` / `Old.charge[0]` is written only by
  `get_random_trial_orientation` (`src_clean/mc_widom.h`, `Mol.pos[0] = NewMol.pos[FirstBeadTrial]`),
  which runs **only for molsize > 1** (chain growth).
- `Initialize_WaveVector_General` then reads `pos = Old.pos[ij]` for `ij < Newsize` — for
  a single bead that is `Old.pos[0]`: **stale data from whatever move used the scratch
  buffer last**.

### Impact

- Affects the **Ewald reciprocal difference** (and therefore the Rosenbluth weight and
  acceptance) of CBMC `INSERTION`, `DELETION`-symmetric paths that reuse this staging,
  and `WIDOM` sampling — but **only** for species that are BOTH single-bead
  (`Moleculesize == 1`) AND charged (`hasPartialCharge`), e.g. bare ions (Na+, Cl-),
  charged single-site probes. The gate is
  `if(!FF.noCharges && SystemComponents.hasPartialCharge[comp])` in
  `mc_swap_utilities.h` (`Insertion_Body`).
- **Not affected:** all multi-bead molecules (CO2, water, alkanes — the chain kernel
  writes `Old.pos[0]` correctly); chargeless single beads (CH4, Kr, Xe — Ewald skipped);
  `SINGLE_INSERTION`/`SINGLE_DELETION` moves (different, correct staging branch that
  copies `New → Old` directly); framework-only quantities.
- The Widom energy-grid feature on this branch reads its **positions** from the correct
  location (`Sims.New.pos[selectedTrial]`), but for single-bead charged species the
  recorded `E_ewald` values would inherit this upstream error (they are whatever
  `Insertion_Body` computed — consistent with gRASPA's own statistics, both wrong).

### Suggested fix (for whoever picks this up)

In `Initialize_Copy_Positions_Together`, INSERTION branch, also stage the first bead:

```cpp
else if(MoveType == INSERTION || MoveType == CBCF_INSERTION)
{
  if(ij == 0) // first bead: selected first-bead trial lives in New (Location for
  {           // single beads = CBMC.selectedTrial; for chains Old.pos[0] is already
              // set by get_random_trial_orientation, rewriting it is harmless ONLY
              // if the correct first-bead index is passed here — needs care!
    ...
  }
  if(ij < chainsize) { Old.pos[ij+1] = New.pos[Location*chainsize+ij]; ... }
}
```

Care needed: for chains, `Location` at this point is `selectedTrialOrientation`, not the
first-bead trial index, so the fix must special-case `chainsize == 0` (use `Location` as
the first-bead index — `Insertion_Body` passes `SelectedTrial = CBMC.selectedTrial` for
single beads) rather than blindly writing `Old.pos[0]` for all sizes. After fixing, run:
1. an A/B against a CPU/analytic Ewald difference for a single charged bead insertion;
2. the INITIAL vs FINAL energy-drift check on an ion-containing example
   (e.g. `Examples/Ionic-MOF-mixtures`);
3. the full surgical sweep (`sweep_compare.sh`) — this fix SHOULD change results for
   single-bead-charged cases (list them in the expected-diff manifest) and change
   nothing else.

### Marker in code

The exact spot carries a `//KNOWN ISSUE (see KNOWN_ISSUES.md)//` comment in
`Ewald_Energy_Functions.h`.

---

## 2. Minor observations (not bugs on this branch)

- **`NumberOfBlocks 1` did not take effect** in a Henry-coefficient run: the run used the
  default `Nblock = 5` (`data_struct.h`), evidenced by nonzero ± error bars on the
  printed Widom averages. Parser behavior worth checking when block counts matter.
- Printed "AVERAGE WIDOM ADSORPTION ENERGY" is a block-averaged estimator
  (mean over blocks of `Σ(E·W)/N_b/<W>_b`); the Widom energy grid's `E_W` is the pooled
  `Σ(E·W)/ΣW`. They differ by O(block scatter²) — both are legitimate; do not "fix" one
  to match the other without understanding this (see `ENERGY_GRIDS.md` §2).
