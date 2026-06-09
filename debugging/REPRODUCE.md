# REPRODUCE.md — reproduce the gRASPA force_field.def bug hunt from a clean clone

This is the **runnable, deterministic** recipe behind `.claude/skills/graspa-debug/SKILL.md`. It
reproduces the headline bug (Bug #4: pairwise `force_field.def` overrides silently dropped) two ways:
a **30-second, GPU-free** parser check, and the **full GPU A/B** with the constant-results gate.

Everything needed is vendored in this repo: the comparator (`debugging/score.py`), the parser repro
(`debugging/repro/`), and the regression input + reference output (`Examples/CO2_NaX_Zeolite/`).

Determinism is the whole point: every comparison uses a **fixed `RandomSeed 0`**, so the same input
on a correct binary yields **bit-comparable** energies and loadings. If a step here is not
reproducible on your machine, that is itself a finding — check seed, build flags, and the `2>&1`
capture before anything else.

## A/B commits (also in ../CLONE_INFO.md / git log)
| Role | SHA | Note |
|------|-----|------|
| good reference | `3a88f7f` | last commit BEFORE the 12-6-4 refactor |
| bug introduced | `186e4d3` | introduced all 4 bugs |
| the fix        | `ea5a825` | parser fix (merged; this branch's HEAD descends from it) |
| regression case added | `cd680e1` | `Examples/CO2_NaX_Zeolite` |

---

## Part 1 — GPU-free parser repro (30 seconds, no toolchain beyond g++)

Proves Bug #4 purely from the parser logic — no GPU, no build of gRASPA.

```bash
cd debugging/repro
g++ -O2 -std=c++17 parse_repro2.cpp -o parse_repro2
./parse_repro2 ../../Examples/CO2_NaX_Zeolite/force_field.def
```

**Expected output:**
- buggy first-pass (last-marker-wins, as in `186e4d3..ea5a825`): **`Nmixrule = 0`** → all 11 pairwise
  overrides dropped → silent fallback to Lorentz-Berthelot.
- fixed first-pass (first-occurrence-wins): **`Nmixrule = 11`** → overrides applied.

Why the trigger works: `Examples/CO2_NaX_Zeolite/force_field.def` has a real overrides block
(11 pairs, lines 3-16) followed by a **trailing** `# mixing rules to overwrite` / `0` (lines 17-18).
The buggy scan keeps reassigning `mixrule_startline` so the LAST marker (count 0) wins. MOF inputs
have a single marker, so they're unaffected — hence "zeolites differ, MOFs don't."

**Specificity check** (prove that trailing marker is the SOLE trigger): delete lines 17-18 from a
copy of the file and re-run — buggy and fixed now both report `Nmixrule = 11`.

---

## Part 2 — Full GPU A/B with the constant-results gate

Run the *same* input on three builds and compare with the vendored gate. Requires the GPU build
toolchain (see `../Cluster-Setup/Quest_April_2026/QUICKSTART.md`).

```bash
# from repo root. $REPO = this checkout.
REPO=$(pwd)

# 1. One worktree per version (keep all three builds alive at once)
git worktree add ../wt_good 3a88f7f      # correct reference
git worktree add ../wt_bug  186e4d3      # buggy
git worktree add ../wt_fix  HEAD         # fixed (this branch)

# 2. Build each. On a Quest LOGIN node you MUST add -tp haswell (else SIGILL/exit 132 on A100).
#    Easiest: build inside a gengpu job (module load works there). See Cluster-Setup/Quest_April_2026/.
for d in ../wt_good ../wt_bug ../wt_fix; do ( cd "$d/src_clean" && bash NVC_COMPILE ); done

# 3. Run the SAME regression input in each, capturing BOTH streams (stderr has the summary!)
for d in wt_good wt_bug wt_fix; do
  run=../$d/run_NaX
  rm -rf "$run"; mkdir -p "$run"
  cp "$REPO"/Examples/CO2_NaX_Zeolite/* "$run"/
  cp ../$d/src_clean/nvc_main.x "$run"/
  ( cd "$run" && ./nvc_main.x > output.txt 2>&1 )   # RandomSeed 0 is already in simulation.input
done

# 4. Compare with the constant-results gate
echo "== good vs fix (expect: NO comparison-class failures — constant/correct) =="
python3 debugging/score.py ../wt_fix/run_NaX/output.txt ../wt_good/run_NaX/output.txt
echo "== good vs bug (expect: comparison-class failures — loadings/energy differ) =="
python3 debugging/score.py ../wt_bug/run_NaX/output.txt ../wt_good/run_NaX/output.txt
```

**Expected:**
- **good vs fix:** no `comparison`-class entries in `failures` (move counts, PseudoAtom counts, and
  final energies match) → results are *constant* → the fix is correct.
- **good vs bug:** comparison-class failures — the buggy build adsorbs **0 CO2** (E ≈ −2.8M) while
  good/fix adsorb **~17 CO2** with energy identical to every digit.

> Reading `score.py`: it reports two classes of check. **comparison** (moves / counts / energies vs
> reference) = "the runs differ." **absolute** (per-component energy drift < 3e-5, structure factor)
> = a property of the single run, and it fires **even on a self-compare** — e.g.
> `score.py Examples/CO2_NaX_Zeolite/output.txt <same file>` exits 1 because that run has an inherent
> `vdw_hh` drift ≈ 1.65, with zero comparison failures. So judge A/B equality by the **absence of
> comparison-class failures**, not by exit code alone.

### Cleanup
```bash
git worktree remove ../wt_good && git worktree remove ../wt_bug && git worktree remove ../wt_fix
```

---

## Definition of done
A change is correct when, on identical `RandomSeed 0` inputs, `debugging/score.py` shows **no
comparison-class failures** between the candidate and a trusted reference across the affected case(s),
AND the full non-MLIP `Examples/` suite shows no comparison-class failures except where the fix was
intended to change behavior. WARN on unrecognized override lines; never silently drop them.
