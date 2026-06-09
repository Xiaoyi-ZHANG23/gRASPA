# SOLUTION — answer key & grading rubric (withhold from the agent under test)

## Root cause (the headline "Bug #4")
In `OverWrite_Mixing_Rule()` (`src_clean/read_data.cpp`, first pass — ~line 960 at the buggy commit
`186e4d3`; the guarded fix lands at lines 965/972 at HEAD `9e94189`), the scan that locates the
override section reassigns the marker line on **every** matching line:

```cpp
if (str.find("mixing rules to overwrite", 0) != string::npos) mixrule_startline = counter; // BUG: last wins
```

So when a `force_field.def` contains **two** `# mixing rules to overwrite` markers, the **last** one
wins. The NaX zeolite force field (`input/force_field.def`) is:

```
# mixing rules to overwrite     <- real block: 11 pairwise overrides (incl. CO2-Na+)
11
# type type2 interaction
...11 lines...
# mixing rules to overwrite     <- trailing duplicate marker
0
```

The trailing marker sets `Nmixrule = 0`, the function early-returns
(`if (Ndefint==0 && Nmixrule==0) return;`), and **all 11 overrides are dropped** → silent fallback to
Lorentz-Berthelot. The dropped pairs include the strong **CO₂–Na⁺** interactions
(`O_co2 Na = 200.831`, `C_co2 Na = 362.292`), so CO₂ barely binds → ~0 adsorbed.

**Why MOFs are fine:** their `force_field.def` has a single marker (or a single trailing one), so
last == first; nothing is dropped. This is why the symptom is input-structure-specific.

**Introduced by:** commit `186e4d3` ("Add 12-6-4 LJ / UseLJ1264"), the refactor right after the
Jan-4 merge `3a88f7f`. **Fixed by:** `ea5a825` (merged to this fork's `main`, PR #81).

## The fix (general, one line per marker)
First occurrence wins — don't let a later marker overwrite the first:

```cpp
if (str.find("mixing rules to overwrite", 0) != string::npos && mixrule_startline == 0) mixrule_startline = counter;
// and symmetrically for the defint marker:
if (str.find("number of defined interactions", 0) != string::npos && defint_startline == 0) defint_startline = counter;
```

Reference implementation: `parser_fixed.cpp` (diff vs `../challenge/parser_under_test.cpp` is exactly
the two `&& *_startline == 0` guards). In the real codebase the same idiom is applied in
`read_data.cpp`; the broader fix also **warns** rather than silently dropping unrecognized override
lines (lesson: *never silently drop*).

## Expected results

### Tier 1 — GPU-free harness
| input | buggy `parser_under_test.cpp` | fixed |
|---|---|---|
| `input/force_field.def` (trigger) | `Nmixrule=0`, early-return, 0 overrides | `Nmixrule=11`, incl. `C_co2 Na 362.292` and `O_co2 Na 200.831` |
| `ff_control_mof.def` | `Nmixrule=3` | `Nmixrule=3` (unchanged) |
| `ff_control_zero.def` | `Nmixrule=0` | `Nmixrule=0` (unchanged) |
| `grade_control_adv.def` (held-out; agent never sees it) | `Nmixrule=2` | `Nmixrule=0` (first marker wins) |

Grade with: `bash grade.sh` — 5 checks: trigger→11 with CO₂–Na present, both controls unchanged, and
the **held-out adversarial** control (a *leading* zero-count marker before a real block) →0. The
adversarial check fails non-general fixes such as "ignore a marker whose count is 0" (which yields 2),
while passing both legitimate general fixes (`&& *_startline==0` and the explicit found-flag variant).

**Whole-suite no-regression (built into `grade.sh`, also runnable alone as
`bash regression_examples.sh <candidate>`):** parses **every** `Examples/*/force_field.def` (34 files)
with the candidate vs the pristine buggy baseline and confirms the fix changes **only**
`CO2_NaX_Zeolite` (0→11) — the other **33 are byte-identical**. Of all examples, only
`CO2_NaX_Zeolite` has the 2-marker structure, so a correct fix is provably surgical across the suite.
(Tier 2 / GPU equivalent: run the full suite head-vs-fix through `debugging/score.py`.)

### Tier 2 — full GPU run
- buggy build: CO₂ loading ≈ **0**, total energy stays huge/positive.
- fixed build (and any pre-`186e4d3` build): CO₂ adsorbs (init reaches **17**, production fluctuates
  ~18–24); final `reference_output.txt` shows `PseudoAtom Type: C_co2 #: 20`, `O_co2 #: 40`
  (i.e. 20 CO₂), total energy ~ −2.9–3.0e5.
- With `RandomSeed 0`, the fixed build vs a pre-bug reference is **bit-comparable**:
  `python3 ../score.py <fixed_output> reference_output.txt` shows **no comparison-class failures**
  (the only `score.py` complaint is the run's inherent absolute `vdw_hh` drift ≈1.65, which is present
  in the reference too — see `../DEBUGGING.md` on comparison-class vs absolute-class checks).

## What the graspa-debug method should surface (process check)
A strong solution uses the skill, not luck:
- **specificity test:** deleting only the trailing `# mixing rules to overwrite` / `0` from
  `force_field.def` makes the buggy harness/build behave correctly → proves that marker is the sole
  trigger.
- **worktree A/B (Tier 2):** `3a88f7f` (good) vs `186e4d3` (bug) vs `HEAD` (fixed), same `RandomSeed
  0`, compared with `score.py`.
- **read the parser:** finds the last-wins marker scan and the `Nmixrule==0` early return.

## Grading rubric
- **Full credit:** identifies last-marker-wins → trailing-`0` → early-return → overrides dropped;
  fixes with the first-occurrence-wins idiom (general); `grade.sh` PASS (trigger→11 with CO₂–Na,
  controls + held-out adversarial unchanged, and the **all-Examples no-regression** clean — only
  CO₂_NaX_Zeolite changes). Bonus: names the introducing/​fixing commits via bisect; notes the
  "never silently drop / warn" principle.
- **Partial:** fixes the trigger but without explaining *why* MOFs were spared, or with no
  specificity/A-B evidence, or via a narrow approach that only happens to pass the held-out
  adversarial control. (Note: the common "ignore a marker whose count is 0" hack is now **auto-caught**
  — it fails the `grade_control_adv.def` check, exit 1.)
- **Fail / red flags:** edits `force_field.def` instead of the code; special-cases this filename or
  these atom types; "fixes" by suppressing the early-return so garbage is read; `grade.sh` still
  FAIL. (Editing the input is auto-caught: `grade.sh` uses its own `grade_*.def` copies.)
