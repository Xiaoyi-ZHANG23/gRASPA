# Debugging task — CO₂ in NaX zeolite adsorbs nothing

> Give this file (plus `input/` and `challenge/`) to the agent under test. **Do not** give it the
> `solution/` directory — that's the answer key.

## The report
We simulate **CO₂ adsorption in the Na-exchanged faujasite zeolite NaX** (303 K, 10 kPa). The 55
extra-framework **Na⁺ cations** are a mobile framework component, and the important pairwise
Lennard-Jones interactions — especially the strong **CO₂–Na⁺** pairs — are supplied as explicit
`# mixing rules to overwrite` entries in `force_field.def`, layered on top of the Lorentz-Berthelot
base rules in `force_field_mixing_rules.def`.

**Symptom:** with the current gRASPA build this case adsorbs essentially **0 CO₂** (total energy
stays huge/positive). It *should* adsorb on the order of **~17–20 CO₂** at equilibrium — that's what
an older gRASPA build (before the recent 12-6-4 / `UseLJ1264` force-field refactor) produced, and it
matches RASPA2. It looks as if the pairwise `force_field.def` overrides — including CO₂–Na⁺ — are
**not being applied** (the run behaves as if it fell back to plain Lorentz-Berthelot mixing).

**Curiously:** our **MOF** examples are unaffected — their loadings are unchanged. The bug seems
specific to inputs like this zeolite. The override-parsing code lives in
`OverWrite_Mixing_Rule()` in `src_clean/read_data.cpp`.

## Your job
Find the **root cause**, fix it, and prove the fix is correct **and surgical** (it must not change
the behavior of inputs that were already working).

### Tier 1 — required, GPU-free (minutes)
A standalone, GPU-free reproduction of the relevant parser is provided in `challenge/`. It reads a
`force_field.def` and prints how many pairwise overrides it would apply. Build it, run it on
`input/force_field.def`, and observe the behavior. Then **fix `challenge/parser_under_test.cpp`** so
the NaX `force_field.def`'s overrides are applied correctly. See `challenge/README.md` for build/run.

### Tier 2 — optional, needs a GPU + the gRASPA build toolchain
Reproduce end-to-end: build gRASPA at the current `HEAD` and at an earlier commit (use `git log` /
`git bisect` to find where the behavior changed), run `input/` (it already pins `RandomSeed 0`), and
compare loadings/energies with the constant-results gate `debugging/score.py`.

## Suggested method
See [`METHOD.md`](METHOD.md) for the debugging method to apply: the determinism harness
(`RandomSeed 0`), the standalone parser repro, the **specificity test** (change one thing in the
input and see whether the symptom flips), and worktree A/B across commits.

## Deliverables
1. **Root cause** — what the code does wrong and *why* it only bites this kind of input.
2. **Location** — `file:line` (and the commit that introduced it, if you bisect).
3. **The fix** — a diff. Apply it to `challenge/parser_under_test.cpp`.
4. **Verification** — how you confirmed it (commands + observed before/after).
5. **No-regression argument** — evidence the normal / MOF path is unchanged.

## Ground rules
- **Do not special-case this file or its contents.** The fix must be general.
- **Do not fix it by editing the input** (`force_field.def`). Fix the parser.
- A correct fix should never *silently drop* valid input — prefer to apply it, or warn.
