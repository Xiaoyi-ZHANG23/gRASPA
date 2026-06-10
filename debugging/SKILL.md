---
name: graspa-debug
description: Debug gRASPA energy/correctness discrepancies and force_field.def parser bugs by proving results are CONSTANT (deterministic, bit-comparable) across versions. Drives a git-worktree A/B method with RandomSeed 0 and a vendored correctness gate (debugging/score.py, from AutoJIT-gRASPA). Use when gRASPA results differ between versions or inputs, when forcefield overrides seem silently ignored, when a "zeolites differ but MOFs don't"-type bug is reported, when validating a gRASPA fix against the Examples suite, or as the PRE-MERGE GATE whenever a new feature/keyword/branch is added to gRASPA (sweep_compare.sh + an expected-diff manifest prove the change is surgical).
---

# graspa-debug: Deterministic Bug-Hunting for gRASPA

You are an agent debugging **gRASPA**, an NVIDIA-GPU Monte Carlo simulator for molecular adsorption.
Your job: find WHY two gRASPA runs disagree (different commits, different inputs, or
suspected-wrong output) and prove the root cause **empirically** — by showing results are *constant*
(deterministic and bit-comparable) when they should be, and pinpointing the exact line where they
stop being constant.

Governing principle: gRASPA is **deterministic for a fixed `RandomSeed`** — same seed ⇒ same RNG ⇒
same move sequence ⇒ bit-for-bit identical energies and loadings. Every technique here exploits
that. "Constant results" is both the test oracle and the definition of done.

> **Install (for Claude Code auto-discovery):** copy this file to `.claude/skills/graspa-debug/SKILL.md`
> (this repo's `.gitignore` excludes `.claude/`, so the tracked/uploadable copy lives here in
> `debugging/`). All paths below are written to run **from the repo root**. Reproducible run:
> [`REPRODUCE.md`](REPRODUCE.md) · portable mirror: [`DEBUGGING.md`](DEBUGGING.md).

## When to Use This Skill
- Reported energy/loading mismatch between two gRASPA versions, commits, or branches.
- "Forcefield overrides are being ignored / results fall back to Lorentz-Berthelot."
- A symptom that depends on input *structure* ("zeolites differ but MOFs don't").
- Validating that a fix changes ONLY the intended case and leaves everything else byte-identical.
- Any "the number looks wrong but the run didn't crash" situation (silent-wrong-answer class).

## Prerequisites Check
1. **The repo** — this gRASPA checkout (the debug/fix target). Build entry: `src_clean/NVC_COMPILE`.
2. **GPU build toolchain** — full Quest recipe in-repo at `Cluster-Setup/Quest_April_2026/`
   (`NVC_COMPILE_QUEST_VANILLA`, `compile_graspa.job`, `QUICKSTART.md`); compiler
   `nvhpc/23.3-gcc-10.4.0`.
   ⚠️ **If you compile on the LOGIN node, add `-tp haswell`** or the binary dies with **SIGILL
   (exit 132)** on the older A100 nodes. Energy kernels run on the GPU, so `-tp` never affects
   results. Login-node `module load` silently no-ops — set PATH directly or build inside the SLURM
   job (`gengpu`, account `p32082`).
3. **The constant-results gate** — vendored at `debugging/score.py` (stdlib python3 only).
4. **(optional) AutoJIT-gRASPA** — the JIT sibling (https://github.com/Zhaoli2042/AutoJIT-gRASPA);
   source of `score.py`/`benchmark.sh` and of vetted `examples/*/baseline_output.txt` oracles.

If a piece is missing, say so before proceeding; don't guess.

## The Core Method: Worktree A/B + RandomSeed 0

```
            good ref commit            suspect HEAD              HEAD + candidate fix
            ───────────────            ────────────             ────────────────────
git worktree add wt_good <ref>   wt_head (HEAD)            wt_fix (HEAD, patch applied)
build each (NVC_COMPILE)         build                    build
run SAME input, RandomSeed 0  →  run SAME input        →  run SAME input
        │                              │                          │
        └─────────── debugging/score.py compares all three ──────┘
   wt_good == wt_fix  (constant ⇒ correct)     wt_head differs  (the bug)
```

1. **Pin the seed.** Every input you compare must have `RandomSeed 0` in `simulation.input`.
2. **One worktree per suspect version** — `git worktree add ../wt_<name> <sha>`.
3. **Build each** with `src_clean/NVC_COMPILE` (add `-tp haswell` on login node).
4. **Run identical input**, redirecting **`> output.txt 2>&1`** (see stderr GOTCHA).
5. **Compare with `debugging/score.py`** (next section).

To localize a regression, bisect with worktrees against the good baseline. (That isolated `186e4d3`,
the commit right after the Jan-4 merge `3a88f7f`.) Full runnable recipe: [`REPRODUCE.md`](REPRODUCE.md).

## The Constant-Results Gate: `debugging/score.py`  ★

Don't hand-grep "Total Energy" and eyeball digits — that misses move-count and loading drift.

```bash
python3 debugging/score.py <test_output.txt> <reference_output.txt>
# exit 0 = CONSTANT (identical within tol)   1 = differs   2 = parse error
```

It runs two *classes* of check — read the JSON `failures` array, not just the exit code:

| Class | Check | Meaning |
|---|---|---|
| **comparison (test vs reference)** | Move statistics exact | same RNG seed ⇒ identical move counts |
| **comparison** | PseudoAtom / component counts exact | same final state (e.g. "17 CO2") |
| **comparison** | Final energies within 1e-5 rel / 1e-3 abs | components match the reference |
| **absolute (test file only)** | Per-component energy drift < 3e-5 | run's own running-vs-recomputed self-check |
| **absolute / relative** | GPU drift (vs reference, else absolute) | GPU/CPU consistency not degraded |
| **absolute** | Structure factor CPU/GPU < 1e-4 | Ewald path consistent |

**Comparison-class** failures ⇒ *the two runs genuinely differ* (your A/B signal). **Absolute-class**
failures (energy drift, structure factor) are a property of the *test run itself* and fire **even on
a self-compare**.

⚠️ **Verified gotcha:** `python3 debugging/score.py Examples/CO2_NaX_Zeolite/output.txt <same file>`
returns **exit 1**, not 0 — that run has an inherent `vdw_hh`/`total` drift ≈ 1.65 (absolute-class),
with *zero* comparison-class failures. So a nonzero exit does **not** by itself mean "the runs
differ." For an A/B equality question, the decisive evidence is **no comparison-class entries in
`failures`** (and both runs showing the *same* drift). Reserve "exit 0" as definition-of-done for
drift-clean systems.

(One-shot compile+run+score in AutoJIT's working-dir layout: `AutoJIT-gRASPA/scripts/benchmark.sh`,
which already runs with `2>&1`.)

## stderr GOTCHA (the silent false-IDENTICAL)
gRASPA writes its end-of-run summary — **`PRODUCTION AVERAGE ENERGY`** and the final
**`PseudoAtom Type: ... #:`** loadings — to **stderr**, not stdout. Capture only stdout and you
compare truncated files and wrongly conclude "IDENTICAL." Always **`> output.txt 2>&1`**. `score.py`
parses these lines, so it needs the merged capture.

## Specificity Test (prove the trigger is the SOLE cause)
Delete *only* the suspected trigger from the input and re-run HEAD: if HEAD then matches the good
reference (no comparison-class failures), the trigger is the sole cause and your fix doesn't touch
the normal path. (Bug #4 trigger: a trailing duplicate `# mixing rules to overwrite` / `0` block —
see `Examples/CO2_NaX_Zeolite/force_field.def` lines 17-18.)

## Standalone Parser Repro (debug without a GPU build)
Most parser bugs live in `src_clean/read_data.cpp` (`OverWrite_Mixing_Rule()`). You don't need the
~20-min GPU build to test parsing — vendored templates reproduce Bug #4 with plain `g++` in seconds:

```bash
cd debugging/repro
g++ -O2 -std=c++17 parse_repro2.cpp -o parse_repro2
./parse_repro2 ../../Examples/CO2_NaX_Zeolite/force_field.def
# buggy (last-marker-wins) → Nmixrule=0 ;  fixed → Nmixrule=11
```

(`read_data.cpp` is large and fragile — repro in isolation first.)

## Recurring Failure Mode: "Never Silently Drop"
The entire `186e4d3` bug class shares one signature: **force-field overrides silently discarded →
wrong-but-quiet results.** When auditing parser/override code, grep for:
- **Last-marker-wins instead of first-occurrence** — reassigns `*_startline` on every matching
  marker, so a trailing duplicate (count 0) wins. Fix idiom: `... && *_startline == 0`.
- **Fixed-offset data windows** — `[startline+3, +3+N)` breaks when comments sit between count and
  data. Fix: scan past the count, skip comments/blanks, consume exactly N real entries.
- **Over-tight token-count guards** — `if(tokens.size()==5)` drops a valid 6-token
  `I J lennard-jones-1264 eps sig C4` override. Accept the documented arity.
- **Tail-correction prefactor** — per-pair tail is `2π·∫U r² dr`; a stray `4π` doubles TailE.

Rule when fixing: **WARN on any unrecognized / unconsumed override line; never drop it silently.**

## Regression Suite (prove the fix is surgical)
After a fix, run the whole non-MLIP `Examples/` suite head-vs-fix and `score.py` each pair: expect
**no comparison-class failures anywhere except the case(s) the fix intentionally changes**, which
should now match the good reference. `Examples/CO2_NaX_Zeolite/` is the Bug #4 regression case (ships
`RandomSeed 0`, the trailing-marker `force_field.def`, and a reference `output.txt` — its first line
is the exe-path banner; drop it when matching shipped format).

## Validating a NEW Feature (run this on EVERY feature branch)

The same machinery, run as a pre-merge gate. **Surgical** means: the feature changes the cases it
is *supposed* to change and **nothing else** (no comparison-class failures anywhere else).

1. **Declare intent first.** Write `expected_diffs.txt` — one `Examples/` case name per line for
   every case the feature SHOULD change. Everything not listed must stay comparison-clean. An
   empty manifest = "pure refactor / new-keyword-only: nothing may change." Commit the manifest
   *before* (or separately from) the feature change, so intent-first is evidenced — not
   self-reported.
2. **Give the feature a regression case.** Add an `Examples/` case that exercises the new
   input/feature with `RandomSeed 0`, and commit its vetted `output.txt` as the reference —
   exactly what `Examples/CO2_NaX_Zeolite/` is for Bug #4. (The vetted reference needs a GPU run;
   until you have one, mark the case "reference pending" in its README — see the caveat in
   step 6.)
3. **Baseline runs:** build the base branch in a worktree; run every non-MLIP example →
   `runs_base/<case>/output.txt` (always `> output.txt 2>&1`).
4. **Feature runs:** same inputs, feature build → `runs_feat/<case>/output.txt`.
5. **The gate:**
   ```bash
   bash debugging/sweep_compare.sh runs_base runs_feat expected_diffs.txt
   # exit 0 = surgical. Verdicts: CONSTANT / EXPECTED-DIFF (ok) ·
   #          UNEXPECTED-DIFF / EXPECTED-MISSING / ERROR (fail)
   ```
   **UNEXPECTED-DIFF** = a regression — the failing check names what drifted (moves? counts? one
   energy component?); localize with the specificity test / standalone repro above.
   **EXPECTED-MISSING** = the feature didn't actually do what it claims — equally a finding.
6. **GPU-free shortcut for input/parser features:** extract the touched reader into a standalone
   harness (pattern: `test_case/challenge/parser_under_test.cpp`) that prints, per input file,
   the fired branch(es) and parsed values; then gate base-vs-feature traces over EVERY
   `Examples/*/simulation.input` (or `force_field.def`) with **`sweep_compare.sh --diff`**
   (byte-exact mode — plain `score.py` mode can only parse full gRASPA outputs). Seconds per
   sweep, before any GPU build. Caveats, all field-tested:
   - ⚠️ **Substring keyword matching.** gRASPA matches input keywords with bare `str.find`, so a
     new keyword that contains — or is contained by — an existing one (`Pressure`,
     `FugacityCoefficient`, …) silently triggers the wrong branch. **Standard fix: match the
     keyword as an exact first token** (tokenize the line, compare `tokens[0]`), never bare
     `find`. (A naive `Fugacity` keyword corrupted 38/44 example inputs via their
     `FugacityCoefficient` lines; the token-exact version swept clean — see `FIELD_TEST.md`.)
   - ⚠️ **Old binaries silently swallow new keywords** — and `Check_Inputs_In_read_data_cpp`
     "validates" keywords by substring-grepping the *source*, so it can falsely accept a new
     keyword on an old build too (which then quietly runs with the old behavior). State the
     minimum gRASPA version wherever the new keyword is documented.
   - The harness is a hand-copied mirror of production code: re-diff it against
     `src_clean/read_data.cpp` whenever either changes. And a clean parse-level sweep proves the
     **reader** is surgical — not the physics. The full-run gate (steps 3–5) is still required
     before merge.

## Known-Bug Catalog — commit `186e4d3` (4 bugs, all fixed; for pattern-matching)
This fork already has the fix merged (branch `fix/lj1264-forcefield-parser`, PR #81). A/B SHAs:
good ref `3a88f7f` · bug introduced `186e4d3` · fix `ea5a825` · regression-case add `cd680e1`.

| # | Sev | Where | Bug | Trigger |
|---|-----|-------|-----|---------|
| 4 | HIGH (default LJ path) | `read_data.cpp` `OverWrite_Mixing_Rule` 1st pass | last `# mixing rules to overwrite` marker wins → trailing `/0` ⇒ `Nmixrule=0` ⇒ ALL overrides dropped (→ Lorentz-Berthelot) | NaX/zeolite `force_field.def` ending in duplicate marker + `0`; MOFs (single marker) unaffected ⇒ "zeolites differ, MOFs don't" |
| 1 | UseLJ1264 | `GetTailCorrectionValue_Coeff()` | `4π` should be `2π` ⇒ TailE 2× too large | `UseLJ1264 yes` + tail on |
| 2 | UseLJ1264 | defint fixed-offset window | `GENERIC2_HC` ignored when comments precede data | TIP4PEW example (3 comment lines) |
| 3 | LOW | mixrule `size()==5` guard | 6-token `lennard-jones-1264` pair override dropped | `I J lennard-jones-1264 e s C4` |

Empirical proof for Bug #4: NaX FAU CO2 → buggy HEAD adsorbs **0 CO2** (E≈−2.8M); good-ref and
HEAD+fix both **~17 CO2, energy identical to every digit**.

## AutoJIT-gRASPA as a Cross-Reference Oracle
`AutoJIT-gRASPA` (Zhaoli2042 — the upstream gRASPA author) is the JIT/optimization sibling. Use it,
**don't debug it**: source of the vendored `score.py`; its `scripts/benchmark.sh` (compile+run+score,
`2>&1`-safe); and `examples/*/baseline_output.txt` as vetted behavioral references. Its `SKILL.md`
documents the same determinism contract ("same RNG = same moves") and an agent-agnostic loop
(Claude Code *or* Codex). ⚠️ AutoJIT's compiler block points at its own host
(`/opt/nvidia/hpc_sdk/.../24.5`) — on Quest use `Cluster-Setup/Quest_April_2026/`.

## Tips
1. Pin `RandomSeed 0` first — no determinism, no comparison.
2. Always `2>&1` — the summary is on stderr.
3. Repro parser bugs standalone before the GPU build.
4. Read `score.py`'s `failures` array; distinguish comparison-class from absolute-class.
5. A "DIFFER" with identical setup+initial energy is MC non-determinism, not a bug.
6. Bisect with worktrees, not serial checkouts.
7. When you fix: WARN, don't drop; then re-run the Examples suite through `score.py`.
