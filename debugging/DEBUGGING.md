# gRASPA Debugging Playbook (portable / tool-agnostic)

> Harness-agnostic mirror of `.claude/skills/graspa-debug/SKILL.md` (the authoritative copy — if the
> two diverge, defer to it). Readable by Codex, other agents, and humans. For a copy-paste
> reproducible run, see [`REPRODUCE.md`](REPRODUCE.md).

Debug **gRASPA** (NVIDIA-GPU Monte Carlo adsorption simulator) energy/correctness discrepancies and
`force_field.def` parser bugs by proving results are **constant** (deterministic, bit-comparable)
when they should be, and pinpointing the exact line where they stop being constant.

Governing principle: gRASPA is **deterministic for a fixed `RandomSeed`** — same seed ⇒ same RNG ⇒
same move sequence ⇒ bit-for-bit identical energies and loadings. "Constant results" is both the
test oracle and the definition of done.

## When to use
- Energy/loading mismatch between two gRASPA versions, commits, or branches.
- "Forcefield overrides ignored / results fell back to Lorentz-Berthelot."
- A symptom that depends on input *structure* ("zeolites differ but MOFs don't").
- Validating that a fix changes ONLY the intended case, byte-identical elsewhere.
- Any "the number looks wrong but it didn't crash" (silent-wrong-answer class).

## Build (Quest)
- `nvhpc/23.3-gcc-10.4.0`; recipe in-repo at `../Cluster-Setup/Quest_April_2026/`. Build entry:
  `../src_clean/NVC_COMPILE`.
- ⚠️ Building on the LOGIN node: add **`-tp haswell`** or the binary dies **SIGILL (exit 132)** on
  older A100 nodes (`-tp` never affects results; energy kernels are CUDA). Login-node `module load`
  no-ops; set PATH directly or build inside the SLURM job (`gengpu`, account `p32082`).

## Core method: worktree A/B + RandomSeed 0
1. **Pin `RandomSeed 0`** in every `simulation.input` you compare — that is what makes results
   constant. Without it, move sequences diverge and the comparison is meaningless.
2. **One git worktree per suspect version:** `git worktree add ../wt_<name> <sha>`. Keep
   good-ref / current-HEAD / HEAD-plus-fix side by side.
3. **Build each** (`../src_clean/NVC_COMPILE`; `-tp haswell` on login node).
4. **Run identical input**, capturing **`> output.txt 2>&1`** (stderr gotcha below).
5. **Compare with `score.py`** (next section).
6. **To localize a regression, bisect** with worktrees against the good baseline. (This isolated
   `186e4d3`, the commit right after Jan-4 merge `3a88f7f`.)

## The constant-results gate: `score.py`  ★
Don't hand-grep "Total Energy" and eyeball digits — that misses move-count and loading drift.

```bash
python3 debugging/score.py <test_output.txt> <reference_output.txt>
# exit 0 = CONSTANT (identical within tol)   1 = differs   2 = parse error
```

`score.py` is vendored here from AutoJIT-gRASPA (Zhaoli2042) so this workflow is self-contained.
It runs two CLASSES of check — read the JSON `failures` array, not just the exit code:
- **comparison (test vs reference):** move statistics exact (same RNG ⇒ identical move counts);
  PseudoAtom/component counts exact (e.g. "17 CO2"); final energies within 1e-5 rel / 1e-3 abs.
  Failures here ⇒ **the two runs genuinely differ** (your A/B signal).
- **absolute (test file only):** per-component energy drift < 3e-5; structure factor CPU/GPU < 1e-4.
  These describe the *test run itself* and fire **even on a self-compare**.

⚠️ Verified gotcha: `score.py Examples/CO2_NaX_Zeolite/output.txt <same file>` returns **exit 1**,
not 0 — that run carries an inherent `vdw_hh`/`total` drift ≈ 1.65 (absolute-class), with zero
comparison-class failures. So a nonzero exit ≠ "the runs differ." For A/B equality, the decisive
evidence is **no comparison-class entries in `failures`** (and both runs showing the *same* drift).

(One-shot compile+run+score in AutoJIT's layout: its `scripts/benchmark.sh`, already `2>&1`-safe.)

## stderr GOTCHA
gRASPA writes its end-of-run summary — `PRODUCTION AVERAGE ENERGY` and final
`PseudoAtom Type: ... #:` loadings — to **stderr**. Capture stdout only and you compare truncated
files and wrongly conclude "IDENTICAL." Always `> output.txt 2>&1`; `score.py` needs the merged
capture.

## Specificity test
Delete ONLY the suspected trigger from the input and re-run HEAD: if HEAD then matches the good
reference (no comparison-class failures), that trigger is the sole cause and the fix doesn't touch
the normal path. (Bug #4 trigger: trailing duplicate `# mixing rules to overwrite` / `0` —
`Examples/CO2_NaX_Zeolite/force_field.def` lines 17-18.)

## Standalone parser repro (no GPU build)
Parser bugs live in `../src_clean/read_data.cpp` (`OverWrite_Mixing_Rule()`). Vendored templates
reproduce Bug #4 with plain `g++` in seconds:

```bash
cd debugging/repro
g++ -O2 -std=c++17 parse_repro2.cpp -o parse_repro2
./parse_repro2 ../../Examples/CO2_NaX_Zeolite/force_field.def   # buggy → Nmixrule=0 ; fixed → 11
```

## Recurring failure mode: "never silently drop"
The `186e4d3` bug class all silently discarded force-field overrides → wrong-but-quiet. Audit for:
- **Last-marker-wins** instead of first-occurrence (fix idiom: `... && *_startline == 0`).
- **Fixed-offset data windows** that break when comments precede data (scan past count, skip
  comments, consume exactly N).
- **Over-tight token-count guards** (`size()==5` drops a valid 6-token
  `I J lennard-jones-1264 eps sig C4` line).
- **Tail-correction prefactor** `4π` where it should be `2π` (TailE 2×).

Rule when fixing: **WARN on any unrecognized/unconsumed override line; never drop it silently.**

## Regression suite
After a fix, run the whole non-MLIP `../Examples/` suite head-vs-fix and `score.py` each pair:
expect no comparison-class failures except the intentionally-changed case(s), which should now match
the good reference. `../Examples/CO2_NaX_Zeolite/` is the Bug #4 regression case (`RandomSeed 0`,
trailing-marker `force_field.def`, reference `output.txt` — drop its first exe-path banner line when
matching shipped format).

## Known-bug catalog — commit `186e4d3` (4 bugs, fixed; for pattern-matching)
Fork already has the fix merged (branch `fix/lj1264-forcefield-parser`, PR #81). A/B SHAs: good ref
`3a88f7f` · bug `186e4d3` · fix `ea5a825` · regression-case add `cd680e1`.

| # | Sev | Where | Bug | Trigger |
|---|-----|-------|-----|---------|
| 4 | HIGH (default path) | `read_data.cpp` `OverWrite_Mixing_Rule` 1st pass | last `# mixing rules to overwrite` marker wins → trailing `/0` ⇒ `Nmixrule=0` ⇒ all overrides dropped (→ Lorentz-Berthelot) | NaX/zeolite `force_field.def` ending in duplicate marker + `0`; MOFs unaffected ⇒ "zeolites differ, MOFs don't" |
| 1 | UseLJ1264 | `GetTailCorrectionValue_Coeff()` | `4π` should be `2π` ⇒ TailE 2× | `UseLJ1264 yes` + tail on |
| 2 | UseLJ1264 | defint fixed-offset window | `GENERIC2_HC` ignored when comments precede data | TIP4PEW example (3 comment lines) |
| 3 | LOW | mixrule `size()==5` guard | 6-token `lennard-jones-1264` override dropped | `I J lennard-jones-1264 e s C4` |

Proof for Bug #4: NaX FAU CO2 → buggy HEAD adsorbs **0 CO2** (E≈−2.8M); good-ref and HEAD+fix both
**~17 CO2, energy identical to every digit.**

## AutoJIT-gRASPA as oracle (reuse, don't debug)
`AutoJIT-gRASPA` (Zhaoli2042 — upstream gRASPA author) is the JIT sibling. Source of the vendored
`score.py`; its `scripts/benchmark.sh` (compile+run+score, `2>&1`-safe); and
`examples/*/baseline_output.txt` as vetted behavioral references. Its `SKILL.md` documents
"same RNG = same moves" and an agent-agnostic loop (Claude Code *or* Codex). ⚠️ AutoJIT's compiler
block points at its own host (`/opt/nvidia/hpc_sdk/.../24.5`); on Quest use
`../Cluster-Setup/Quest_April_2026/`.

## Tips
1. Pin `RandomSeed 0` first — no determinism, no comparison.
2. Always `2>&1` — the summary is on stderr.
3. Repro parser bugs standalone before the GPU build.
4. Read `score.py`'s `failures` array; distinguish comparison-class from absolute-class.
5. A "DIFFER" with identical setup+initial energy is MC non-determinism, not a bug.
6. Bisect with worktrees, not serial checkouts.
7. When you fix: WARN, don't drop; then re-run the Examples suite through `score.py`.
