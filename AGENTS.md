# AGENTS.md

Entry point for coding agents (Codex CLI, Claude Code, etc.) working in this gRASPA checkout.

## Debugging gRASPA correctness? Start here
New here? Read **[`debugging/TUTORIAL.md`](debugging/TUTORIAL.md)** first (quick step-by-step).
Then **[`debugging/DEBUGGING.md`](debugging/DEBUGGING.md)** (the playbook) and
**[`debugging/REPRODUCE.md`](debugging/REPRODUCE.md)** (a copy-paste reproducible run). The matching
Claude Code skill is **[`debugging/SKILL.md`](debugging/SKILL.md)** (install it by copying to
`.claude/skills/graspa-debug/SKILL.md`; `.claude/` is gitignored here). A ready-to-run challenge for
testing another agent is in **[`debugging/test_case/`](debugging/test_case/)**.

## The one rule that matters: results must be CONSTANT
gRASPA is deterministic for a fixed `RandomSeed`. To prove a bug or a fix, compare two runs of the
**same input with `RandomSeed 0`** and check they are bit-comparable, using the vendored gate:

```bash
python3 debugging/score.py <test_output.txt> <reference_output.txt>
```

Read the JSON `failures` array, not just the exit code: **comparison-class** failures (moves /
counts / energies vs reference) mean the runs differ; **absolute-class** failures (energy drift,
structure factor) are a property of the single run and fire even on a self-compare.

GOTCHA: gRASPA prints its final energy + loadings to **stderr** — always run with
`> output.txt 2>&1`, or `score.py` sees a truncated file and falsely reports "identical."

## Build
`src_clean/NVC_COMPILE`. Quest recipe in `Cluster-Setup/Quest_April_2026/`. On a Quest LOGIN node add
`-tp haswell` or the binary SIGILLs (exit 132) on A100 nodes; energy kernels are CUDA so `-tp` never
affects results.

## When you fix parser/override code
**WARN on unrecognized lines; never silently drop them**, then re-run the whole `Examples/` suite
through `debugging/score.py` to prove the fix is surgical. See the `186e4d3` bug catalog in
`debugging/DEBUGGING.md`.
