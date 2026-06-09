# debugging/ — gRASPA deterministic bug-hunting kit

A self-contained, reproducible toolkit for debugging gRASPA **correctness** (energy/loading
discrepancies, `force_field.def` parser bugs). Built around one invariant: gRASPA is **deterministic
for a fixed `RandomSeed`**, so a correct binary produces **constant (bit-comparable)** results — and
any deviation localizes a bug.

## Contents
| File | What |
|------|------|
| [`TUTORIAL.md`](TUTORIAL.md) | **Start here** — quick step-by-step: run the test case on an agent, or use the playbook on a real bug. |
| [`SKILL.md`](SKILL.md) | **The Claude Code skill** (with frontmatter). Tracked/uploadable copy. |
| [`DEBUGGING.md`](DEBUGGING.md) | The full playbook (tool-agnostic; Codex/human-readable; mirrors SKILL.md). |
| [`REPRODUCE.md`](REPRODUCE.md) | Copy-paste recipe: GPU-free parser repro + full GPU A/B. |
| [`score.py`](score.py) | **Constant-results gate.** Vendored from AutoJIT-gRASPA. `python3 score.py <test> <ref>`. |
| [`repro/`](repro/) | Standalone `g++` reproductions of the parser bug (no GPU build). |
| [`test_case/`](test_case/) | **Ready-to-run debugging challenge** for testing another agent (Codex/Claude) — symptom prompt, GPU-free repro, answer key + automated grader. See its `QUICKSTART.md`. |

**Activate the skill for Claude Code:** copy `SKILL.md` to `.claude/skills/graspa-debug/SKILL.md`.
This repo's `.gitignore` excludes `.claude/`, which is why the tracked/uploadable copy lives here in
`debugging/` (same pattern AutoJIT-gRASPA uses for its root `SKILL.md`). For **Codex**, see
`../AGENTS.md`.

## 60-second start
```bash
# GPU-free proof of the headline bug:
cd repro && g++ -O2 -std=c++17 parse_repro2.cpp -o parse_repro2
./parse_repro2 ../../Examples/CO2_NaX_Zeolite/force_field.def   # buggy → Nmixrule=0 ; fixed → 11
```
Full GPU A/B (worktrees + `score.py`): see [`REPRODUCE.md`](REPRODUCE.md).

## Two things that will bite you
1. **stderr** — gRASPA prints its end-of-run energy + loadings to **stderr**. Always run with
   `> output.txt 2>&1`, or `score.py` sees a truncated file and falsely reports "identical."
2. **`score.py` exit code ≠ "runs differ."** It mixes *comparison-vs-reference* checks (moves,
   counts, energies — these mean the runs differ) with *absolute self-checks* (energy drift,
   structure factor — these fire even on a self-compare). Read the JSON `failures` array and judge
   A/B equality by the **absence of comparison-class failures**. (e.g. `CO2_NaX_Zeolite` self-compare
   exits 1 due to an inherent `vdw_hh` drift ≈ 1.65, with zero comparison failures.)

## Attribution
`score.py` is vendored unmodified from **AutoJIT-gRASPA** (https://github.com/Zhaoli2042/AutoJIT-gRASPA),
the JIT/optimization sibling project, whose correctness gate defines "results are constant" for
gRASPA. Bundled here so this debugging workflow runs from this repo alone.
