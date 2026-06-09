# test_case/ — a self-contained gRASPA debugging challenge

A ready-to-run **debugging exercise** for testing another agent (Claude Code, Codex, …) on the
`graspa-debug` workflow. It reproduces a real, already-fixed gRASPA bug (the headline "Bug #4" — see
`../DEBUGGING.md`): pairwise `force_field.def` overrides silently dropped for inputs whose force
field ends in a trailing `# mixing rules to overwrite` / `0` block, so a CO₂/NaX-zeolite case adsorbs
**0** instead of **~17–20** CO₂, while MOF cases are unaffected.

The challenge is **deterministic** and **GPU-free** at its core (Tier 1 is a `g++` one-liner), with
an optional full-GPU tier.

> **In a hurry?** [`QUICKSTART.md`](QUICKSTART.md) has copy-paste commands to launch this on Codex
> CLI or Claude Code and grade the result.

## How to run a test on another agent
1. **Give the agent** (and nothing else from here):
   - [`TASK.md`](TASK.md) — the prompt / symptom report (no spoilers)
   - [`input/`](input/) — the full gRASPA input set (the "sample debugging input"), `RandomSeed 0`
   - [`challenge/`](challenge/) — the GPU-free reproduction harness it will fix
   - it may also use the playbook `../DEBUGGING.md` / skill `../SKILL.md`
2. **Withhold** [`solution/`](solution/) — that's the answer key + automated grader.
3. **Grade** the agent's edited `challenge/parser_under_test.cpp`:
   ```bash
   bash solution/grade.sh                 # grades ../challenge/parser_under_test.cpp
   # or: bash solution/grade.sh <path-to-candidate.cpp>
   ```
   Exit 0 = PASS. The grader compiles the candidate and checks observable behavior on its **own**
   copies of the trigger + three controls (including a **held-out adversarial** input the agent never
   sees, which fails non-general "cheat" fixes), so the agent can't pass by editing inputs or by
   special-casing.
4. Read [`solution/SOLUTION.md`](solution/SOLUTION.md) for the root cause, the fix, the expected
   results for both tiers, and a grading rubric (the grader is necessary, not sufficient — the rubric
   covers fix *quality*).

## Layout
```
test_case/
├── TASK.md                      # AGENT-VISIBLE prompt (symptom only)
├── input/                       # AGENT-VISIBLE gRASPA inputs (RandomSeed 0); the trigger force_field.def
├── challenge/                   # AGENT-VISIBLE GPU-free repro the agent fixes
│   ├── parser_under_test.cpp    #   contains the bug; agent edits this
│   ├── ff_control_mof.def       #   single-marker control (with overrides)
│   ├── ff_control_zero.def      #   single-marker control (no overrides)
│   └── README.md
└── solution/                    # HUMAN-ONLY answer key — withhold from the agent
    ├── SOLUTION.md              #   root cause, fix, expected results, rubric
    ├── parser_fixed.cpp         #   reference fix
    ├── parser_buggy_baseline.cpp#   pristine unfixed baseline (for the no-regression sweep)
    ├── grade.sh                 #   automated grader (5 checks + all-Examples no-regression)
    ├── regression_examples.sh   #   no-regression sweep over every Examples/*/force_field.def
    ├── grade_*.def              #   grader's own copies of trigger + controls (incl. held-out adversarial)
    └── reference_output.txt     #   correct full-GPU output (for Tier-2 score.py grading)
```

## Quick self-check (it works out of the box)
```bash
bash solution/grade.sh challenge/parser_under_test.cpp   # FAIL (the unfixed bug)
bash solution/grade.sh solution/parser_fixed.cpp         # PASS (the reference fix)
```
