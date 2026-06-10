# QUICKSTART — run this debugging test on another agent

Goal: hand the challenge to an agent (Codex CLI or Claude Code), let it fix
`challenge/parser_under_test.cpp`, then grade it. Tier 1 is GPU-free — just needs `g++`.

```bash
REPO=$(git rev-parse --show-toplevel)   # run from inside this checkout (or set the absolute path)
TC="$REPO/debugging/test_case"
```

## 1. Make a sandbox with ONLY the agent-visible parts
Keeps the pristine challenge clean and keeps the answer key (`solution/`) out of the agent's reach.

```bash
rm -rf /tmp/graspa_challenge && mkdir -p /tmp/graspa_challenge
cp "$TC/TASK.md" "$TC/METHOD.md" /tmp/graspa_challenge/
cp -r "$TC/input" "$TC/challenge" /tmp/graspa_challenge/
cd /tmp/graspa_challenge
# NOTE: METHOD.md is the spoiler-free method. Do NOT copy ../DEBUGGING.md or ../SKILL.md into the
# sandbox — their known-bug catalog gives away the answer.
```

## 2a. Run it with **Codex CLI**
Non-interactive (scriptable):
```bash
codex exec --cd /tmp/graspa_challenge "$(cat <<'EOF'
Read TASK.md and challenge/README.md, then debug and fix the issue.
Tier 1 (required, GPU-free): edit challenge/parser_under_test.cpp so the NaX case
(input/force_field.def) applies its pairwise overrides correctly. Build/test with:
  g++ -O2 -std=c++17 challenge/parser_under_test.cpp -o /tmp/put && /tmp/put input/force_field.def
Report the root cause (file:line), the fix as a diff, and how you verified it.
Rules: don't edit the input files; keep the fix general (no special-casing).
EOF
)"
```
Interactive (works across Codex versions): `cd /tmp/graspa_challenge && codex`, then paste the same
prompt. (Flags vary by Codex version — if `--cd` isn't recognized, just `cd` first and drop it.)

## 2b. ...or with **Claude Code**
```bash
cd /tmp/graspa_challenge
claude -p "Read TASK.md, METHOD.md, and challenge/README.md and fix the bug per Tier 1. Edit
challenge/parser_under_test.cpp, keep the fix general, and report root cause (file:line), the diff,
and your verification."
```
Or interactive: `cd /tmp/graspa_challenge && claude`, then paste the prompt.
⚠️ For a fair test, tell the agent **not** to use the globally-installed `graspa-debug` skill or read
files outside the sandbox — its known-bug catalog spoils this exact answer. (Use `METHOD.md` instead.)

## 3. Grade the result
The grader compiles the candidate and checks it against its **own** copies of the inputs (so editing
inputs can't cheat). Exit 0 = PASS.
```bash
bash "$TC/solution/grade.sh" /tmp/graspa_challenge/challenge/parser_under_test.cpp
```
Expected: **FAIL** before the agent fixes it, **PASS** after. A passing run shows: `trigger
Nmixrule=11`, `CO2-Na present`, both controls + a held-out adversarial control unchanged, **and** a
whole-`Examples/` no-regression sweep (only `CO2_NaX_Zeolite` changes; the other 33 force fields are
byte-identical). Run that sweep alone with:
```bash
bash "$TC/solution/regression_examples.sh" /tmp/graspa_challenge/challenge/parser_under_test.cpp
```
Then skim the agent's written explanation against the rubric in
[`solution/SOLUTION.md`](solution/SOLUTION.md).

## 4. Reset between runs
```bash
rm -rf /tmp/graspa_challenge        # then re-do step 1
```

## Sanity check (no agent needed)
```bash
bash "$TC/solution/grade.sh" "$TC/challenge/parser_under_test.cpp"   # FAIL (unfixed bug)
bash "$TC/solution/grade.sh" "$TC/solution/parser_fixed.cpp"         # PASS (reference fix)
```

## Optional Tier 2 (full GPU)
Build gRASPA at `HEAD` vs `186e4d3`, run `input/` (it pins `RandomSeed 0`), and compare with the
constant-results gate: `python3 "$REPO/debugging/score.py" <output> "$TC/solution/reference_output.txt"`.
See [`../REPRODUCE.md`](../REPRODUCE.md). Buggy ≈ 0 CO₂; fixed ≈ 17–20 CO₂.
