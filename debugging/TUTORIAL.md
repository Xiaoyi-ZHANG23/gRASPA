# Quick Tutorial — gRASPA debugging kit

This branch adds a small, self-contained **debugging kit** for gRASPA under `debugging/`:

- a **skill / playbook** for hunting gRASPA correctness bugs deterministically, and
- a ready-to-run **debugging test case** you can hand to another AI agent (Codex / Claude Code) to
  practice on, with an automated grader.

Everything is GPU-free at its core (a `g++` one-liner); a full-GPU tier is optional.

---

## What's here
```
debugging/
├── TUTORIAL.md        <- you are here
├── README.md          index of the kit
├── SKILL.md           the graspa-debug skill (Claude Code) — install: copy to .claude/skills/graspa-debug/
├── DEBUGGING.md       the same playbook, tool-agnostic (Codex/human)
├── REPRODUCE.md       copy-paste recipe to reproduce the real bug (parser repro + GPU A/B)
├── score.py           constant-results comparator (vendored from AutoJIT-gRASPA)
├── selftest.sh        one-command health check of the GPU-free kit
├── repro/             standalone g++ reproductions of the parser bug
└── test_case/         the agent-testable debugging challenge (see below)
```
Repo root `AGENTS.md` points Codex/agents here.

---

## Part A — Run the debugging test case on another agent (5 minutes)

The challenge reproduces a real, already-fixed gRASPA bug: pairwise `force_field.def` overrides being
silently dropped, so a CO₂/NaX-zeolite case adsorbs **0** CO₂ instead of ~17–20 (MOFs unaffected).

### 1. Verify it works (no agent)
```bash
REPO=$(git rev-parse --show-toplevel)
TC=$REPO/debugging/test_case

bash "$REPO/debugging/selftest.sh"   # health-checks the whole GPU-free kit, or just the grader:
bash "$TC/solution/grade.sh" "$TC/challenge/parser_under_test.cpp"   # RESULT: FAIL (the unfixed bug)
bash "$TC/solution/grade.sh" "$TC/solution/parser_fixed.cpp"         # RESULT: PASS (reference fix)
```

### 2. Give the agent a sandbox (only the agent-visible parts; withholds the answer key)
```bash
rm -rf /tmp/graspa_challenge && mkdir -p /tmp/graspa_challenge
cp "$TC/TASK.md" "$TC/METHOD.md" /tmp/graspa_challenge/
cp -r "$TC/input" "$TC/challenge" /tmp/graspa_challenge/
cd /tmp/graspa_challenge
```

### 3. Launch the agent

**Codex CLI**
```bash
codex exec --cd /tmp/graspa_challenge "$(cat <<'EOF'
Read TASK.md and challenge/README.md, then debug and fix the issue (Tier 1, GPU-free).
Build/test with: g++ -O2 -std=c++17 challenge/parser_under_test.cpp -o /tmp/put && /tmp/put input/force_field.def
Edit challenge/parser_under_test.cpp so the NaX overrides are applied. Report root cause (file:line),
the fix as a diff, and your verification. Don't edit input files; keep the fix general.
EOF
)"
```
(If your Codex version rejects `--cd`, `cd` into the folder first and drop the flag, or run `codex`
interactively and paste the prompt.)

**Claude Code**
```bash
cd /tmp/graspa_challenge
claude -p "Read TASK.md, METHOD.md, and challenge/README.md and fix the bug (Tier 1). Edit
challenge/parser_under_test.cpp, keep the fix general, and report root cause (file:line), the diff,
and your verification. Do NOT use the installed graspa-debug skill or read files outside this folder."
```
> For a fair test, tell the agent not to use the globally-installed `graspa-debug` skill — its
> known-bug catalog spoils the answer. `METHOD.md` is the spoiler-free method to use instead.

### 4. Grade
```bash
bash "$TC/solution/grade.sh" /tmp/graspa_challenge/challenge/parser_under_test.cpp
echo "exit=$?"   # 0 = PASS, 1 = FAIL
```
A PASS means: trigger applies all 11 overrides (incl. CO₂–Na⁺), the controls + a held-out adversarial
control are unchanged, and a sweep over **every** `Examples/*/force_field.def` shows the fix changed
**only** `CO2_NaX_Zeolite` (the other 33 are byte-identical). Then compare the agent's written
explanation to the rubric in `test_case/solution/SOLUTION.md`.

### 5. Reset for the next agent
```bash
rm -rf /tmp/graspa_challenge   # repeat from step 2
```

---

## Part B — Use the playbook on a real gRASPA bug

Read `debugging/DEBUGGING.md` (or install `debugging/SKILL.md` as a Claude Code skill). The method in
one line: gRASPA is **deterministic for a fixed `RandomSeed`**, so compare two runs that should be
identical and find where they stop being identical. Tools:

```bash
# constant-results comparison of two runs (capture BOTH streams — the summary is on stderr):
./nvc_main.x > a.txt 2>&1
python3 debugging/score.py a.txt reference.txt   # 0 = constant; read the JSON "failures" array

# A/B across commits, kept side-by-side:
git worktree add ../wt_good <good_sha> && git worktree add ../wt_head HEAD
```
Full runnable example (the bug this kit is built around): `debugging/REPRODUCE.md`.

---

## Optional: full-GPU tier
Build gRASPA (`src_clean/NVC_COMPILE`; on Quest see `Cluster-Setup/Quest_April_2026/`), run
`test_case/input/` (already `RandomSeed 0`), and compare against the shipped reference with
`debugging/score.py`. Buggy build ≈ 0 CO₂; fixed ≈ 17–20 CO₂. Details in `debugging/REPRODUCE.md`.
