# METHOD — the gRASPA debugging method (spoiler-free)

A condensed, **answer-free** version of the `graspa-debug` playbook, safe to give to the agent under
test. It describes *how* to debug a gRASPA correctness discrepancy — not what this particular bug is.
(The full playbook in `../DEBUGGING.md` contains a known-bug catalog; don't read that during the
exercise — it spoils the answer.)

## Core principle: results must be CONSTANT
gRASPA is **deterministic for a fixed `RandomSeed`**: same seed ⇒ same RNG ⇒ same move sequence ⇒
bit-for-bit identical energies and loadings. So any correctness investigation reduces to: *make two
runs that should be identical, and find where/why they stop being identical.* The provided input pins
`RandomSeed 0` for exactly this reason.

## Techniques (apply the ones that fit)
1. **Read the relevant parser.** The override-parsing path is `OverWrite_Mixing_Rule()` in
   `src_clean/read_data.cpp`. Trace what it does with the input you were given.
2. **Standalone repro (GPU-free).** Don't rebuild the whole GPU code to test parsing — a small
   standalone program that mimics the parse lets you iterate in seconds. (One is provided in
   `challenge/`.) Run it on the failing input and compare what it *does* to what the input *says*.
3. **Specificity test.** Change exactly one thing in the input and see whether the symptom flips.
   This isolates the precise trigger and proves your fix targets it (and only it).
4. **Compare against a known-good version (A/B).** Use `git worktree add <dir> <commit>` to build an
   older "good" commit alongside the current one, run the *same* input on both with `RandomSeed 0`,
   and diff the results. `git log`/`git bisect` narrows which commit changed the behavior.
5. **Compare rigorously, not by eye.** A constant-results comparator (`debugging/score.py`) checks
   move counts, molecule counts, and per-component energies — far better than grepping one number.
   (Capture gRASPA output with `> output.txt 2>&1`; its end-of-run summary goes to **stderr**.)

## Principle for the fix
Prefer a **general** fix over a special case. A correct parser should never *silently drop* valid
input — apply it, or warn. Always show that inputs which already worked are **unchanged**
(no regression).
