---
name: graspa-1264-tail-and-parse-bugs
description: "Two bugs in gRASPA's 12-6-4 LJ commit (186e4d3) and their fixes"
metadata: 
  node_type: memory
  type: project
  originSessionId: 696c45e2-6284-4fb0-a4fb-4e9f8e33f9be
---

gRASPA commit **186e4d3** (Mar 11 2026, "Add 12-6-4 LJ potential / UseLJ1264") — the
commit right after the Jan 4 2026 print-flag merge `3a88f7f` — introduced FOUR real bugs in
the refactored `OverWrite_Mixing_Rule()` / tail code.

**Bug 4 is THE reported "zeolites differ, MOFs don't" bug — it is in the DEFAULT LJ path
(no UseLJ1264).** (Earlier mis-diagnosis: I first concluded the regression required
UseLJ1264+tail, because the ~10 example sweep showed the default path bit-identical — but those
examples lack the trailing-marker force_field.def structure. The user's real input
(graspa_test_zeolite.zip, NaX FAU zeolite + Na cations + CO2) exposed Bug 4.)

4. **All pairwise `force_field.def` overrides silently dropped when the file has a SECOND/trailing
   `# mixing rules to overwrite` section.** The 12-6-4 refactor removed the `break` in the first
   pass and re-assigns `mixrule_startline`/`Nmixrule` on EVERY matching line (read_data.cpp:966),
   so the LAST marker wins. Real NaX/zeolite force_field.def = `# mixing rules to overwrite\n11\n
   <11 lines>\n# mixing rules to overwrite\n0` → Nmixrule read as 0 → `if(Ndefint==0&&Nmixrule==0)
   return` → all 11 overrides (incl. O-Na, O_co2-Na, C_co2-Na cation pairs) dropped → very wrong
   energies. OLD code broke after the FIRST marker → Noverwrite=11 → applied. Proven GPU-free with a
   standalone parser repro (NEW: Nmixrule=0; OLD/FIX: Nmixrule=11). MOFs unaffected (single marker,
   0 overrides anyway). Fix: first-occurrence wins — add `&& mixrule_startline == 0` (and `&&
   defint_startline == 0`) guards so a trailing empty section can't clobber the real one.
   Severity: HIGH (default path, silent wrong answer).

Bugs 1-3 are real but only trigger under `UseLJ1264 yes` and did NOT affect the user's reported case
(which uses default LJ, no tail, 5-token lennard-jones overrides):

1. **Tail correction 2× too large** — `GetTailCorrectionValue_Coeff()` in `src_clean/read_data.cpp`
   used `4.0*pi` instead of `2.0*pi`. Fix: change to `2.0*pi`. Verified: 12-6-4 tail then
   matches the standard LJ tail exactly. Only bites when tail corrections are ON (zeolite FFs
   often enable them; the shipped MOF 12-6-4 example uses `tailcorrections no`, hence
   "zeolites differ, MOFs don't").
2. **GENERIC2_HC interaction silently ignored** — the `force_field.def` parser assumed the data
   line sat at a fixed offset `[defint_startline+3, +3+Ndefint)`; comment lines before the data
   (as in the shipped `TIP4PEW-MgMOF-LJ1264` example) pushed it out of the window, so the 12-6-4
   interaction was never applied (proved: shipped example == deleting the line entirely). Fix:
   scan all lines after the count line, skip comments/blanks, consume exactly `Ndefint` real
   `GENERIC2_HC` lines (counter `Ndefint_done`). The C12=-p5/C6=p3/C4=-p2 mapping itself is correct.
3. **6-token `lennard-jones-1264` override in `force_field.def` silently ignored** (confirmed Jun 2026
   by source trace + 3 independent agents; not yet patched). `OverWrite_Mixing_Rule()` "mixing rules
   to overwrite" block guards on `if(termsScannedLined.size() == 5)` (read_data.cpp:1036). The
   `lennard-jones-1264` pair-override format documented in gRASPA/CLAUDE.md:94 has 6 tokens (extra C4)
   → fails `==5`, no else, GENERIC2_HC block needs `>=9 && [2]=="GENERIC2_HC"` → dropped, pair keeps
   default LB mixing. Corollary: no working way to set per-pair C4 in force_field.def at all. The
   per-atom `lennard-jones-1264` in force_field_mixing_rules.def IS handled (line 820). Fix: accept
   size==6 with token[2]=="lennard-jones-1264", read C4=stod(token[5])/1.20272430057 into Mix_Z — or
   fix the CLAUDE.md doc. Lower severity than #1/#2.

Consequence: the commit message's "verified to machine precision" and the shipped example's
committed `output.txt` are stale/misleading — the feature never actually ran end-to-end.
Patch + worktrees under `/home/xiaoyi/graspa_bughunt/` (FIX_1264_tailcorrection.patch).
Not yet applied to user's repo or PR'd as of 2026-06-08. See [[graspa-build-and-test]].

**Artifacts (portable handoff):** `/home/xiaoyi/graspa_bughunt/BUGFIX_HANDOFF.md` (full writeup,
ready to copy to another machine), `FIX_graspa_1264_and_override_bugs.patch` (consolidated #1/#2/#4 vs
HEAD 3fc256d; #3 documented-not-applied), `wt_fix/` (built binary), `usercase/graspa_test_zeolite/`
(NaX repro inputs), `usercase/parse_repro2.cpp` (GPU-free parser old/new/fix check). Open items:
full GCMC empirical confirm of Bug #4 (blocked by GPU CUDA-init issue), decide Bug #3, regenerate
stale TIP4PEW output.txt, open PR.

**Why it slipped through (dev-time records, checked Jun 2026):** No surviving transcript of the
Mar 11 dev session (Claude Code history has only a May 14 Ewald-printout session + Cursor snapshots
that all predate the 12-6-4 code). The surviving memory = commit msg + the example's `readme.md`.
Neither flags the bug; both claim correctness. The readme reveals the blind spots: (a) the example
is explicitly "truncated, no tail corrections" (stated 3x) → Bug #1's code path was never run;
(b) validation was "formula verified to machine precision in isolation" (correct U(r) math, but not
the parsing) + "loading vs RASPA2 within error bars" (473.9 vs 487.4 — insensitive, since one Mg-Ow
pair barely moves total loading even with the term inactive); the host-guest energy diff that would
have exposed it was attributed to "different Ewald implementations."
