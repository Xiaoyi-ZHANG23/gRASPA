# gRASPA 12-6-4 / force_field.def Bug Investigation — Handoff

**Date:** 2026-06-08  ·  **Repo:** snurr-group/gRASPA  ·  **Latest verified:** HEAD `3fc256d` (== origin/main, up to date)
**Workspace:** `/home/xiaoyi/graspa_bughunt/` (self-contained; nothing outside it was modified)

## TL;DR
Commit **`186e4d3`** ("Add 12-6-4 LJ potential / UseLJ1264", Mar 11 2026 — the commit right after the
Jan 4 print-flag merge `3a88f7f`) refactored `OverWrite_Mixing_Rule()` and the tail-correction code and
introduced **4 distinct bugs**. The **reported "zeolites differ, MOFs don't"** symptom is **Bug #4**, which is
in the **default LJ path** (no `UseLJ1264` needed). Bugs #1–#3 are real but only trigger under `UseLJ1264 yes`.

## The bugs

### Bug #4 — pairwise `force_field.def` overrides silently dropped  ★ THE REPORTED BUG ★
- **File/line:** `src_clean/read_data.cpp:966` (`OverWrite_Mixing_Rule`, first pass).
- **Root cause:** the refactor removed the old `break` and re-assigns `mixrule_startline = counter` on
  **every** `# mixing rules to overwrite` line, so the **LAST** marker wins. A `force_field.def` that ends
  with a trailing `# mixing rules to overwrite` / `0` section makes `Nmixrule = 0` →
  `if(Ndefint==0 && Nmixrule==0) return;` → **all overrides dropped** (falls back to Lorentz-Berthelot).
- **Trigger:** real NaX/zeolite input (see test case) whose `force_field.def` is
  `# mixing rules to overwrite \n 11 \n <11 lines> \n # mixing rules to overwrite \n 0`.
  Dropped pairs include the Na⁺ cation interactions (O–Na, O_co2–Na=200.8, C_co2–Na=362.3) and all CO₂ pairs.
- **Why MOFs are fine:** their `force_field.def` has a single marker with 0 overrides → unchanged.
- **Default path** — does NOT require UseLJ1264. This is why the earlier 10-example sweep missed it (none had
  the trailing-marker structure).
- **Fix (APPLIED):** first-occurrence wins — add `&& mixrule_startline == 0` (and `&& defint_startline == 0`).
- **Verified:** standalone parser repro → buggy `Nmixrule=0`, fixed `Nmixrule=11`; MOF unchanged. Patched source compiles.
- **Severity: HIGH** (default path, silent wrong answer).

### Bug #1 — 12-6-4 tail correction 2× too large
- **File/line:** `read_data.cpp:856`, `GetTailCorrectionValue_Coeff()`.
- **Cause:** prefactor `4.0 * pi` should be `2.0 * pi` (per-pair tail = `2π·∫U r² dr`).
- **Trigger:** `UseLJ1264 yes` AND tail corrections on. **Verified empirically:** TailE was exactly 2× the LJ value (−31490.73 vs −15745.36); after fix they match exactly.
- **Fix (APPLIED):** `4.0 * pi` → `2.0 * pi`.

### Bug #2 — `GENERIC2_HC` interaction silently ignored when comments precede the data
- **File/line:** `read_data.cpp:989`, defint section fixed-offset window `[defint_startline+3, +3+Ndefint)`.
- **Cause:** comment lines between the count and the data push the real line out of the window.
- **Trigger:** `UseLJ1264 yes` + a `force_field.def` GENERIC2_HC block with comments (the shipped TIP4PEW example has 3). **Verified:** shipped example == deleting the GENERIC2_HC line entirely (12-6-4 never applied).
- **Fix (APPLIED):** scan after the count line, skip comments/blanks, consume exactly `Ndefint` real `GENERIC2_HC` lines (counter `Ndefint_done`).

### Bug #3 — 6-token `lennard-jones-1264` pair override in `force_field.def` ignored (NOT YET APPLIED)
- **File/line:** `read_data.cpp:1036`, mixrule block guard `if(termsScannedLined.size() == 5)`.
- **Cause:** the documented (CLAUDE.md:94) 6-token `Mg O_co2 lennard-jones-1264 e s C4` override fails `==5`;
  the GENERIC2_HC block needs `>=9 && [2]=="GENERIC2_HC"`. No branch matches → silently dropped. (The per-atom
  `lennard-jones-1264` form in `force_field_mixing_rules.def` IS handled, line 820.)
- **Verified:** source trace + 3 independent agents, unanimous "ignored".
- **Severity:** LOW (less-common route, doc-driven). **Decision needed:** support the format OR fix the doc.
- **Suggested code fix (apply on the other machine if you want the format supported):** in the mixrule block,
  accept `size()==5 || size()==6`; when `size()==6 && [2]=="lennard-jones-1264"` read
  `C4 = std::stod(termsScannedLined[5]) / 1.20272430057` and, in the `Input.Use1264` branch, write it to
  `Input.Mix_Z[typeI*FFsize+typeJ]` / `[typeJ*FFsize+typeI]` (instead of 0.0) and pass it to
  `Get_Shifted_Value_Coeff(...)`. Note: if `UseLJ1264` is off the LJ kernel ignores C4 — decide whether to warn.

## Files in this workspace
- `FIX_graspa_1264_and_override_bugs.patch` — consolidated patch for **#1, #2, #4** (against HEAD `3fc256d`).
- `wt_fix/` — git worktree of HEAD with #1/#2/#4 applied (binary built: `wt_fix/src_clean/nvc_main.x`).
- `wt_head/`, `wt_1264/`, `wt_pre1264/`, `wt_preJan4/` — worktrees at HEAD / `186e4d3` / `3a88f7f` / `0580d6c` for A/B.
- `usercase/graspa_test_zeolite/` — the user's NaX FAU zeolite reproduction inputs (from `graspa_test_zeolite.zip`).
- `usercase/parse_repro2.cpp` — **GPU-free** standalone reproduction of the parser (old vs new vs fixed). Run:
  `g++ -O2 -std=c++17 parse_repro2.cpp -o parse_repro2 && ./parse_repro2 <force_field.def>`

## Apply the fix on another machine
```bash
cd <your gRASPA clone>            # should be on/near snurr-group main 3fc256d
git apply /path/to/FIX_graspa_1264_and_override_bugs.patch
# or copy the 3 hunks manually into src_clean/read_data.cpp (lines ~856, ~959-967, ~991-998)
```

## Build & A/B test recipe
- Compiler: `nvc++` from NVIDIA HPC SDK 24.5 (`/opt/nvidia/hpc_sdk/Linux_x86_64/24.5/compilers/bin`).
- Build: `cd src_clean && bash NVC_COMPILE` → `nvc_main.x` (~20 s).
- Determinism: `RandomSeed 0` → bit-comparable runs. Compare the `END OF SIMULATION` energy line
  (`HGVDW/HGEwaldE/...`) or md5 the numeric block. Reduce cycles for quick A/B.
- A/B across commits: `git worktree add ../wt_<tag> <commit>`, build each, run same input, diff.

## OPEN ITEMS (continue here)
1. **Full GCMC empirical confirmation of Bug #4** on the NaX case: run `usercase/graspa_test_zeolite/` at
   `wt_preJan4` (or `wt_pre1264`) vs `wt_head` vs `wt_fix`, same seed; expect pre-12-6-4 ≈ fixed (overrides
   applied) and HEAD wrong (overrides dropped). **Blocked here** by a CUDA init failure (see below); parsing
   proof is already conclusive.
2. **Decide Bug #3:** apply the suggested code fix (above) or correct `CLAUDE.md:94`.
3. **Regenerate** the stale `Examples/TIP4PEW-MgMOF-LJ1264/output.txt` reference (it was produced with 12-6-4
   silently inactive due to Bug #2).
4. **Open a PR** to snurr-group/main with the consolidated fixes once empirically confirmed.

## Environment note (this machine only)
GPU hit a transient CUDA init error mid-session (`cudaErrorInitializationError`) while `nvidia-smi` works and
driver/lib versions all match 575.57.08 — i.e. NOT a version mismatch. Fix needs a sudo module reload
(`~/fix_nvidia_driver.sh`, or just `sudo modprobe -r nvidia_uvm && sudo modprobe nvidia_uvm`, which is
CUDA-only and shouldn't disturb the display). Not required on a healthy machine.

## Key commits
- `0580d6c` (Dec 5) pre-Jan-4 · `3a88f7f` (Jan 4) print-flag merge · `186e4d3` (Mar 11) **12-6-4 — introduced all 4 bugs** · `3fc256d` (May 13) HEAD.
