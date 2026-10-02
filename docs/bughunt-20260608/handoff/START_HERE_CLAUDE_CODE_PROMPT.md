# Kickoff prompt for Claude Code (paste the block below)

Unzip this bundle on the new machine, `cd` into it, start Claude Code, and paste the prompt
between the lines. Adjust the two PATH placeholders if needed.

----------------------------------------------------------------------------------------------------
I'm the developer/maintainer of gRASPA (https://github.com/snurr-group/gRASPA), a GPU GCMC code for
adsorption in zeolites/MOFs. In a prior session we diagnosed bugs introduced by commit `186e4d3`
("Add 12-6-4 LJ potential / UseLJ1264", the commit right after the Jan-4 print-flag merge `3a88f7f`).
A complete writeup is in this folder — **read `BUGFIX_HANDOFF.md` first**, then `claude_memory_notes/`.

Verified findings (4 distinct bugs, all in `src_clean/read_data.cpp`, all from commit `186e4d3`):
- **Bug #4 (THE reported "zeolites differ, MOFs don't", DEFAULT LJ path):** the refactor removed the
  `break` in `OverWrite_Mixing_Rule` first pass, so the LAST `# mixing rules to overwrite` marker wins.
  A `force_field.def` with a trailing `# mixing rules to overwrite / 0` (real NaX/zeolite inputs have
  this) → `Nmixrule=0` → early return → ALL pairwise overrides silently dropped. Fix = first-occurrence
  wins (`&& mixrule_startline == 0`). Proven GPU-free with `reproduction/parse_repro2.cpp`.
- **Bug #1:** 12-6-4 tail correction prefactor `4.0*pi` should be `2.0*pi` (TailE was 2× too large;
  only with `UseLJ1264 yes` + tail corrections on).
- **Bug #2:** `GENERIC2_HC` interaction silently ignored when comment lines precede the data (fixed-offset
  parse window); only with `UseLJ1264 yes`.
- **Bug #3 (not yet applied):** 6-token `lennard-jones-1264` pair override in `force_field.def` ignored
  (`size()==5` guard). Decide: support it (suggested code in `BUGFIX_HANDOFF.md`) or fix `CLAUDE.md:94`.

Fixes for #1/#2/#4 are in `patch/FIX_graspa_1264_and_override_bugs.patch` (base commit `3fc256d`) and as a
full file in `fixed_source/read_data.cpp`. The reported test case is in `test_case_NaX/`.

Please do this, working in a NEW folder (don't modify anything outside it), and verify empirically
(build + run + bit-level diff), not just by reading code:
1. Clone gRASPA fresh (or point me at an existing clone at/near `3fc256d`). Confirm the 4 buggy lines are
   present. Apply `patch/FIX_graspa_1264_and_override_bugs.patch`.
2. Build with nvc++ (NVIDIA HPC SDK): `cd src_clean && bash NVC_COMPILE` → `nvc_main.x`.
   GPU build/run env: <FILL IN: e.g. /opt/nvidia/hpc_sdk/.../compilers/bin on PATH, CUDA_VISIBLE_DEVICES=0>.
3. **Empirically confirm Bug #4** on `test_case_NaX/`: run the NaX case at (a) pre-12-6-4 commit `3a88f7f`,
   (b) HEAD `3fc256d`, (c) HEAD+patch — same `RandomSeed 0`. Show that HEAD silently drops the 11
   `force_field.def` overrides (wrong energy/loading) while pre-12-6-4 and HEAD+patch apply them and agree.
   Use a short cycle count for the A/B; compare the `END OF SIMULATION` energy line (HGVDW/HGEwaldE/…).
4. Also run a couple of MOF examples (e.g. `Examples/CO2_MgMOF74_UFF`) to confirm the patch is a no-op there
   (no regression).
5. Decide and implement Bug #3 (code fix or doc fix).
6. Regenerate the stale `Examples/TIP4PEW-MgMOF-LJ1264/output.txt` (it was produced with 12-6-4 silently
   inactive due to Bug #2).
7. Prepare a PR to snurr-group/main with the consolidated fixes + a short description per bug.

GPU note: if a gRASPA run dies immediately with `CUDA Error: Error allocating Malloc: initialization error`
while `nvidia-smi` works, it's a transient CUDA/UVM init state (not a version mismatch) — reload the module
(`sudo modprobe -r nvidia_uvm && sudo modprobe nvidia_uvm`) and retry.

Start by reading `BUGFIX_HANDOFF.md`, confirming the buggy lines on this machine's gRASPA, then doing step 3
(the empirical confirmation of Bug #4) first.
----------------------------------------------------------------------------------------------------
