---
name: graspa-build-and-test
description: "How to build and A/B-test gRASPA on this machine (nvc++, GPUs, examples)"
metadata: 
  node_type: memory
  type: reference
  originSessionId: 696c45e2-6284-4fb0-a4fb-4e9f8e33f9be
---

Building/running gRASPA (`snurr-group/gRASPA`) on this server:

- Compiler: `nvc++` from the NVIDIA HPC SDK at
  `/opt/nvidia/hpc_sdk/Linux_x86_64/24.5/compilers/bin` (add to PATH). System
  also has CUDA 11.5 + g++ 11.4. GPUs: RTX 3080 Ti (gpu0) + RTX 3090 (gpu1).
- Build: `cd src_clean && bash NVC_COMPILE` → produces `nvc_main.x` in ~20 s.
- Run an example: copy an `Examples/<case>/` dir, drop in `nvc_main.x`, `./nvc_main.x > run.out`.
  Each run ~3 s for a few hundred cycles. `RandomSeed 0` makes runs deterministic, so two
  builds on the same input are bit-comparable.
- A/B-testing commits: use `git worktree add` per commit, build each, run the same input,
  diff the `END OF SIMULATION` energy line (`HGVDW/HGEwaldE/...`) or md5 the numeric block.
- Zeolite example: `CO2-MFI` (uses default LJ path, charges from CIF, no `force_field.def`).
  MOF examples: `CO2_MgMOF74_UFF`, etc. 12-6-4 example: `TIP4PEW-MgMOF-LJ1264`.
- Workspace for the Jan-2026 12-6-4 bug hunt: `/home/xiaoyi/graspa_bughunt/` (worktrees +
  fix + `BUGFIX_HANDOFF.md` + `FIX_graspa_1264_and_override_bugs.patch` + `usercase/` NaX repro +
  `usercase/parse_repro2.cpp` GPU-free parser check). See [[graspa-1264-tail-and-parse-bugs]].
- GPU gotcha (seen 2026-06-08): mid-session `cudaErrorInitializationError` ("Error allocating
  Malloc: initialization error") while nvidia-smi works and driver/lib versions match (575.57.08).
  Not a version mismatch — needs a sudo module reload: `~/fix_nvidia_driver.sh`, or minimally
  `sudo modprobe -r nvidia_uvm && sudo modprobe nvidia_uvm` (CUDA-only, won't disturb the display).
  The gRASPA binary crashes at the first cudaMalloc *before* any input parsing, so no output appears.
