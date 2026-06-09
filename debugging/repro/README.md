# GPU-free parser repro

Standalone reproductions of the `force_field.def` parser bug (Bug #4) — no GPU, no gRASPA build.
Each file inlines the exact control flow of `src_clean/read_data.cpp`'s `OverWrite_Mixing_Rule()`
first pass so you can test parsing in seconds.

```bash
g++ -O2 -std=c++17 parse_repro2.cpp -o parse_repro2
./parse_repro2 ../../Examples/CO2_NaX_Zeolite/force_field.def
```

**Expected:** buggy (last-marker-wins, `186e4d3..ea5a825`) → `Nmixrule = 0` (all 11 overrides
dropped → silent Lorentz-Berthelot fallback); fixed (first-occurrence-wins) → `Nmixrule = 11`.

`parse_repro.cpp` is the longer/annotated variant; `parse_repro2.cpp` is the minimal one used in
`../REPRODUCE.md`. Both take the `force_field.def` path as `argv[1]`.

Origin: distilled during the 2026-06-08 bug hunt. The trigger lives in
`../../Examples/CO2_NaX_Zeolite/force_field.def` (overrides block + a trailing
`# mixing rules to overwrite` / `0`, lines 17-18).
