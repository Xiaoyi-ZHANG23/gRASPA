# challenge/ — GPU-free reproduction

`parser_under_test.cpp` is a standalone, GPU-free extract of gRASPA's `force_field.def` override
parser (the first pass of `OverWrite_Mixing_Rule()` in `src_clean/read_data.cpp`). It reproduces the
parsing behavior of the full GPU run, so you can investigate and fix the problem in seconds.

## Build & run
```bash
g++ -O2 -std=c++17 parser_under_test.cpp -o parser_under_test

# The reported case (CO2 in NaX zeolite):
./parser_under_test ../input/force_field.def

# Two control inputs for comparison / no-regression checking:
./parser_under_test ff_control_mof.def     # a single-marker, MOF-style force field with overrides
./parser_under_test ff_control_zero.def    # a single-marker force field with no overrides
```

## Output format
```
Ndefint=<n>
Nmixrule=<n>            # how many pairwise overrides will be applied
overrides_applied=<n>
<one normalized override line per applied entry>   # e.g.  C_co2 Na lennard-jones 362.292 3.320
```

`Nmixrule` / `overrides_applied` is the number of `# mixing rules to overwrite` pairwise entries the
parser decided to apply. Compare what you get on `../input/force_field.def` against what that file
actually contains, and against the control inputs.

## What "fixed" means
After your fix, the NaX case must apply its pairwise overrides (the strong CO₂–Na⁺ pairs must be
among them), **and** the control inputs must be unchanged. Keep the fix general — see the ground
rules in `../TASK.md`.
