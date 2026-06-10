# FIELD_TEST.md — the new-feature workflow, field-tested end to end

On 2026-06-09 the "Validating a NEW feature" workflow was field-tested by an AI agent that knew
the symptom-free task only: *add a mock `Fugacity <Pa>` keyword to `simulation.input` (input =
gas-phase fugacity, internally converted to the pressure gRASPA uses; `Fugacity` overrides
`Pressure` when both are present), then validate it with the kit.* An independent adversarial
grader then re-ran every claim from scratch. Verdict: **PASS** — the workflow caught the real
pitfalls. The mock feature itself was never merged; this file records what the exercise proved.

## Why `Fugacity` is a treacherous keyword (and why the workflow caught it)

gRASPA matches input keywords by **substring** (`str.find`), and `FugacityCoefficient` — present
in **26 of 29** `Examples/*/simulation.input` files — *contains* `Fugacity`. So the obvious
house-style implementation

```cpp
if (str.find("Fugacity", 0) != std::string::npos)  // WRONG
```

fires on every `FugacityCoefficient` line:

| Input line | Naive result |
|---|---|
| `FugacityCoefficient 1.0` | pressure silently becomes **1 Pa** (silent corruption) |
| `FugacityCoefficient PR-EOS` | `std::stod` throws → crash |

## What the workflow produced

Following the skill's steps (manifest first → regression case → standalone reader harness →
all-Examples sweep, gated like `sweep_compare.sh`):

- **Correct implementation:** match the keyword as an **exact first token** (tokenize the line,
  compare `tokens[0] == "Fugacity"`), apply the override after the read loop so keyword order
  doesn't matter.
- **Gate result over all 44 `Examples/**/simulation.input`:** token-exact version =
  **43 CONSTANT + 1 EXPECTED-DIFF** (only the new regression case changes — surgical, exit 0).
  Naive substring version = **38/44 FAIL** in exactly the two flavors above.
- The grader independently rebuilt the harnesses and reproduced both sweeps byte-identically,
  and confirmed the token-exact guard lives in the real `read_data.cpp` diff, not just the
  harness.

## What the exercise added to this kit

1. `sweep_compare.sh --diff` — byte-exact gate mode, because `score.py` can only parse full
   gRASPA outputs, not harness traces.
2. The **token-exact first-token rule** as the prescribed fix for keyword collisions (skill
   step 6).
3. A warning that **old binaries silently swallow new keywords** — `Check_Inputs_In_read_data_cpp`
   "validates" keywords by substring-grepping the *source*, so an old build can falsely accept a
   new keyword and quietly run with the old behavior. New keywords must document their minimum
   gRASPA version.
4. Caveats that the parse-level sweep proves the *reader* is surgical (not the physics), and that
   a hand-mirrored harness must be re-diffed against production.

## Bonus pre-existing finding

`Examples/Simplify_Force_Field_Files/simulation.input` contains a placeholder line
`Pressure pres` on which a production `std::stod` would throw. Behavior is identical on base and
feature builds (so it's outside any feature gate), but worth knowing.
