#!/bin/bash
# grade.sh — automated grader for the gRASPA Bug #4 debugging challenge.
#
# Compiles a candidate parser_under_test.cpp and checks its OBSERVABLE behavior
# against the trigger and two control force fields. Uses the grader's OWN copies
# of the .def inputs (solution/grade_*.def), so a candidate cannot pass by editing
# the input files instead of fixing the code.
#
# Usage:
#   bash grade.sh [path/to/parser_under_test.cpp]
# Default candidate: ../challenge/parser_under_test.cpp
#
# Exit: 0 = PASS (correct AND surgical), 1 = FAIL.

set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CAND="${1:-$HERE/../challenge/parser_under_test.cpp}"
TRIG="$HERE/grade_trigger.def"
CMOF="$HERE/grade_control_mof.def"
CZERO="$HERE/grade_control_zero.def"
CADV="$HERE/grade_control_adv.def"   # held-out adversarial input (not given to the agent)

TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT
bin="$TMP/cand"

echo "== Candidate: $CAND =="
if [ ! -f "$CAND" ]; then echo "FAIL: candidate file not found"; exit 1; fi
if ! g++ -O2 -std=c++17 "$CAND" -o "$bin" 2>"$TMP/cc.log"; then
  echo "FAIL: candidate does not compile:"; sed 's/^/    /' "$TMP/cc.log"; exit 1
fi

nmix(){ "$bin" "$1" 2>/dev/null | sed -n 's/^Nmixrule=//p' | head -1; }
out(){  "$bin" "$1" 2>/dev/null; }

pass=1
check(){ # name expected actual
  if [ "$2" = "$3" ]; then printf "  PASS  %-26s %s\n" "$1" "$3";
  else printf "  FAIL  %-26s got '%s' expected '%s'\n" "$1" "$3" "$2"; pass=0; fi
}

echo "== Checks =="
# 1. Trigger (the NaX zeolite case): all 11 pairwise overrides must be applied.
check "trigger Nmixrule"      "11" "$(nmix "$TRIG")"

# 2. The strong CO2-Na+ pairs must actually be among the applied overrides.
if out "$TRIG" | grep -q "C_co2 Na" && out "$TRIG" | grep -q "O_co2 Na"; then
  printf "  PASS  %-26s %s\n" "trigger CO2-Na present" "yes"
else
  printf "  FAIL  %-26s %s\n" "trigger CO2-Na present" "missing"; pass=0
fi

# 3. Single-marker MOF-like control with overrides must be UNCHANGED (= 3).
check "control_mof Nmixrule"  "3"  "$(nmix "$CMOF")"

# 4. Single-marker control with zero overrides must stay 0 (no spurious overrides).
check "control_zero Nmixrule" "0"  "$(nmix "$CZERO")"

# 5. Held-out adversarial control: a LEADING zero-count marker followed by a real
#    block. First-occurrence-wins (the general fix) selects the leading marker -> 0.
#    A non-general "ignore a marker whose count is 0" cheat selects the later block
#    -> 2, so this check fails it. (The agent never sees this input.)
check "adversarial Nmixrule"  "0"  "$(nmix "$CADV")"

# 6. Whole-suite no-regression: the fix must change ONLY the trigger example across all of Examples/
#    (skipped automatically if Examples/ isn't reachable, e.g. when grading a sandbox copy).
echo "== Examples no-regression (all force_field.def) =="
if bash "$HERE/regression_examples.sh" "$CAND" > "$TMP/reg.log" 2>&1; then
  sed -n 's/^/  /p' "$TMP/reg.log" | grep -E 'scanned=|RESULT|SKIP' || true
  grep -q '^SKIP' "$TMP/reg.log" && echo "  (Examples/ not reachable — suite check skipped)"
else
  sed 's/^/  /' "$TMP/reg.log"; pass=0
fi

echo
if [ "$pass" = "1" ]; then
  echo "RESULT: PASS — overrides applied correctly on the trigger, normal path preserved, no regressions across Examples/."
  exit 0
else
  echo "RESULT: FAIL — see checks above. (Buggy/unfixed code drops the trigger's overrides → Nmixrule=0.)"
  exit 1
fi
