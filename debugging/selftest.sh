#!/bin/bash
# selftest.sh — one-command health check of the GPU-free core of this debugging kit.
#
# Runs in ~1 minute, needs only g++ and python3 (no GPU, no gRASPA build). Checks:
#   1. the standalone parser repro reproduces Bug #4 (buggy pass drops all 11
#      overrides; fixed pass applies them) on Examples/CO2_NaX_Zeolite/
#   2. the test_case grader FAILs the unfixed challenge and PASSes the reference
#      fix (including the whole-Examples/ no-regression sweep)
#   3. score.py demonstrates its documented gotcha: a self-compare of the shipped
#      CO2_NaX_Zeolite output exits nonzero on ABSOLUTE-class drift only, with
#      zero comparison-class failures — judge A/B equality by the failures array,
#      not the exit code.
#
# Usage:  bash debugging/selftest.sh        (from anywhere; paths are self-anchored)
# Exit:   0 = kit is healthy on this machine, 1 = something is broken.

set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"   # debugging/
REPO="$(cd "$HERE/.." && pwd)"
TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT

pass=1
ok(){  printf '  PASS  %s\n' "$1"; }
bad(){ printf '  FAIL  %s\n' "$1"; pass=0; }

echo "== 0. toolchain =="
for t in g++ python3; do
  if command -v "$t" >/dev/null 2>&1; then ok "$t found"; else bad "$t not found"; fi
done
[ "$pass" -eq 1 ] || { echo "RESULT: FAIL — missing tools."; exit 1; }

echo "== 1. GPU-free parser repro (Bug #4) =="
if g++ -O2 -std=c++17 "$HERE/repro/parse_repro2.cpp" -o "$TMP/pr2" 2>"$TMP/cc.log"; then
  ok "repro/parse_repro2.cpp compiles"
  out="$("$TMP/pr2" "$REPO/Examples/CO2_NaX_Zeolite/force_field.def")"
  printf '%s\n' "$out" | sed 's/^/      /'
  if printf '%s\n' "$out" | grep -q 'NEW(buggy):.*Nmixrule= 0.*WRONG'; then
    ok "buggy pass drops all overrides (Nmixrule=0)"
  else bad "buggy pass did not reproduce the bug"; fi
  if printf '%s\n' "$out" | grep -q 'NEW+FIX:.*Nmixrule=11.*CORRECT'; then
    ok "fixed pass applies all 11 overrides"
  else bad "fixed pass gave the wrong answer"; fi
else
  bad "repro/parse_repro2.cpp does not compile:"; sed 's/^/      /' "$TMP/cc.log"
fi

echo "== 2. test_case grader (must FAIL unfixed / PASS reference fix) =="
if bash "$HERE/test_case/solution/grade.sh" "$HERE/test_case/challenge/parser_under_test.cpp" >"$TMP/g1.log" 2>&1; then
  bad "grader PASSed the unfixed challenge (it should FAIL)"
else
  ok "grader FAILs the unfixed challenge"
fi
if bash "$HERE/test_case/solution/grade.sh" "$HERE/test_case/solution/parser_fixed.cpp" >"$TMP/g2.log" 2>&1; then
  ok "grader PASSes the reference fix (incl. whole-Examples/ no-regression sweep)"
else
  bad "grader FAILed the reference fix:"; tail -15 "$TMP/g2.log" | sed 's/^/      /'
fi

echo "== 3. score.py gotcha (self-compare: absolute-class drift only) =="
ref="$REPO/Examples/CO2_NaX_Zeolite/output.txt"
python3 "$HERE/score.py" "$ref" "$ref" >"$TMP/score.json" 2>"$TMP/score.err"; rc=$?
absonly="$(python3 -c '
import json,sys
f = json.load(open(sys.argv[1])).get("failures", [])
print(1 if f and all(x.startswith("Energy drift") for x in f) else 0)
' "$TMP/score.json" 2>/dev/null || echo 0)"
if [ "$rc" -eq 1 ] && [ "$absonly" = "1" ]; then
  ok "self-compare exits 1 with ONLY absolute-class (Energy drift) failures"
  echo "      (this is the documented gotcha: read the JSON failures array, not the exit code)"
else
  bad "score.py self-compare behaved unexpectedly (exit=$rc):"
  head -20 "$TMP/score.json" 2>/dev/null | sed 's/^/      /'
fi

echo "== 4. sweep_compare.sh (pre-merge feature gate) =="
mkdir -p "$TMP/base/caseA" "$TMP/feat/caseA" "$TMP/base/caseB" "$TMP/feat/caseB"
cp "$ref" "$TMP/base/caseA/output.txt"; cp "$ref" "$TMP/feat/caseA/output.txt"
cp "$ref" "$TMP/base/caseB/output.txt"
sed 's/PseudoAtom Type: Al\[0\], #: 55/PseudoAtom Type: Al[0], #: 54/' "$ref" > "$TMP/feat/caseB/output.txt"
if bash "$HERE/sweep_compare.sh" "$TMP/base" "$TMP/feat" >"$TMP/sw1.log" 2>&1; then
  bad "sweep with an injected diff (caseB) PASSed — it should FAIL"
else
  if grep -q 'UNEXPECTED-DIFF' "$TMP/sw1.log" && grep -q 'caseA *CONSTANT' "$TMP/sw1.log"; then
    ok "no-manifest sweep: caseA CONSTANT, injected caseB flagged UNEXPECTED-DIFF"
  else
    bad "no-manifest sweep gave unexpected verdicts:"; sed 's/^/      /' "$TMP/sw1.log"
  fi
fi
echo "caseB" > "$TMP/manifest.txt"
if bash "$HERE/sweep_compare.sh" "$TMP/base" "$TMP/feat" "$TMP/manifest.txt" >"$TMP/sw2.log" 2>&1; then
  ok "manifest sweep: declared diff (caseB) accepted as EXPECTED-DIFF — surgical"
else
  bad "manifest sweep should PASS when the diff is declared:"; sed 's/^/      /' "$TMP/sw2.log"
fi
mkdir -p "$TMP/db/caseT" "$TMP/dt/caseT"
printf 'branch=Pressure P=10000\n' > "$TMP/db/caseT/output.txt"
printf 'branch=Fugacity P=25000\n' > "$TMP/dt/caseT/output.txt"
if bash "$HERE/sweep_compare.sh" --diff "$TMP/db" "$TMP/dt" >"$TMP/sw3.log" 2>&1; then
  bad "--diff mode missed a byte-level trace difference"
else
  ok "--diff mode (harness traces): byte diff flagged UNEXPECTED-DIFF"
fi

echo
if [ "$pass" -eq 1 ]; then
  echo "RESULT: PASS — debugging kit is healthy on this machine."
else
  echo "RESULT: FAIL — see checks above."; exit 1
fi
