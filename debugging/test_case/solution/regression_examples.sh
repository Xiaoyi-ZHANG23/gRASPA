#!/bin/bash
# regression_examples.sh — GPU-free no-regression check across the WHOLE Examples/ suite.
#
# Runs the parser on every Examples/*/force_field.def and confirms a candidate fix changes the
# override-parsing behavior on ONLY the trigger case (CO2_NaX_Zeolite) — i.e. it "did not affect
# other parts". For every other example the candidate must match the unfixed (buggy) baseline.
#
# Usage:  bash regression_examples.sh [candidate.cpp]
# Default candidate: ../challenge/parser_under_test.cpp (after the agent has fixed it).
# (Use solution/parser_fixed.cpp to see the reference pass.)
#
# Exit: 0 = surgical (only the trigger changed, and it changed correctly), 1 = regression / wrong.

set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CAND="${1:-$HERE/../challenge/parser_under_test.cpp}"
BUGGY="$HERE/parser_buggy_baseline.cpp"            # pristine unfixed baseline
EXAMPLES="$(cd "$HERE/../../.." && pwd)/Examples"  # repo Examples/ (test_case/solution -> repo root)
EXPECTED_CHANGED="CO2_NaX_Zeolite"                 # the only example that should change

TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT
cand="$TMP/cand"; bug="$TMP/bug"

if [ ! -d "$EXAMPLES" ]; then
  echo "SKIP: Examples/ not found at $EXAMPLES (run from inside the repo to enable this check)."; exit 0
fi
echo "== Compiling candidate + buggy baseline =="
g++ -O2 -std=c++17 "$CAND"  -o "$cand" 2>"$TMP/c.log" || { echo "FAIL: candidate does not compile"; sed 's/^/    /' "$TMP/c.log"; exit 1; }
g++ -O2 -std=c++17 "$BUGGY" -o "$bug"  2>/dev/null   || { echo "FAIL: buggy baseline does not compile (bug in test_case)"; exit 1; }

nmix(){ "$1" "$2" 2>/dev/null | sed -n 's/^Nmixrule=//p' | head -1; }

total=0; changed=0; unexpected=0; trigger_seen=0; trigger_ok=0
changed_list=""
while IFS= read -r f; do
  total=$((total+1))
  rel="${f#"$EXAMPLES"/}"
  b="$(nmix "$bug" "$f")"; c="$(nmix "$cand" "$f")"
  if [ "$b" != "$c" ]; then
    changed=$((changed+1)); changed_list="$changed_list $rel"
    case "$f" in
      *"/$EXPECTED_CHANGED/"*) trigger_seen=1; [ "$c" = "11" ] && trigger_ok=1
        printf "  changed (expected)  buggy=%s candidate=%s  %s\n" "$b" "$c" "$rel" ;;
      *) unexpected=$((unexpected+1))
        printf "  CHANGED (UNEXPECTED) buggy=%s candidate=%s  %s\n" "$b" "$c" "$rel" ;;
    esac
  fi
done < <(find "$EXAMPLES" -name force_field.def | sort)

echo
echo "scanned=$total  unchanged=$((total-changed))  changed=$changed  unexpected=$unexpected"
ok=1
[ "$unexpected" -ne 0 ] && { echo "FAIL: the fix changed examples it should not have:$changed_list"; ok=0; }
[ "$trigger_seen" -ne 1 ] && { echo "FAIL: the trigger example ($EXPECTED_CHANGED) did not change — candidate didn't fix the bug."; ok=0; }
[ "$trigger_seen" -eq 1 ] && [ "$trigger_ok" -ne 1 ] && { echo "FAIL: trigger changed but candidate Nmixrule != 11 there."; ok=0; }

echo
if [ "$ok" = "1" ]; then
  echo "RESULT: PASS — surgical. Only $EXPECTED_CHANGED changed (0 -> 11); all $((total-1)) other examples unchanged."
  exit 0
else
  echo "RESULT: FAIL — see above."
  exit 1
fi
