#!/bin/bash
# sweep_compare.sh — sweep a comparator over two directories of run outputs (base vs feature).
#
# The pre-merge gate for NEW FEATURES: prove a change altered ONLY the cases it was
# supposed to alter. Pairs <case>/output.txt by case name and judges each pair.
#
# Usage:
#   bash sweep_compare.sh [--diff] BASE_DIR TEST_DIR [EXPECTED_DIFF_FILE]
#
#   BASE_DIR/, TEST_DIR/   each contains <case>/output.txt (case = subdirectory name)
#   EXPECTED_DIFF_FILE     optional manifest: one case name per line ('#' comments OK;
#                          names must not contain '#'). Lists every case the change is
#                          SUPPOSED to alter. Omitted/empty = "nothing may change".
#   --diff                 judge by byte-identity (cmp) instead of score.py. Use this
#                          for GPU-free harness traces or any non-gRASPA output that
#                          score.py cannot parse.
#
# Default (score.py) mode judges by COMPARISON-class failures only (absolute-class
# entries — "Energy drift", "GPU drift", "Structure factor" — are reported but never
# fail a case: they are properties of a single run and fire even on a self-compare).
# NOTE score.py checks the REFERENCE's quantities against the test run — a quantity
# present ONLY in the test output is not flagged. For strict symmetric judgement
# (harness traces, text reports) use --diff, which is byte-exact in both directions.
#
# Per-case verdicts:
#   CONSTANT          no (comparison-class) difference, not in manifest   -> OK
#   EXPECTED-DIFF     differs and is listed in the manifest               -> OK
#   UNEXPECTED-DIFF   differs but is NOT in the manifest                  -> FAIL (regression)
#   EXPECTED-MISSING  in the manifest but did NOT differ                  -> FAIL (feature inert)
#   ERROR             an output could not be parsed/read                  -> FAIL
# A case that cannot be compared (missing dir or output.txt on either side) is
# SKIPPED — a warning, NOT a pass — unless it is declared in the manifest, in which
# case it is a FAIL (a declared diff that was never verified must not pass the gate).
# Comparing zero pairs is also a FAIL (wrong directories?).
#
# Exit: 0 = surgical (every case OK), 1 = at least one FAIL / nothing compared,
#       2 = bad invocation.

set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCORE="$HERE/score.py"

usage(){ sed -n '2,/bad invocation/p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//' >&2; exit 2; }

DIFFMODE=0
if [ "${1:-}" = "--diff" ]; then DIFFMODE=1; shift; fi
[ $# -ge 2 ] && [ $# -le 3 ] || usage
BASE="${1%/}"; TEST="${2%/}"; MANIFEST="${3:-}"
[ -d "$BASE" ] || { echo "ERROR: base dir not found: $BASE" >&2; exit 2; }
[ -d "$TEST" ] || { echo "ERROR: test dir not found: $TEST" >&2; exit 2; }
if [ -n "$MANIFEST" ] && [ ! -f "$MANIFEST" ]; then
  echo "ERROR: manifest not found: $MANIFEST" >&2; exit 2
fi

TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT

# expected-diff set: strip '#' comments, trim leading/trailing whitespace only
declare -A EXPECT=(); nexpect=0
if [ -n "$MANIFEST" ]; then
  while IFS= read -r line || [ -n "$line" ]; do
    line="${line%%#*}"
    line="${line#"${line%%[![:space:]]*}"}"
    line="${line%"${line##*[![:space:]]}"}"
    [ -n "$line" ] && { EXPECT["$line"]=1; nexpect=$((nexpect+1)); }
  done < "$MANIFEST"
fi

# classify one score.py JSON: prints "<n_comparison_failures>" then indented details
classify(){
  python3 -c '
import json, sys
ABS = ("Energy drift", "GPU drift", "Structure factor")
try:
    j = json.load(open(sys.argv[1]))
except Exception:
    print(-1); sys.exit(0)
comp = [f for f in j.get("failures", []) if not f.startswith(ABS)]
print(len(comp))
for f in comp[:8]:
    print("        " + f)
if len(comp) > 8:
    print(f"        ... and {len(comp)-8} more")
' "$1"
}

pass=1; n_ok=0; n_fail=0; n_skip=0
declare -A SEEN=()

# skip_or_fail <case> <reason> — a manifest-declared case may never be silently skipped
skip_or_fail(){
  if [ -n "${EXPECT[$1]:-}" ]; then
    printf '  FAIL  %-34s declared in manifest but %s — diff never verified\n' "$1" "$2"
    pass=0; n_fail=$((n_fail+1))
  else
    printf '  SKIP  %-34s %s\n' "$1" "$2"
    n_skip=$((n_skip+1))
  fi
}

echo "== sweep_compare: $TEST vs baseline $BASE =="
[ "$DIFFMODE" -eq 1 ] && echo "   mode: --diff (byte-exact)"
[ -n "$MANIFEST" ] && echo "   expected-diff manifest: $MANIFEST ($nexpect case(s))"

for bdir in "$BASE"/*/; do
  [ -d "$bdir" ] || continue
  case="$(basename "$bdir")"; SEEN["$case"]=1
  bout="${bdir%/}/output.txt"; tout="$TEST/$case/output.txt"
  if [ ! -f "$bout" ]; then skip_or_fail "$case" "no output.txt in base dir"; continue; fi
  if [ ! -f "$tout" ]; then skip_or_fail "$case" "missing in test dir"; continue; fi

  if [ "$DIFFMODE" -eq 1 ]; then
    if [ ! -r "$bout" ] || [ ! -r "$tout" ]; then
      printf '  FAIL  %-34s ERROR: output not readable\n' "$case"
      pass=0; n_fail=$((n_fail+1)); continue
    fi
    if cmp -s "$bout" "$tout"; then ncomp=0; detail=""
    else
      ncomp=1
      detail="$(diff "$bout" "$tout" 2>/dev/null | head -6 | sed 's/^/        /')"
    fi
  else
    # score.py validates only the TEST side; guard the baseline ourselves
    if ! grep -q "Work took" "$bout" 2>/dev/null; then
      printf '  FAIL  %-34s ERROR: baseline output unparseable (no end-of-run summary — was it captured with 2>&1?)\n' "$case"
      pass=0; n_fail=$((n_fail+1)); continue
    fi
    python3 "$SCORE" "$tout" "$bout" > "$TMP/s.json" 2>"$TMP/s.err"; rc=$?
    if [ "$rc" -eq 2 ]; then
      printf '  FAIL  %-34s ERROR: score.py could not parse (exit 2)\n' "$case"
      sed 's/^/        /' "$TMP/s.err" | head -3
      pass=0; n_fail=$((n_fail+1)); continue
    fi
    out="$(classify "$TMP/s.json")"
    ncomp="$(printf '%s\n' "$out" | head -1)"
    detail="$(printf '%s\n' "$out" | tail -n +2)"
    if [ "$ncomp" = "-1" ]; then
      printf '  FAIL  %-34s ERROR: unreadable score.py output\n' "$case"
      pass=0; n_fail=$((n_fail+1)); continue
    fi
  fi

  if [ "$ncomp" -gt 0 ]; then
    if [ -n "${EXPECT[$case]:-}" ]; then
      printf '  OK    %-34s EXPECTED-DIFF\n' "$case"
      n_ok=$((n_ok+1))
    else
      printf '  FAIL  %-34s UNEXPECTED-DIFF:\n' "$case"
      [ -n "$detail" ] && printf '%s\n' "$detail"
      pass=0; n_fail=$((n_fail+1))
    fi
  else
    if [ -n "${EXPECT[$case]:-}" ]; then
      printf '  FAIL  %-34s EXPECTED-MISSING (manifest says it should differ; it did not)\n' "$case"
      pass=0; n_fail=$((n_fail+1))
    else
      printf '  OK    %-34s CONSTANT\n' "$case"
      n_ok=$((n_ok+1))
    fi
  fi
done

# cases only in the test dir (e.g. a new regression case added by the feature)
for tdir in "$TEST"/*/; do
  [ -d "$tdir" ] || continue
  case="$(basename "$tdir")"
  if [ -z "${SEEN[$case]:-}" ]; then
    if [ -f "$tdir/output.txt" ]; then
      skip_or_fail "$case" "new in test dir (no baseline — commit a vetted reference)"
    else
      skip_or_fail "$case" "test dir entry without output.txt"
    fi
    SEEN["$case"]=1
  fi
done

# manifest entries that never appeared at all (guarded expansion: safe on bash <= 4.3)
for case in ${EXPECT[@]+"${!EXPECT[@]}"}; do
  if [ -z "${SEEN[$case]:-}" ]; then
    printf '  FAIL  %-34s in manifest but no such case in either dir\n' "$case"
    pass=0; n_fail=$((n_fail+1))
  fi
done

echo
echo "Summary: $n_ok OK, $n_fail FAIL, $n_skip skipped."
if [ $((n_ok + n_fail)) -eq 0 ]; then
  echo "RESULT: NOTHING COMPARED — no <case>/output.txt pairs found. Wrong directories?"
  exit 1
fi
if [ "$pass" -eq 1 ]; then
  echo "RESULT: SURGICAL — the change alters exactly what the manifest declares."
else
  echo "RESULT: NOT SURGICAL — see FAIL lines above."; exit 1
fi
