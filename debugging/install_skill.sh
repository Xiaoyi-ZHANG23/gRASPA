#!/bin/bash
# install_skill.sh — install the graspa-debug skill for Claude Code (one command).
#
# The tracked copy of the skill lives at debugging/SKILL.md (this repo gitignores
# .claude/, so it cannot ship there directly). This script copies it to where
# Claude Code auto-discovers skills:
#
#   bash debugging/install_skill.sh             # user-level:    ~/.claude/skills/graspa-debug/
#   bash debugging/install_skill.sh --project   # this clone:    <repo>/.claude/skills/graspa-debug/
#
# User-level follows you across every project on this machine; project-level stays
# inside this clone (each person installs their own — .claude/ is never committed).
# Codex / other agents need NO install: the repo-root AGENTS.md already points them
# at debugging/DEBUGGING.md.
#
# Exit: 0 = installed, 1 = error.

set -euo pipefail
SELF="$(readlink -f "${BASH_SOURCE[0]}" 2>/dev/null || echo "${BASH_SOURCE[0]}")"
HERE="$(cd "$(dirname "$SELF")" && pwd)"   # debugging/ (symlink-safe)
SRC="$HERE/SKILL.md"
[ -f "$SRC" ] || { echo "ERROR: $SRC not found" >&2; exit 1; }
[ $# -le 1 ] || { echo "ERROR: too many arguments (use --user, --project, or --help)" >&2; exit 1; }

case "${1:---user}" in
  --user)    : "${HOME:?HOME is not set}"; DEST="$HOME/.claude/skills/graspa-debug" ;;
  --project) DEST="$(cd "$HERE/.." && pwd)/.claude/skills/graspa-debug" ;;
  -h|--help) sed -n '2,/^# Exit:/p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
  *) echo "ERROR: unknown option '$1' (use --user, --project, or --help)" >&2; exit 1 ;;
esac

mkdir -p "$DEST"
if [ -f "$DEST/SKILL.md" ] && cmp -s "$SRC" "$DEST/SKILL.md"; then
  echo "Already up to date: $DEST/SKILL.md"
else
  cp "$SRC" "$DEST/SKILL.md"
  echo "Installed: $DEST/SKILL.md"
fi
echo "Claude Code will auto-discover the 'graspa-debug' skill in your next session."
echo "(Sanity-check the kit itself anytime with: bash $HERE/selftest.sh)"
