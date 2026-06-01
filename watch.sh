#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VENV_PYTHON="${ROOT}/.venv/bin/python"
SCRIPT="${ROOT}/printableFiles/hexagon.py"

if [[ ! -x "${VENV_PYTHON}" ]]; then
  echo "Virtual env not found. Run: ${ROOT}/setup.sh" >&2
  exit 1
fi

if ! command -v entr >/dev/null 2>&1; then
  echo "entr is required (e.g. brew install entr)" >&2
  exit 1
fi

echo "Watching printableFiles/*.py — runs ${SCRIPT} on change (Ctrl+C to stop)"
ls "${ROOT}"/printableFiles/*.py | entr "${VENV_PYTHON}" "${SCRIPT}"
