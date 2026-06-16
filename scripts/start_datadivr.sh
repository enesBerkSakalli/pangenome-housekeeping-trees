#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATADIVR_DIR="$ROOT_DIR/external/DataDiVR_WebApp"
PORT="${PORT:-3000}"

if [[ ! -f "$DATADIVR_DIR/app.py" ]]; then
  echo "DataDiVR checkout not found at $DATADIVR_DIR" >&2
  echo "Clone https://github.com/menchelab/DataDiVR_WebApp into external/DataDiVR_WebApp first." >&2
  exit 1
fi

if [[ -x "$DATADIVR_DIR/venv/bin/python" ]]; then
  PYTHON="$DATADIVR_DIR/venv/bin/python"
elif [[ -x "$ROOT_DIR/.venv/bin/python" ]]; then
  PYTHON="$ROOT_DIR/.venv/bin/python"
else
  PYTHON="python3"
fi

cd "$DATADIVR_DIR"
export FLASK_APP=app.py
export FLASK_DEBUG="${FLASK_DEBUG:-0}"
export MPLCONFIGDIR="${MPLCONFIGDIR:-${TMPDIR:-/tmp}/pangenome-matplotlib}"
mkdir -p "$MPLCONFIGDIR"

exec "$PYTHON" -m flask run --with-threads --host=127.0.0.1 --port="$PORT"
