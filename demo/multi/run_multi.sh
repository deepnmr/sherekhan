#!/bin/bash
# Generate the N-individual (m1..m4) and K-global-group (g2, g3) scenarios,
# run the CLI single-global fit on each, then evaluate identification of the
# individuals (Family A) and the K-global structure selection (Family B).
set -e

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
ROOT_DIR="$(dirname "$(dirname "$SCRIPT_DIR")")"
PY="${PYTHON:-python3}"

echo "Generating scenarios..."
"$PY" "$SCRIPT_DIR/make_multi.py"

for s in m1 m2 m3 m4 g2 g3; do
    echo "Running $s..."
    ( cd "$SCRIPT_DIR/$s" && "$PY" "$ROOT_DIR/sk_run.py" "$s.conf" > out.txt 2>/dev/null )
done

echo
"$PY" "$SCRIPT_DIR/evaluate_multi.py"
