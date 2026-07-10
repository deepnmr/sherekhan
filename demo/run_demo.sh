#!/bin/bash
# Generate the global / individual / mixed scenarios, run the full ShereKhan
# fit + per-residue comparison + jackknife on each, and evaluate the tool's
# verdict against ground truth.
set -e

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
ROOT_DIR="$(dirname "$SCRIPT_DIR")"
PY="${PYTHON:-python3}"

echo "Generating scenarios..."
"$PY" "$SCRIPT_DIR/make_scenarios.py"

for s in global individual mixed; do
    echo "Running $s..."
    ( cd "$SCRIPT_DIR/$s" && "$PY" "$ROOT_DIR/sk_run.py" "$s.conf" > out.txt 2>/dev/null )
done

echo "Running mixed recovery (outlier R5 removed)..."
( cd "$SCRIPT_DIR/mixed" && "$PY" "$ROOT_DIR/sk_run.py" mixed-recover.conf > out-recover.txt 2>/dev/null )

echo
"$PY" "$SCRIPT_DIR/evaluate.py"
