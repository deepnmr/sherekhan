#!/bin/bash
set -e

# Get the directory of this script
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
ROOT_DIR="$(dirname "$SCRIPT_DIR")"

echo "Running tests..."

# Update PATH to include root dir where scripts are
export PATH=$ROOT_DIR:$PATH

echo "Root dir: $ROOT_DIR"
echo "Script dir: $SCRIPT_DIR"

echo "Running unit self-checks..."
for t in "$SCRIPT_DIR"/test_*.py; do
    echo "  $(basename "$t")"
    "${PYTHON:-python3}" "$t"
done

cd "$SCRIPT_DIR/fast"
echo "Running fast tests..."
"${PYTHON:-python3}" ../../sk_run.py fast.conf
"${PYTHON:-python3}" ../../sk_run.py fast-matrix.conf
"${PYTHON:-python3}" ../../sk_run.py fast-globalfit.conf   # individual-vs-global + jackknife

cd "$SCRIPT_DIR/slow"
echo "Running slow tests..."
"${PYTHON:-python3}" ../../sk_run.py slow.conf
"${PYTHON:-python3}" ../../sk_run.py slow-matrix.conf

echo "Test execution finished."
