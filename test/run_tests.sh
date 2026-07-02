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

cd "$SCRIPT_DIR/fast"
echo "Running fast tests..."
python3 ../../sk_run.py fast.conf
python3 ../../sk_run.py fast-matrix.conf

cd "$SCRIPT_DIR/slow"
echo "Running slow tests..."
python3 ../../sk_run.py slow.conf
python3 ../../sk_run.py slow-matrix.conf

echo "Test execution finished."
