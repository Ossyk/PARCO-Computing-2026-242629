#!/bin/bash
# run.sh - Linux shell script to run spMV_parall
# This script is in: project/scripts/

# Get the script's directory and go to project root
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR/.."

# Check if matrix file is provided
if [ $# -eq 0 ]; then
    echo "Usage: scripts/run.sh <matrix.csr>"
    echo "Example: scripts/run.sh matrices/m1.csr"
    exit 1
fi

if [ $# -ne 1 ]; then
    echo "Usage: $0 <matrix_path>"
    echo "Example: $0 matrices/m4.csr"
    exit 1
fi

MATRIX="$1"
EXE="src/spMV_parall"

# Check if executable exists
if [ ! -f "$EXE" ]; then
    echo "Error: Executable not found at $EXE"
    echo "Please compile first with: gcc -fopenmp -o src/spMV_parall src/spMV_parall.c"
    exit 1
fi

# Check if matrix file exists
if [ ! -f "$MATRIX" ]; then
    echo "Error: Matrix file not found: $MATRIX"
    exit 1
fi

echo "==================================="
echo "Running SpMV OpenMP Benchmarks"
echo "Matrix: $MATRIX"
echo "==================================="
echo

# Test different thread counts with static schedule
echo "--- Testing thread scaling (static, chunk=100) ---"
"$EXE" 1 static 100 "$MATRIX"
"$EXE" 2 static 100 "$MATRIX"
"$EXE" 4 static 100 "$MATRIX"
"$EXE" 8 static 100 "$MATRIX"
"$EXE" 16 static 100 "$MATRIX"
"$EXE" 32 static 100 "$MATRIX"
echo

# Test different schedules with 4 threads
echo "--- Testing schedules (4 threads, chunk=100) ---"
"$EXE" 4 static 100 "$MATRIX"
"$EXE" 4 dynamic 100 "$MATRIX"
"$EXE" 4 guided 100 "$MATRIX"
"$EXE" 4 auto 0 "$MATRIX"
echo

# Test different chunk sizes with dynamic schedule
echo "--- Testing chunk sizes (4 threads, dynamic) ---"
"$EXE" 4 dynamic 1 "$MATRIX"
"$EXE" 4 dynamic 10 "$MATRIX"
"$EXE" 4 dynamic 100 "$MATRIX"
"$EXE" 4 dynamic 1000 "$MATRIX"
"$EXE" 4 dynamic 0 "$MATRIX"
echo

echo "==================================="
echo "Benchmarks Complete"
echo "==================================="