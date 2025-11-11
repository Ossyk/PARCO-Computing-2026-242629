#!/bin/bash
# run.sh - Linux shell script to run spMV_parall
# This script is in: project/scripts/

# Get the script's directory and go to project root
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR/.."

echo "Usage $0 <threads> <scheduler> <chunks> <matrix>"


THREADS="$1"
SCHEDULER="$2"
CHUNKS="$3"
MATRIX="$4"
EXE="src/spMV_parall"


if [ $# -ne 4 ]; then
    echo "Usage: $0 <threads> <scheduler> <chunks> <matrix_path>"
    echo "Example: $0 4 dynamic 100 matrices/m4.csr"
    exit 1
fi

if [[ "$SCHEDULER" != "static" && "$SCHEDULER" != "dynamic" && "$SCHEDULER" != "guided" && "$SCHEDULER" != "auto" ]]; then
    echo "Error: scheduler must be one of {static, dynamic, guided, auto}."
    exit 1
fi

if [ ! -f "$MATRIX" ]; then
    echo "Error: Matrix file not found: $MATRIX"
    exit 1
fi

# Check if executable exists
if [ ! -f "$EXE" ]; then
    echo "Error: Executable not found at $EXE"
    echo "Please compile first with: gcc -fopenmp -std=c99 -o src/spMV_parall src/spMV_parall.c"
    exit 1
fi



echo "==================================="
echo "Running SpMV OpenMP Benchmarks"
echo "Matrix: $MATRIX"
echo "==================================="
echo

# Test different thread counts with static schedule
echo "--- Testing thread scaling ---"
"$EXE" 1 "$SCHEDULER" "$CHUNKS" "$MATRIX"
"$EXE" 2 "$SCHEDULER" "$CHUNKS" "$MATRIX"
"$EXE" 4 "$SCHEDULER" "$CHUNKS" "$MATRIX"
"$EXE" 8 "$SCHEDULER" "$CHUNKS" "$MATRIX"
echo

# Test different schedules 
echo "--- Testing schedules ---"
"$EXE" "$THREADS" static "$CHUNKS" "$MATRIX"
"$EXE" "$THREADS" dynamic "$CHUNKS" "$MATRIX"
"$EXE" "$THREADS" guided "$CHUNKS" "$MATRIX"
"$EXE" "$THREADS" auto 0 "$MATRIX"
echo

# Test different chunk sizes 
echo "--- Testing chunk sizes ---"
"$EXE" "$THREADS" "$SCHEDULER" 1 "$MATRIX"
"$EXE" "$THREADS" "$SCHEDULER" 10 "$MATRIX"
"$EXE" "$THREADS" "$SCHEDULER" 100 "$MATRIX"
"$EXE" "$THREADS" "$SCHEDULER" 1000 "$MATRIX"

echo

echo "==================================="
echo "Benchmarks Complete"
echo "==================================="


