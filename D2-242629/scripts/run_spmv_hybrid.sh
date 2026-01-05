#!/bin/bash


# Compile (not timed)
mpicc -O3 -fopenmp ./src/spmv_mpi_parallelIO.c -o ./src/spmv_mpi_parallelIO

EXEC=./src/spmv_mpi_parallelIO
MATDIR=./matrices

export OMP_PROC_BIND=close
export OMP_PLACES=cores

configs=(
  "128 1"
  "64 2"
  "32 4"
  "16 8"
  "8 16"
)
echo "********************************"
echo "Hybrid MPI+OpenMP Strong scaling"
echo "********************************"

for mat in "$MATDIR/webbase-1M.mtx" \
           "$MATDIR/Flan_1565.mtx" \
           "$MATDIR/bcsstk18.mtx"; do
  echo
  echo "Matrix: $mat"
  echo "--------------------------------"
  for cfg in "${configs[@]}"; do
    set -- $cfg
    export OMP_NUM_THREADS=$2
    echo
    echo "MPI=$1 | OMP=$2"
    mpirun -np $1 $EXEC $mat
  done
done
