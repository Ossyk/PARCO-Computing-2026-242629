#!/bin/bash


EXEC=./src/spmv_mpi_parallelIO
MATDIR=./matrices

echo "***********************************"
echo "         STRONG SCALING — 1D"
echo "***********************************"
matrices=(
  "$MATDIR/webbase-1M.mtx"
  "$MATDIR/Flan_1565.mtx"
  "$MATDIR/bcsstk18.mtx"
)
prs_strong=(1 2 4 8 16 32 64 128)

for mat in "${matrices[@]}"; do
  echo
  echo "Matrix: $mat"
  echo "-----------------------------------"
  for p in "${prs_strong[@]}"; do
    echo "MPI processes: $p"
    mpirun -np $p $EXEC $mat
    echo
  done
done

echo "***********************************"
echo "         WEAK SCALING — 1D"
echo "***********************************"
for i in 1 2 4 8 16 32; do
  name=$((i * 20))
  echo "MPI processes: $i"
  echo "------------------"
  
  perf stat mpirun -np $i $EXEC $MATDIR/randomN${name}k.mtx
  echo
done
echo "MPI processes: 64"
echo "------------------"
perf stat mpirun -np 64  $EXEC $MATDIR/random1M280k.mtx
echo "MPI processes: 128"
echo "------------------"
perf stat mpirun -np 128 $EXEC $MATDIR/random2M560k.mtx