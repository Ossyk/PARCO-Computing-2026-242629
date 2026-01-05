#!/bin/bash

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 matrix.mtx"
  exit 1
fi

matrix=$1



echo "*************************"
echo "1D Strong scaling"
echo "*************************"

prs=("1" "2" "4" "8" "16" "32" "64" "128") 
for p in "${prs[@]}"; do
  echo "Nb of processes : $p" 
  mpirun -np $p ./src/spmv_mpi_parallelIO ./matrices/$matrix
  echo 
done

echo "*************************"
echo "2D Strong scaling"
echo "*************************"

prs=("1" "2" "4" "8" "16" "32" "64" "128") 
for p in "${prs[@]}"; do
  echo "Nb of processes : $p" 
  mpirun -np $p ./src/spmv_mpi_2d_parallelIO ./matrices/$matrix
  echo 
done



