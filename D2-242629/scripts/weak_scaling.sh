#!/bin/bash



echo "*************************"
echo "1D weak scaling"
echo "*************************"

for i in 1 2 4 8 16 32; do
  echo "Nb of processes: $i"
  name=$((i*20))
  
  mpirun -np $i ./src/spmv_mpi_parallelIO ./matrices/randomN${name}k.mtx
  echo
  
done

echo "Nb of processes: 64"

mpirun -np 64 ./src/spmv_mpi_parallelIO ./matrices/random1M280k.mtx
echo

echo "Nb of processes: 128"

mpirun -np 128 ./src/spmv_mpi_parallelIO ./matrices/random2M560k.mtx
echo 


echo "*************************"
echo "2D Strong scaling"
echo "*************************"

for i in 1 2 4 8 16 32; do
  echo "Nb of processes: $i"
  name=$((i*20))
  

  mpirun -np $i ./src/spmv_mpi_2d_parallelIO ./matrices/randomN${name}k.mtx
  echo
  
done

echo "Nb of processes: 64"

mpirun -np 64 ./src/spmv_mpi_2d_parallelIO ./matrices/random1M280k.mtx
echo

echo "Nb of processes: 128"
mpirun -np 128 ./src/spmv_mpi_2d_parallelIO ./matrices/random2M560k.mtx
echo