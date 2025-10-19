#!/bin/bash
#SBATCH -J parco_d1
#SBATCH -o results/output_%j.txt
#SBATCH -e results/error_%j.txt
#SBATCH -c 8
#SBATCH --time=00:05:00

module load gcc

echo "==== Dense Sequential ===="
for s in 0.1 0.5 0.8 0.9; do
    ./dense_seq 1000 1000 $s 5 >> results/dense_cluster.csv
done

echo "==== CSR Sequential ===="
for s in 0.1 0.5 0.8 0.9; do
    ./csr_seq 1000 1000 $s 5 >> results/csr_cluster.csv
done
