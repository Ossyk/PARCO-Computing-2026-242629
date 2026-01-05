# Distributed Sparse Matrix–Vector Multiplication (SpMV)
The project implements SpMV using **MPI** with both **1D cyclic row-wise partitioning**
and **2D block partitioning**, and includes a **single-node OpenMP vs MPI comparison**
as well as an exploratory **hybrid MPI+OpenMP evaluation**.


---

###  Import from Git

Clone the repository and access the project directory :
```bash
git clone https://github.com/Ossyk/PARCO-Computing-2026-242629.git
cd PARCO-Computing-2026-242629/D2-242629/
```

### Requirements

- GCC (tested with `gcc91`)
- MPI implementation (`mpich-3.2.1--gcc-9.1.0`)
- OpenMP support
```bash
module load gcc91
module load mpich-3.2.1--gcc-9.1.0
module load perf
```
- Linux environment (HPC cluster recommended, requesting a queue with two nodes, providing 128 cores in total)
```bash
qsub -I -q short_cpuQ -l select=2:ncpus=64:mem=64gb,walltime=00:45:00
```


---

### Missing Large Matrices (Required Setup)

Some large matrices used in the experiments are **not included in the Git repository**
due to size limitations.

Before running any MPI or OpenMP experiments, these matrices must be generated
locally using the provided scripts.

#### Step 1: Compile the matrix generator
```bash
gcc src/random_matrix_generator.c -o src/random_matrix_generator
```
#### Step 2: Generate the missing matrices
```bash
scripts/add_missing_matrices.sh
```

This script generates all required large synthetic matrices in the matrices/
directory.

Note: This setup step must be executed only once.
PBS job scripts assume that all required matrices already exist and do not
perform compilation or matrix generation.

---
### Compilation 
#### MPI versions 
```bash
mpicc -fopenmp  src/spmv_mpi.c -o src/spmv_mp
mpicc src/spmv_mpi_parallelIO.c -o spmv_mpi
mpicc src/spmv_mpi_2d_parallelIO.c -o spmv_mpi_2d
```
For detailed debug, add the `--debug` flag

#### OpenMP version
```bash
gcc -fopenmp src/spmv_parallel_openMP.c -o spmv_omp
```


### Running the Code

#### MPI 1D SpMV

Run the 1D cyclic row-wise MPI implementation:

```bash
mpirun -np <P> ./spmv_mpi_io <matrix>
//Example:
mpirun -np 8 ./spmv_mpi_io matrices/webbase-1M.mtx
```

#### MPI 2D SpMV
Run the 2D MPI implementation using Cartesian process grids:

```bash
mpirun -np <P> ./spmv_mpi_2d <matrix>
//Example:
mpirun -np 16 ./spmv_mpi_2d matrices/Flan_1565.mtx
```

#### OpenMP SpMV (Single Node)
OpenMP execution is intended for single-node runs only.
```bash
qsub -I -q short_cpuQ -l select=1:ncpus=64:mem=64gb,walltime=00:45:00
```

Run the shared-memory OpenMP implementation:

```bash
export OMP_NUM_THREADS=<T>
src/spmv_parallel_openMP <T> <schedule> <chunk> <matrix>
Example:
export OMP_NUM_THREADS=32
src/spmv_parallel_openMP 32 static 1 matrices/webbase-1M.mtx
```
- schedule ∈ {static, dynamic, guided, auto}
- chunk = chunk size (use 0 for default scheduling)

