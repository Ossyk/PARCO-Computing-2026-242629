# Distributed Sparse Matrix–Vector Multiplication (SpMV)
The project implements SpMV using **MPI** with both **1D cyclic row-wise partitioning**
and **2D block partitioning**, and includes a **single-node OpenMP vs MPI comparison**
as well as an exploratory **hybrid MPI+OpenMP evaluation**.


---

###  Import from Git

Clone the repository:
```bash
git clone https://github.com/Ossyk/PARCO-Computing-2026-242629.git
cd PARCO-Computing-2026-242629/D2-242629/
```

### Requirements

- Linux environment (HPC cluster recommended, requesting a queue with two nodes, providing 128 cores in total)
```bash
qsub -I -q short_cpuQ -l select=2:ncpus=64:mem=64gb,walltime=00:45:00
```
Inside the node:
- GCC (tested with `gcc91`)
- MPI implementation (`mpich-3.2.1--gcc-9.1.0`)
- OpenMP support
```bash
module load gcc91
module load mpich-3.2.1--gcc-9.1.0
module load perf
```

Access the project directory:
```bash
cd PARCO-Computing-2026-242629/D2-242629/
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
PBS job scripts assume that all required matrices already exist.

---
### Compilation 
For detailed debug, add the `--debug` flag.\
To enable color printing in the terminal, add the `-DUSE_COLOR` flag.

#### MPI versions 
```bash
mpicc -fopenmp  src/spmv_mpi.c -o src/spmv_mpi
mpicc -fopenmp src/spmv_mpi_parallelIO.c -o src/spmv_mpi_parallelIO
mpicc -fopenmp src/spmv_mpi_2d_parallelIO.c -o src/spmv_mpi_2d_parallelIO
```


#### OpenMP version
```bash
gcc -fopenmp src/spmv_parallel_openMP.c -o src/spmv_parallel_openMP
```
---

### Running the Code

#### MPI 1D SpMV

Run the 1D cyclic row-wise MPI implementation:

```bash
mpirun -np <P> src/spmv_mpi_parallelIO <matrix>
//Example:
mpirun -np 8 src/spmv_mpi_parallelIO matrices/webbase-1M.mtx
```

#### MPI 2D SpMV
Run the 2D MPI implementation using Cartesian process grids:

```bash
mpirun -np <P> src/spmv_mpi_2d_parallelIO <matrix>
//Example:
mpirun -np 16 src/spmv_mpi_2d_parallelIO matrices/Flan_1565.mtx
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
---
### Provided Experiment Scripts

| Script              | Description                                   | Execution Example |
|---------------------|-----------------------------------------------|-------------------|
| `strong_scaling.sh` | Strong scaling (1–128 MPI ranks)              | `scripts/strong_scaling.sh matrices/webbase-1M.mtx` |
| `weak_scaling.sh`   | Weak scaling with synthetic matrices          | `scripts/weak_scaling.sh` |
| `run_spmv_1d.sh`    | MPI 1D SpMV execution                         | `scripts/run_spmv_1d.sh` |
| `run_spmv_2d.sh`    | MPI 2D SpMV execution                         | `scripts/run_spmv_2d.sh` |
| `run_spmv_hybrid.sh`| Hybrid MPI + OpenMP execution                 | `scripts/run_spmv_hybrid.sh` |

To run the comparison between omp and mpi, request the a node with 64 cores.
```bash
qsub -I -q short_cpuQ -l select=1:ncpus=64:mem=64gb,walltime=00:45:00
```
| Script              | Description                                   | Execution Example |
|---------------------|-----------------------------------------------|-------------------|
| `omp_vs_mpi.sh`     | OpenMP vs MPI comparison on a single node     | `scripts/omp_vs_mpi.sh matrices/webbase-1M.mtx` |

---
### PBS Execution (Cluster)
PBS scripts are provided for reproducibility.
```bash
qsub run_spmv_1d.pbs  //strong and weak scaling using 1D SpMV
qsub run_spmv_2d.pbs   //strong and weak scaling using 2D SpMV
qsub run_spmv_hybrid.pbs   //Hybrid MPI+OpenMP strong scaling
qsub run_spmv_omp_vs_mpi.pbs  //Hybrid MPI vs OpenMP Strong scaling
```
Each PBS script:
- requests 2 nodes (except for run_spmv_omp_vs_mpi.pbs)
- allocates 64 cores per node
- loads required modules
- runs the corresponding experiment

---
### Metrics Reported
For each configuration:
- SpMV execution time (90th percentile of 15 runs)
- Speedup and efficiency
- GFLOP/s (2 × NNZ per SpMV)
- Communication vs computation breakdown
- Load balance (NNZ min/avg/max per rank)
---
###  Author
**Oussema Kasraoui**  
University of Trento – *Intro to Parallel Computing (I2PP_D2)*  
Academic Year **2025 / 2026**









