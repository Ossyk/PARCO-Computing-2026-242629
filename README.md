#  OpenMP Sparse Matrix–Vector Multiplication (SpMV) Benchmark

##  Overview
This project implements and benchmarks **Sparse Matrix–Vector Multiplication (SpMV)** using the **Compressed Sparse Row (CSR)** format and **OpenMP** parallelization.  
Specifically, the directive **`#pragma omp parallel for`** for multi-threading and **`#pragma omp simd`** for SIMD vectorization.  

It measures runtime performance across:
- Different thread counts (`1-64 cores`)
- Scheduling strategies (`static`, `dynamic`, `guided`, `auto`) with various chunk sizes (`1`, `10`, `100`, `1000`)
- Optional **SIMD** vectorization (`#pragma omp simd`)

---

##  Prerequisites
Make sure GCC is installed:
```bash
sudo apt-get update
sudo apt install gcc
```

If running on a **cluster**, request a node with at least **64 cores**:
```bash
qsub -I -q short_cpuQ -l select=1:ncpus=64:mem=64gb,walltime=01:00:00
```

---
##  Option 1 — Manual Configuration

###  Import from Git

Clone the repository:
```bash
git clone https://github.com/Ossyk/PARCO-Computing-2026-242629.git
```

Access the project directory:
```bash
cd PARCO-Computing-2026-242629/
```

---

###  Compilation

Compile with GCC inside the project directory:

```bash
gcc src/matrix_generator.c src/csr_utils.c -o src/matrix_generator
gcc -std=c99 -lm src/spMV_seq.c src/csr_utils.c -o src/spMV_seq
gcc -fopenmp -std=c99 -lm src/spMV_parall.c src/csr_utils.c -o src/spMV_parall
gcc -fopenmp -march=native -std=c99 -lm src/simd_vs_parallfor.c src/csr_utils.c -o src/simd_vs_parallfor
```

---

##  Run Instructions

### 1. Create a new CSR matrix
```bash
src/matrix_generator <rows> <cols> <sparsity> <matrix_dest_path>
```
Example: 
 ```bash
src/matrix_generator 2000 2000 0.7 matrices/test1.csr
```

---

###  2. Run the sequential code
```bash
src/spMV_seq <matrix_path>
```
Example:
```bash
src/spMV_seq matrices/m1.csr
```

---

###  3. Run the parallel code
```bash
src/spMV_parall <threads> <scheduler> <chunk> <matrix_path>
```
Example:
```bash
src/spMV_parall 4 static 1 matrices/m4.csr
```

---

###  4. Compare Parallel For vs Parallel For + SIMD
```bash
src/simd_vs_parallfor <matrix_path> <threads> <scheduler>
```
Example:
```bash
src/simd_vs_parallfor matrices/m3.csr 4 dynamic
```

---


###  5. Run full benchmark (custom settings)
```bash
scripts/run_linux_personalized.sh <threads> <scheduler> <chunk> <matrix_path>
```
Example:
```bash
scripts/run_linux_personalized.sh 4 dynamic 10 matrices/m4.csr
```

---

##  Option 2 — Interactive Mode (Recommended)

For an easy experience, clone the project from git, access the project directory and run the interactive menu:
```bash
scripts/interactive_script_l.sh
```

This menu allows you to:
- Create a new matrix  
- View available matrices  
- Run sequential and parallel versions  
- Execute default or personalized benchmarks  


![Interactive Menu](plots/interactive%20menu.jpg)



---

###  Notes
- The code **flushes the CPU cache** between runs to avoid timing bias.  
- Each multiplication is executed **10 times**, and the **90th percentile** runtime is reported for stable measurement.  
- The benchmark automatically detects and displays the number of non-zero elements in the matrix.  

---

###  Author
**Oussema Kasraoui**  
University of Trento – *Parallel Computing (I2PP_D1)*  
Academic Year **2025 / 2026**





