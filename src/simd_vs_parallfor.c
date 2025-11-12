// spmv_full_benchmark_cli.c
#define _POSIX_C_SOURCE 200112L   // enable POSIX setenv()
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <time.h>

#include "csr_utils.h"

double now_sec(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

/* --- Variant 1: Parallel For --- */
void spmv_parallel_for(const Sparse_CSR* A, const double* x, double* y) {
    #pragma omp parallel for schedule(runtime)
    for (size_t r = 0; r < A->n_rows; r++) {
        double sum = 0.0;
        for (size_t i = A->row_ptrs[r]; i < A->row_ptrs[r+1]; i++)
            sum += A->values[i] * x[A->col_indices[i]];
        y[r] = sum;
    }
}

/* --- Variant 2: Parallel For + SIMD --- */
void spmv_parallel_simd(const Sparse_CSR* A, const double* x, double* y) {
    #pragma omp parallel for schedule(runtime)
    for (size_t i = 0; i < A->n_rows; i++) {
        double sum = 0.0;

        
        #pragma omp simd reduction(+:sum)
        for (size_t j = A->row_ptrs[i]; j < A->row_ptrs[i+1]; j++)
            sum += A->values[j] * x[A->col_indices[j]];
        y[i] = sum;


    }
}

/* --- Benchmark helper --- */
typedef void (*spmv_func)(const Sparse_CSR*, const double*, double*);

double benchmark(spmv_func func, const Sparse_CSR* A, const double* x, double* y) {
    const int iters = 10;
    double times[iters];

    func(A, x, y); // warm-up

    for (int t = 0; t < iters; t++) {
        double t0 = now_sec();
        func(A, x, y);
        double t1 = now_sec();
        times[t] = (t1 - t0) * 1000.0;
    }

    // Sort for 90th percentile
    for (int i = 0; i < iters - 1; i++)
        for (int j = i + 1; j < iters; j++)
            if (times[j] < times[i]) {
                double tmp = times[i];
                times[i] = times[j];
                times[j] = tmp;
            }

    int idx = (int)(0.9 * iters);
    if (idx >= iters) idx = iters - 1;
    return times[idx];
}

/* --- Utility to set environment variables --- */
void set_env(const char *var, const char *val) {
#ifdef _WIN32
    _putenv_s(var, val);
#else
    setenv(var, val, 1);
#endif
}

/* --- Main --- */
int main(int argc, char **argv) {
    if (argc < 4) {
        fprintf(stderr, "Usage: %s <matrix.csr> <threads> <schedule>\n", argv[0]);
        fprintf(stderr, "Example: %s matrices/m2.csr 4 dynamic,8\n", argv[0]);
        return 1;
    }

    const char *fname = argv[1];
    const char *threads = argv[2];
    const char *schedule = argv[3];

    set_env("OMP_NUM_THREADS", threads);
    set_env("OMP_SCHEDULE", schedule);

    /* --- Read CSR matrix --- */
    FILE *f = fopen(fname, "r");
    if (!f) { perror("fopen"); return 1; }

    size_t rows, cols, nnz;
    fscanf(f, "%zu %zu %zu", &rows, &cols, &nnz);

    Sparse_CSR A;
    A.n_rows = rows; A.n_cols = cols; A.n_nz = nnz;
    A.row_ptrs = malloc((rows + 1) * sizeof(size_t));
    A.col_indices = malloc(nnz * sizeof(size_t));
    A.values = malloc(nnz * sizeof(double));

    for (size_t i = 0; i <= rows; i++) fscanf(f, "%zu", &A.row_ptrs[i]);
    for (size_t i = 0; i < nnz; i++) fscanf(f, "%zu", &A.col_indices[i]);
    for (size_t i = 0; i < nnz; i++) fscanf(f, "%lf", &A.values[i]);
    fclose(f);

    double *x = malloc(cols * sizeof(double));
    double *y = malloc(rows * sizeof(double));
    for (size_t j = 0; j < cols; j++) x[j] = 1.0;

    printf("===================================\n");
    printf("Running SpMV Benchmark (CLI)\n");
    printf("Matrix: %s\n", fname);
    printf("Threads: %s | Schedule: %s\n", threads, schedule);
    printf("Rows=%zu, Cols=%zu, NNZ=%zu\n", rows, cols, nnz);
    printf("===================================\n\n");

    /* --- Benchmark parallel for --- */
    double t_for = benchmark(spmv_parallel_for, &A, x, y);
    /* --- Benchmark simd variant --- */
    double t_simd = benchmark(spmv_parallel_simd, &A, x, y);

    printf("[parallel for]        90th percentile = %.3f ms\n", t_for);
    printf("[parallel for + simd] 90th percentile = %.3f ms\n", t_simd);
    printf("Speedup (simd/for) = %.2fx\n", t_for / t_simd);

    printf("\n===================================\n");
    printf("Benchmark Complete\n");
    printf("===================================\n");

    free_sparse_csr(&A);
    free(x); free(y);
    return 0;
}
