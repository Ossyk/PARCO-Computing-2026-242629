// spmv_seq.c
#define _POSIX_C_SOURCE 199309L

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "csr_utils.h"

double now_sec(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

int spMV_csr_multiplication(const Sparse_CSR* A_csr, const double* vec, double* res) {
    
    for (size_t row = 0; row < A_csr->n_rows; row++) {
        double sum = 0.0;
        for (size_t i = A_csr->row_ptrs[row]; i < A_csr->row_ptrs[row + 1]; i++)
            sum += A_csr->values[i] * vec[A_csr->col_indices[i]];
        res[row] = sum;
    }


    return EXIT_SUCCESS;
}

void flush_cache() {
    const size_t size = 8 * 1024 * 1024; // 8 MB
    static char *buffer = NULL;

    if (!buffer) buffer = (char *)malloc(size);

    for (size_t i = 0; i < size; i++) {
        buffer[i] = (char)(i);
    }
}


int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <matrix.csr>\n", argv[0]);
        return 1;
    }

    const char *fname = argv[1];
    FILE *f = fopen(fname, "r");
    if (!f) { perror("fopen"); return 1; }

    size_t rows, cols, nnz;
    fscanf(f, "%zu %zu %zu", &rows, &cols, &nnz);

    Sparse_CSR A_csr;
    A_csr.n_rows = rows; A_csr.n_cols = cols; A_csr.n_nz = nnz;
    A_csr.row_ptrs = malloc((rows + 1) * sizeof(size_t));
    A_csr.col_indices = malloc(nnz * sizeof(size_t));
    A_csr.values = malloc(nnz * sizeof(double));

    for (size_t i = 0; i <= rows; i++) fscanf(f, "%zu", &A_csr.row_ptrs[i]);
    for (size_t i = 0; i < nnz; i++) fscanf(f, "%zu", &A_csr.col_indices[i]);
    for (size_t i = 0; i < nnz; i++) fscanf(f, "%lf", &A_csr.values[i]);
    fclose(f);

    double *x = malloc(cols * sizeof(double));
    double *y = malloc(rows * sizeof(double));
    for (size_t j = 0; j < cols; j++) x[j] = 1.0;

    /* Warm-up */
    spMV_csr_multiplication(&A_csr, x, y);

    /* Benchmark */
    int iters = 10;
    double times[iters];
    for (int t = 0; t < iters; t++) {
        
        flush_cache();

        double t0 = now_sec();
        spMV_csr_multiplication(&A_csr, x, y);
        double t1 = now_sec();
        times[t] = (t1 - t0) * 1000.0;
    }

    // Sort and find 90th percentile
    for (int i = 0; i < iters - 1; i++)
        for (int j = i + 1; j < iters; j++)
            if (times[j] < times[i]) {
                double tmp = times[i];
                times[i] = times[j];
                times[j] = tmp;
            }

    int index = (int)(0.9 * iters);
    if (index >= iters) index = iters - 1;
    double p90 = times[index];

    printf("Sequential CSR %zux%zu, nnz=%zu, 90th percentile: %.3f ms\n",
           rows, cols, nnz, p90);

    free_sparse_csr(&A_csr);
    free(x); free(y);
    return 0;
}
