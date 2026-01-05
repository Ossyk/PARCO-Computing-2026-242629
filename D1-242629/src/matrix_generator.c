// matrix_generator.c
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "csr_utils.h"

void generate_sparse(double *A, int n_rows, int n_cols, double sparsity) {
    for (int i = 0; i < n_rows * n_cols; i++) {
        double r = (double)rand() / RAND_MAX;
        A[i] = (r < sparsity) ? 0.0 : (double)(rand() % 100) / 10.0;
    }
}

int main(int argc, char **argv) {
    if (argc < 5) {
        fprintf(stderr, "Usage: %s <rows> <cols> <sparsity> <output_file>\n", argv[0]);
        return 1;
    }

    int n_rows = atoi(argv[1]);
    int n_cols = atoi(argv[2]);
    double sparsity = atof(argv[3]);
    const char *filename = argv[4];

    srand(42);

    double *A = malloc(n_rows * n_cols * sizeof(double));
    generate_sparse(A, n_rows, n_cols, sparsity);

    // Count non-zeros
    size_t n_nz = 0;
    for (int i = 0; i < n_rows * n_cols; i++)
        if (A[i] != 0.0) n_nz++;

    Sparse_CSR A_csr;
    create_sparse_csr(A, n_rows, n_cols, n_nz, &A_csr);

    FILE *f = fopen(filename, "w");
    if (!f) {
        perror("fopen");
        return 1;
    }

    fprintf(f, "%zu %zu %zu\n", A_csr.n_rows, A_csr.n_cols, A_csr.n_nz);
    for (size_t i = 0; i <= A_csr.n_rows; i++)
        fprintf(f, "%zu ", A_csr.row_ptrs[i]);
    fprintf(f, "\n");
    for (size_t i = 0; i < A_csr.n_nz; i++)
        fprintf(f, "%zu ", A_csr.col_indices[i]);
    fprintf(f, "\n");
    for (size_t i = 0; i < A_csr.n_nz; i++)
        fprintf(f, "%lf ", A_csr.values[i]);
    fprintf(f, "\n");

    fclose(f);
    free_sparse_csr(&A_csr);
    free(A);

    printf("Generated CSR matrix: %s (%d x %d, sparsity %.2f, nnz=%zu)\n",
           filename, n_rows, n_cols, sparsity, n_nz);

    return 0;
}
