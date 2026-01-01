// csr_utils.c
#include "csr_utils.h"
#include <stdlib.h>

int create_sparse_csr(const double* A, size_t n_rows, size_t n_cols,
                      size_t n_nz, Sparse_CSR* A_csr) {

    A_csr->n_rows = n_rows;
    A_csr->n_cols = n_cols;
    A_csr->n_nz = n_nz;
    A_csr->row_ptrs   = calloc(n_rows + 1, sizeof(size_t));
    A_csr->col_indices = calloc(n_nz, sizeof(size_t));
    A_csr->values      = calloc(n_nz, sizeof(double));

    size_t nz_id = 0;

    for (size_t i = 0; i < n_rows; ++i) {
        A_csr->row_ptrs[i] = nz_id;
        for (size_t j = 0; j < n_cols; ++j) {
            double val = A[i * n_cols + j];
            if (val != 0.0) {
                A_csr->col_indices[nz_id] = j;
                A_csr->values[nz_id] = val;
                nz_id++;
            }
        }
    }
    A_csr->row_ptrs[n_rows] = nz_id;
    return EXIT_SUCCESS;
}

int free_sparse_csr(Sparse_CSR* A_csr) {
    free(A_csr->row_ptrs);
    free(A_csr->col_indices);
    free(A_csr->values);
    return EXIT_SUCCESS;
}
