// csr_utils.h
#ifndef CSR_UTILS_H
#define CSR_UTILS_H

#include <stdlib.h>

typedef struct Sparse_CSR {
    size_t n_rows;
    size_t n_cols;
    size_t n_nz;
    size_t* row_ptrs;
    size_t* col_indices;
    double* values;
} Sparse_CSR;

int create_sparse_csr(const double* A, size_t n_rows, size_t n_cols,
                      size_t n_nz, Sparse_CSR* A_csr);

int free_sparse_csr(Sparse_CSR* A_csr);

#endif
