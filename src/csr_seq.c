#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>

typedef struct Sparse_CSR {
	size_t n_rows;
	size_t n_cols;
	size_t n_nz;
	size_t* row_ptrs;
	size_t* col_indices;
	double* values;
} Sparse_CSR;


int create_sparse_csr(
    const double* A,
    size_t n_rows,
    size_t n_cols,
    size_t n_nz,
    Sparse_CSR* A_csr
);

int print_sparse_csr(const Sparse_CSR* A_csr);

int matrix_vector_sparse_csr(
    const Sparse_CSR* A_coo,
    const double* vec,
    double* res
);

void generate_sparse(double *A, int n_rows, int n_cols, double sparsity) {
    for (int i=0; i<n_rows*n_cols; i++) {
        double r = (double)rand()/RAND_MAX;
        A[i] = (r < sparsity) ? 0.0 : (double)(rand()%100)/10.0;
    }
}


double now_sec(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}


int free_sparse_csr(Sparse_CSR* A_csr);


int main (int argc, char** argv) {
   
    int n_rows = (argc > 1) ? atoi(argv[1]) : 1000;
    int n_cols = (argc > 2) ? atoi(argv[2]) : 1000;
    double sparsity = (argc > 3) ? atof(argv[3]) : 0.8;
    int iters = (argc > 4) ? atoi(argv[4]) : 5;

    srand(42);

    double *A = malloc(n_rows * n_cols * sizeof(double));
    double *x = malloc(n_cols * sizeof(double));
    double *Ax = malloc(n_rows * sizeof(double));

    generate_sparse(A, n_rows, n_cols, sparsity);
    for (int i=0; i<n_cols; i++) x[i] = (double)(rand()%100)/10.0;

    size_t n_nz = 0;
    for (int i=0; i<n_rows*n_cols; i++)
        if (A[i] != 0.0) n_nz++;

    Sparse_CSR A_csr;
    create_sparse_csr(A, n_rows, n_cols, n_nz, &A_csr);

    /* Warm-up */
    matrix_vector_sparse_csr(&A_csr, x, Ax);

    /* Benchmark */
    double best = 1e9;
    for (int t=0; t<iters; t++) {
        double t0 = now_sec();
        matrix_vector_sparse_csr(&A_csr, x, Ax);
        double t1 = now_sec();
        double elapsed = (t1 - t0) * 1000.0;
        if (elapsed < best) best = elapsed;
    }

    printf("CSR %dx%d, sparsity=%.2f, nnz=%zu, time=%.3f ms\n",
           n_rows, n_cols, sparsity, A_csr.n_nz, best);
    
    /*
    FILE *f = fopen("results/results_csr.csv", "a");
    if (f) {
        fprintf(f, "CSR,%zu,%zu,%.2f,%zu,%.3f\n",
                A_csr.n_rows, A_csr.n_cols, sparsity, A_csr.n_nz, best);
        fclose(f);
    }
    */

    free_sparse_csr(&A_csr);
    free(A); free(x); free(Ax);
    return EXIT_SUCCESS;


}


int create_sparse_csr(
    const double* A,
    size_t n_rows,
    size_t n_cols,
    size_t n_nz,
    Sparse_CSR* A_csr
) {

	A_csr->n_rows = n_rows;
    A_csr->n_cols = n_cols;
    A_csr->n_nz = n_nz;
    A_csr->row_ptrs = calloc(n_rows+1, sizeof(size_t));
    A_csr->col_indices = calloc(n_nz, sizeof(size_t));
    A_csr->values = calloc(n_nz, sizeof(double));
	
	size_t nz_id=0;
	
	for (size_t i=0; i<n_rows; ++i){
		 A_csr->row_ptrs[i] = nz_id; //under some conditions
		 for (size_t j=0; j<n_cols; ++j){
		 	if (A[i*n_cols + j] != 0.0) {
                A_csr->col_indices[nz_id] = j;
                A_csr->values[nz_id] = A[i*n_cols + j];
                nz_id++;
            }
		 }
	}
	A_csr->row_ptrs[n_rows]=nz_id;
		
    return EXIT_SUCCESS;
}

int print_sparse_csr(const Sparse_CSR* A_csr) {
    
    printf("IP\t");
    for (size_t i=0; i<A_csr->n_rows+1;++i){
        printf("%d\t ",A_csr->row_ptrs[i]);
    }
    printf("\nJ\t");
    for (size_t i=0; i<A_csr->n_nz;++i){
            printf("%d\t ",A_csr->col_indices[i]);
        }
    printf("\nV\t");
    for (size_t i=0; i<A_csr->n_nz;++i){
                printf("%02.2f\t ",A_csr->values[i]);
    }

    return EXIT_SUCCESS;
}

int matrix_vector_sparse_csr(
    const Sparse_CSR* A_csr,
    const double* vec,
    double* res
) {
    size_t nz_id=0;
    for (size_t index=0; index<A_csr->n_rows;index++){
        size_t iterations=A_csr->row_ptrs[index+1]-A_csr->row_ptrs[index];
        res[index]=0.0;
        for (size_t i=0; i<iterations; i++){
            res[index]+=A_csr->values[nz_id]* vec[A_csr->col_indices[nz_id]];
            nz_id++;
        }    
    }
    

    return EXIT_SUCCESS;
}

int free_sparse_csr(Sparse_CSR* A_csr) {
	free(A_csr->row_ptrs);
    free(A_csr->col_indices);
    free(A_csr->values);    
	return EXIT_SUCCESS;
}