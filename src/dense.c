#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/* ---------- Dense Matrix-Vector Multiplication ---------- */
void dense_matvec(int n_rows, int n_cols, double *M, double *V, double *res) {
    for (int i = 0; i < n_rows; i++) {
        double sum = 0.0;
        for (int j = 0; j < n_cols; j++) {
            sum += M[i * n_cols + j] * V[j];
        }
        res[i] = sum;
    }
}

/* ---------- Sparse Matrix Generator (by sparsity %) ---------- */
void generate_sparse_matrix(double *M, int n_rows, int n_cols, double sparsity) {
    for (int i = 0; i < n_rows * n_cols; i++) {
        double r = (double) rand() / RAND_MAX;
        if (r < sparsity)
            M[i] = 0.0;                   // make this element zero
        else
            M[i] = (double)(rand() % 100) / 10.0; // some random non-zero
    }
}

/* ---------- Current Time (seconds) ---------- */
double now_sec(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

/* ---------- Main ---------- */
int main(int argc, char **argv) {

    /* Default parameters or from command line */
    int n_rows = (argc > 1) ? atoi(argv[1]) : 1000;
    int n_cols = (argc > 2) ? atoi(argv[2]) : 1000;
    double sparsity = (argc > 3) ? atof(argv[3]) : 0.8;   // 0.8 → 80% zeros
    int iters = (argc > 4) ? atoi(argv[4]) : 5;           // repeat for averaging

    if (argc<2){
        printf("Usage [row_size] [col_size] [sparsity] [iters]\n");
        printf("Default settings: 1000 1000 0.8 5\n");
    }

    srand(42); // fixed seed for reproducibility

    /* Allocate */
    double *M = malloc(n_rows * n_cols * sizeof(double));
    double *V = malloc(n_cols * sizeof(double));
    double *res = malloc(n_rows * sizeof(double));

    if (!M || !V || !res) {
        fprintf(stderr, "Memory allocation failed\n");
        return EXIT_FAILURE;
    }

    /* Generate random sparse matrix and vector */
    generate_sparse_matrix(M, n_rows, n_cols, sparsity);
    for (int i = 0; i < n_cols; i++)
        V[i] = (double)(rand() % 100) / 10.0;

    /* Warm-up run */
    dense_matvec(n_rows, n_cols, M, V, res);

    /* Benchmark */
    double best = 1e9;
    for (int t = 0; t < iters; t++) {
        double start = now_sec();
        dense_matvec(n_rows, n_cols, M, V, res);
        double end = now_sec();
        double elapsed = (end - start) * 1000.0; // ms
        if (elapsed < best) best = elapsed;
    }

    printf("Matrix %dx%d, sparsity=%.2f, time=%.3f ms\n",
           n_rows, n_cols, sparsity, best);
    
           /*
           
    FILE *f = fopen("results/results_dense.csv", "a");
    if (f) {
        fprintf(f, "Matrix,%zu,%zu,%.2f,%.3f\n",
                n_rows, n_cols, sparsity, best);
        fclose(f);
    }*/


    /* Free */
    free(M);
    free(V);
    free(res);
    return 0;
}
