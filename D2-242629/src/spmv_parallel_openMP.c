
#define _POSIX_C_SOURCE 199309L

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <omp.h>
  

typedef struct {
    size_t row;
    size_t col;
    double val;
} COOEntry;

typedef struct {
    size_t n_rows;
    size_t n_cols;
    size_t n_nz;
    size_t *row_ptrs;
    size_t *col_indices;
    double *values;
} Sparse_CSR;

double now_sec(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

/* Parallel CSR SpMV */
int spMV_csr_parallel(const Sparse_CSR* A_csr, const double* vec, double* res) {


    #pragma omp parallel for schedule(runtime)
    for (size_t row = 0; row < A_csr->n_rows; row++) {
        double sum = 0.0;
        for (size_t i = A_csr->row_ptrs[row]; i < A_csr->row_ptrs[row + 1]; i++)
            sum += A_csr->values[i] * vec[A_csr->col_indices[i]];
        res[row] = sum;
    }

    
    return EXIT_SUCCESS;
}

void flush_cache() {
    const size_t size = 16 * 1024 * 1024; // 8 MB
    static char *buffer = NULL;

    if (!buffer) buffer = (char *)malloc(size);

    for (size_t i = 0; i < size; i++) {
        buffer[i] = (char)(i);
    }
}

void coo_to_csr(size_t nrows, size_t ncols, size_t nnz,
                COOEntry *coo, Sparse_CSR *csr)
{
    csr->n_rows = nrows;
    csr->n_cols = ncols;
    csr->n_nz   = nnz;

    csr->row_ptrs    = calloc(nrows + 1, sizeof(size_t));
    csr->col_indices = malloc(nnz * sizeof(size_t));
    csr->values      = malloc(nnz * sizeof(double));

    /* Count nnz per row */
    for (size_t i = 0; i < nnz; i++)
        csr->row_ptrs[coo[i].row + 1]++;

    /* Prefix sum */
    for (size_t i = 0; i < nrows; i++)
        csr->row_ptrs[i + 1] += csr->row_ptrs[i];

    size_t *cursor = malloc(nrows * sizeof(size_t));
    memcpy(cursor, csr->row_ptrs, nrows * sizeof(size_t));

    for (size_t i = 0; i < nnz; i++) {
        size_t r = coo[i].row;
        size_t pos = cursor[r]++;
        csr->col_indices[pos] = coo[i].col;
        csr->values[pos]      = coo[i].val;
    }

    free(cursor);
}

void free_sparse_csr(Sparse_CSR *A) {
    free(A->row_ptrs);
    free(A->col_indices);
    free(A->values);
}


int main(int argc, char **argv) {
    if (argc < 5) {
        fprintf(stderr, "Usage: %s <num_threads> <schedule> <chunk_size> <matrix.csr>\n", argv[0]);
        fprintf(stderr, "  schedule: static, dynamic, guided, auto\n");
        fprintf(stderr, "  chunk_size: positive integer (0 for default)\n");
        fprintf(stderr, "Example: %s 4 static 100 matrices/m1.csr\n", argv[0]);
        fprintf(stderr, "Example: %s 4 dynamic 0 matrices/m1.csr  (default chunk)\n", argv[0]);
        return 1;
    }

    int num_threads = atoi(argv[1]);
    const char *schedule = argv[2];
    int chunk_size = atoi(argv[3]);
    const char *fname = argv[4];

    if (num_threads < 1) {
        fprintf(stderr, "Error: num_threads must be >= 1\n");
        return 1;
    }
    omp_set_num_threads(num_threads);

    if (chunk_size < 0) {
        fprintf(stderr, "Error: chunk_size must be >= 0\n");
        return 1;
    }

    // Validate and set schedule with chunk size
    if (strcmp(schedule, "static") == 0) {
        omp_set_schedule(omp_sched_static, chunk_size);
    } else if (strcmp(schedule, "dynamic") == 0) {
        omp_set_schedule(omp_sched_dynamic, chunk_size);
    } else if (strcmp(schedule, "guided") == 0) {
        omp_set_schedule(omp_sched_guided, chunk_size);
    } else if (strcmp(schedule, "auto") == 0) {
        omp_set_schedule(omp_sched_auto, chunk_size);
    } else {
        fprintf(stderr, "Error: Invalid schedule '%s'. Use: static, dynamic, guided, or auto\n", schedule);
        return 1;
    }

    /* ---------- Read mtx ---------- */
    FILE *f = fopen(fname, "r");
    if (!f) { perror("fopen"); return 1; }

    char line[512];
    do { fgets(line, sizeof(line), f); }
    while (line[0] == '%');

    size_t rows, cols, nnz;
    sscanf(line, "%zu %zu %zu", &rows, &cols, &nnz);

    COOEntry *coo = malloc(nnz * sizeof(COOEntry));

    for (size_t i = 0; i < nnz; i++) {
        size_t r, c;
        double v;
        fscanf(f, "%zu %zu %lf", &r, &c, &v);
        coo[i].row = r - 1;
        coo[i].col = c - 1;
        coo[i].val = v;
    }
    fclose(f);

    Sparse_CSR A;
    coo_to_csr(rows, cols, nnz, coo, &A);
    free(coo);

    double *x = malloc(cols * sizeof(double));
    double *y = malloc(rows * sizeof(double));

    for (size_t i = 0; i < cols; i++)
        x[i] = (double)rand() / RAND_MAX;

    /* Warm-up */
    spMV_csr_parallel(&A, x, y);

    /* Benchmark */
    const int iters = 15;
    double times[iters];

    for (int i = 0; i < iters; i++) {
        flush_cache();
        double t0 = now_sec();
        spMV_csr_parallel(&A, x, y);
        double t1 = now_sec();
        times[i] = (t1 - t0) * 1000.0;
    }

    /* p90 */
    for (int i = 0; i < iters - 1; i++)
        for (int j = i + 1; j < iters; j++)
            if (times[j] < times[i]) {
                double tmp = times[i];
                times[i] = times[j];
                times[j] = tmp;
            }

    double p90 = times[(int)(0.9 * iters)];
    double sec = p90 / 1000.0;
    double gflops = (2.0 * nnz) / (sec * 1e9);

    printf("OpenMP CSR %zux%zu nnz=%zu | threads=%d | p90=%.3f ms | GFLOP/s=%.3f\n",
           rows, cols, nnz, num_threads, p90, gflops);

    free_sparse_csr(&A);
    free(x);
    free(y);
    return 0;
}