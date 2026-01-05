#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define COLOR_GREEN  "\033[32m"
#define COLOR_RESET  "\033[0m"


// ==== Data structures
  

typedef struct {
    int row;
    int col;
    double val;
} Entry;

typedef struct {
    int nrows;
    int nnz;
    int* rowptr;
    int* colind;
    double* values;
} CSR;

// ==== Utilities

int owner(int idx, int size) {
    return idx % size;
}

int compute_local_nrows(int global_nrows, int rank, int size) {
    int count = 0;
    for (int i = rank; i < global_nrows; i += size)
        count++;
    return count;
}



// ==== Local SpMV kernel
  

void local_spmv(const CSR* A, const double* x, double* y) {

    //pragma added for when needed, in default the threads number is 1.
    //and when we will explore the hybrid omp+mpi num_th will be set using export num_threads <T>
    
    #pragma omp parallel for schedule(runtime)
    for (int r = 0; r < A->nrows; r++) {
        double sum = 0.0;
        for (int k = A->rowptr[r]; k < A->rowptr[r + 1]; k++) {
            sum += A->values[k] * x[A->colind[k]];
        }
        y[r] = sum;
    }
}

void flush_cache() {
    const size_t size = 16 * 1024 * 1024; // 16 MB
    static char *buffer = NULL;

    if (!buffer) buffer = (char *)malloc(size);

    for (size_t i = 0; i < size; i++) {
        buffer[i] = (char)(i);
    }
}


// ==== MAIN

int main(int argc, char** argv) {

    int debug = 0;
    
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--debug") == 0) {
            debug = 1;
        }
    }

    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc < 2) {
        if (rank == 0)
            printf("Usage: %s matrix.mtx\n", argv[0]);
        MPI_Finalize();
        return 0;
    }
    
    double t_ghost_start, t_ghost_end;
    double t_comp_start,  t_comp_end;
    double t_gather_start, t_gather_end;
    
    double ghost_time, comp_time, gather_time;
    double local_ghost_time, local_comp_time, local_gather_time;



    int M = 0, N = 0, NNZ = 0;
    Entry* local_entries = NULL;
    int local_nnz = 0;

    // ====  Read & distribute

    if (rank == 0) {
        FILE* f = fopen(argv[1], "r");
        if (!f) MPI_Abort(MPI_COMM_WORLD, 1);

        char line[256];
        do { fgets(line, sizeof(line), f); } while (line[0] == '%');
        sscanf(line, "%d %d %d", &M, &N, &NNZ);

        Entry** buffers = malloc(size * sizeof(Entry*));
        int* counts = calloc(size, sizeof(int));
        for (int p = 0; p < size; p++)
            buffers[p] = malloc(NNZ * sizeof(Entry));

        for (int k = 0; k < NNZ; k++) {
            int i, j; double v;
            fscanf(f, "%d %d %lf", &i, &j, &v);
            i--; j--;
            int p = owner(i, size);
            buffers[p][counts[p]++] = (Entry){i, j, v};
        }
        fclose(f);

        for (int p = 1; p < size; p++) {
            MPI_Send(&counts[p], 1, MPI_INT, p, 0, MPI_COMM_WORLD);
            MPI_Send(buffers[p], counts[p]*sizeof(Entry),
                     MPI_BYTE, p, 1, MPI_COMM_WORLD);
        }

        local_nnz = counts[0];
        local_entries = buffers[0];

        for (int p = 1; p < size; p++) free(buffers[p]);
        free(buffers); free(counts);
    }
    else {
        MPI_Recv(&local_nnz, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        local_entries = malloc(local_nnz * sizeof(Entry));
        MPI_Recv(local_entries, local_nnz*sizeof(Entry),
                 MPI_BYTE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    MPI_Bcast(&M, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // ==== transformation COO → CSR
 

    CSR csr;
    csr.nnz = local_nnz;
    csr.nrows = compute_local_nrows(M, rank, size);

    csr.rowptr = calloc(csr.nrows + 1, sizeof(int));
    csr.colind = malloc(csr.nnz * sizeof(int));
    csr.values = malloc(csr.nnz * sizeof(double));

    for (int k = 0; k < local_nnz; k++) {
        int lr = local_entries[k].row / size;
        csr.rowptr[lr + 1]++;
    }

    for (int i = 0; i < csr.nrows; i++)
        csr.rowptr[i + 1] += csr.rowptr[i];

    int* offset = calloc(csr.nrows, sizeof(int));
    for (int k = 0; k < local_nnz; k++) {
        int lr = local_entries[k].row / size;
        int pos = csr.rowptr[lr] + offset[lr];
        csr.colind[pos] = local_entries[k].col;
        csr.values[pos] = local_entries[k].val;
        offset[lr]++;
    }

    free(offset);
    free(local_entries);
    
    if (debug) {
    printf("Rank %d owns %d rows and %d nonzeros\n",rank, csr.nrows, csr.nnz);
    }
    
// ==== Ghost exchange

    // Build local vector x (owned entries only) 
    int local_vec_n = csr.nrows;
    double* local_x = malloc(local_vec_n * sizeof(double));
    for (int i = 0; i < local_vec_n; i++)
        local_x[i] = 1.0;  // filled with value=1

    // Prepare counts/displs for Allgatherv 
    int* vec_counts = malloc(size * sizeof(int));
    int* vec_displs = malloc(size * sizeof(int));

    for (int p = 0; p < size; p++)
        vec_counts[p] = compute_local_nrows(N, p, size);

    vec_displs[0] = 0;
    for (int p = 1; p < size; p++)
        vec_displs[p] = vec_displs[p - 1] + vec_counts[p - 1];

    double* global_x = malloc(N * sizeof(double));

    MPI_Barrier(MPI_COMM_WORLD);
    t_ghost_start = MPI_Wtime();
    
    MPI_Allgatherv(
        local_x,
        local_vec_n,
        MPI_DOUBLE,
        global_x,
        vec_counts,
        vec_displs,
        MPI_DOUBLE,
        MPI_COMM_WORLD
    );
    
    MPI_Barrier(MPI_COMM_WORLD);
    t_ghost_end = MPI_Wtime();
    
    local_ghost_time = t_ghost_end - t_ghost_start;
    
    MPI_Reduce(&local_ghost_time, &ghost_time,1, MPI_DOUBLE, MPI_MAX, 0,MPI_COMM_WORLD);



    // Now global_x[j] is available → ghost elements resolved
if (debug) {
    printf("Rank %d completed ghost exchange\n", rank);
    }
// ==== Local y_i computation
    
       
    double* local_y = malloc(csr.nrows * sizeof(double));
    
    /* Each rank computes its local y_i */
    
    const int iters = 15;
    double times[iters];
    
    for (int t = 0; t < iters; t++) {
      flush_cache();
    
      MPI_Barrier(MPI_COMM_WORLD);   // sync before timing
      t_comp_start = MPI_Wtime();
      
      local_spmv(&csr, global_x, local_y);
      
      MPI_Barrier(MPI_COMM_WORLD);   // sync after timing
      t_comp_end = MPI_Wtime();
      
      times[t]=t_comp_end - t_comp_start;
    }
    
    // ===Compute 90th percentile 
    for (int i = 0; i < iters - 1; i++) {
        for (int j = i + 1; j < iters; j++) {
            if (times[j] < times[i]) {
                double tmp = times[i];
                times[i] = times[j];
                times[j] = tmp;
            }
        }
    }
    int index = (int)(0.9 * iters);
    if (index >= iters) index = iters - 1;
    local_comp_time = times[index];
    if (debug) {
    printf("Rank %d computed local y_i\n", rank);
}
    MPI_Reduce(&local_comp_time, &comp_time,1, MPI_DOUBLE, MPI_MAX,0, MPI_COMM_WORLD);
   
    
    long long local_flops = 2LL * csr.nnz;
    long long total_flops;
    MPI_Reduce(&local_flops, &total_flops,1, MPI_LONG_LONG, MPI_SUM,0, MPI_COMM_WORLD);
    ///=== print flops and gflops
    if (rank == 0) {
      double gflops = (double)total_flops / comp_time / 1e9;
      printf(COLOR_GREEN "Total FLOPs = %lld  &  " COLOR_RESET, total_flops);
      printf(COLOR_GREEN "Performance = %.3f GFLOP/s\n" COLOR_RESET, gflops);
    }


    //====Balance metrics
    int min_nnz, max_nnz, sum_nnz;
    
    MPI_Reduce(&local_nnz, &min_nnz, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_nnz, &max_nnz, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_nnz, &sum_nnz, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    
    if (rank == 0) {
        printf(COLOR_GREEN "NNZ per rank: min=%d avg=%.1f max=%d\n" COLOR_RESET,min_nnz, (double)sum_nnz/size, max_nnz);
    }


    
    int* y_counts = NULL;
    int* y_displs = NULL;
    double* y_global = NULL;
    
    if (rank == 0) {
        y_counts = malloc(size * sizeof(int));
        y_displs = malloc(size * sizeof(int));
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    t_gather_start = MPI_Wtime();
    
    MPI_Gather(&csr.nrows, 1, MPI_INT,
               y_counts, 1, MPI_INT,
               0, MPI_COMM_WORLD);
    
    if (rank == 0) {
        y_displs[0] = 0;
        for (int p = 1; p < size; p++)
            y_displs[p] = y_displs[p - 1] + y_counts[p - 1];
    
        y_global = malloc((y_displs[size - 1] + y_counts[size - 1]) * sizeof(double));
    }
    
    MPI_Gatherv(local_y, csr.nrows, MPI_DOUBLE,
                y_global, y_counts, y_displs,
                MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    MPI_Barrier(MPI_COMM_WORLD);
    t_gather_end = MPI_Wtime();
    
    local_gather_time = t_gather_end - t_gather_start;

    MPI_Reduce(&local_gather_time, &gather_time,
           1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);


    if (rank == 0) {
    if (debug) {
        printf( "Rank 0 gathered full y vector\n");
    }
    }
    
    //====Printing times results
    
    if (rank == 0) {
    double total_time = ghost_time + comp_time + gather_time;

    printf(COLOR_GREEN "Time per SpMV: %f s (%.1f%%)\n" COLOR_RESET,
           comp_time, 100.0 * comp_time / total_time);

    printf(COLOR_GREEN "Communication time: %f s (%.1f%%)\n" COLOR_RESET,
           ghost_time + gather_time,
           100.0 * (ghost_time + gather_time) / total_time);

    printf(COLOR_GREEN "  ├─ Ghost exchange: %f s (%.1f%%)\n" COLOR_RESET,
           ghost_time, 100.0 * ghost_time / total_time);

    printf(COLOR_GREEN "  └─ Result gather: %f s (%.1f%%)\n" COLOR_RESET,
           gather_time, 100.0 * gather_time / total_time);
    }

// ==== Cleanup

    free(local_x);
    free(global_x);
    free(vec_counts);
    free(vec_displs);

    free(csr.rowptr);
    free(csr.colind);
    free(csr.values);
    free(local_y);

    if (rank == 0) {
        free(y_counts);
        free(y_displs);
        free(y_global);
    }


    MPI_Finalize();
    return 0;
}
