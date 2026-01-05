#include <mpi.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef USE_COLOR
#define COLOR_GREEN "\033[32m"
#define COLOR_RESET "\033[0m"
#else
#define COLOR_GREEN ""
#define COLOR_RESET ""
#endif


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

// ==== 2D Grid utilities

typedef struct {
    MPI_Comm cart_comm;      // Cartesian communicator
    MPI_Comm row_comm;       // Row communicator
    MPI_Comm col_comm;       // Column communicator
    int dims[2];             // Grid dimensions [rows, cols]
    int coords[2];           // Process's coordinates [row, col]
    int proc_row;            // Row index in grid
    int proc_col;            // Column index in grid
} Grid2D;

void setup_2d_grid(Grid2D* grid, int size, int rank) {
    // Create 2D process grid (as square as possible)
    grid->dims[0] = grid->dims[1] = 0;
    MPI_Dims_create(size, 2, grid->dims);
    
    int periods[2] = {0, 0};  
    MPI_Cart_create(MPI_COMM_WORLD, 2, grid->dims, periods, 1, &grid->cart_comm);
    
    // Get process's coordinates
    MPI_Cart_coords(grid->cart_comm, rank, 2, grid->coords);
    grid->proc_row = grid->coords[0];
    grid->proc_col = grid->coords[1];
    
    // Create row and column communicators
    int remain_dims[2];
    
    // Row communicator: same row, varies column
    remain_dims[0] = 0; remain_dims[1] = 1;
    MPI_Cart_sub(grid->cart_comm, remain_dims, &grid->row_comm);
    
    // Column communicator: varies row, same column
    remain_dims[0] = 1; remain_dims[1] = 0;
    MPI_Cart_sub(grid->cart_comm, remain_dims, &grid->col_comm);
}

void free_2d_grid(Grid2D* grid) {
    MPI_Comm_free(&grid->row_comm);
    MPI_Comm_free(&grid->col_comm);
    MPI_Comm_free(&grid->cart_comm);
}

// distribute blocks between processes
int owner_2d_block(int i, int j, const Grid2D* grid, int M, int N) {
    int Pr = grid->dims[0];
    int Pc = grid->dims[1];
    int br = (M + Pr - 1) / Pr;
    int bc = (N + Pc - 1) / Pc;

    int pr = i / br;
    int pc = j / bc;

    if (pr >= Pr) pr = Pr - 1;
    if (pc >= Pc) pc = Pc - 1;

    return pr * Pc + pc;
}

void get_block_bounds(const Grid2D* grid, int M, int N,
                                    int* r0, int* r1, int* c0, int* c1,
                                    int* br_out, int* bc_out) {
    int Pr = grid->dims[0];
    int Pc = grid->dims[1];
    int br = (M + Pr - 1) / Pr;
    int bc = (N + Pc - 1) / Pc;

    int rr0 = grid->proc_row * br;
    int rr1 = rr0 + br; if (rr1 > M) rr1 = M;

    int cc0 = grid->proc_col * bc;
    int cc1 = cc0 + bc; if (cc1 > N) cc1 = N;

    *r0 = rr0; *r1 = rr1;
    *c0 = cc0; *c1 = cc1;
    if (br_out) *br_out = br;
    if (bc_out) *bc_out = bc;
}


// ==== Transform from coo to csr
void build_csr(const Entry* local_entries,
                     int local_nnz,
                     int r0, int r1,
                     int c0, int c1,
                     int rank,
                     CSR* csr_out,
                     int debug)
{
    int local_rows = r1 - r0;

    csr_out->nnz   = local_nnz;
    csr_out->nrows = local_rows;

    csr_out->rowptr = (int*)calloc((size_t)local_rows + 1, sizeof(int));
    csr_out->colind = (int*)malloc((size_t)local_nnz * sizeof(int));
    csr_out->values = (double*)malloc((size_t)local_nnz * sizeof(double));


    
    int valid_nnz = 0;

    for (int k = 0; k < local_nnz; k++) {
        int gi = local_entries[k].row;
        int gj = local_entries[k].col;

        if (gi >= r0 && gi < r1 && gj >= c0 && gj < c1) {
            int lr = gi - r0;
            csr_out->rowptr[lr + 1]++;
            valid_nnz++;
        }
    }

   
    for (int lr = 0; lr < local_rows; lr++) {
        csr_out->rowptr[lr + 1] += csr_out->rowptr[lr];
    }

    
    if (debug) {
        if (csr_out->rowptr[local_rows] != valid_nnz) {
            printf("[Rank %d][WARN] CSR nnz mismatch: rowptr end=%d, valid_nnz=%d\n",
                   rank, csr_out->rowptr[local_rows], valid_nnz);
        }
    }

    // fill CSR
    int* offset = (int*)calloc((size_t)local_rows, sizeof(int));

    for (int k = 0; k < local_nnz; k++) {
        int gi = local_entries[k].row;
        int gj = local_entries[k].col;

        if (gi >= r0 && gi < r1 && gj >= c0 && gj < c1) {
            int lr  = gi - r0;
            int pos = csr_out->rowptr[lr] + offset[lr]++;

            csr_out->colind[pos] = gj;              
            csr_out->values[pos] = local_entries[k].val;
        }
    }

    
    if (debug) {
        for (int lr = 0; lr < local_rows; lr++) {
            int row_nnz = csr_out->rowptr[lr + 1] - csr_out->rowptr[lr];
            if (row_nnz < 0) {
                printf("[Rank %d][ERROR] Negative nnz in local row %d\n",
                       rank, lr);
            }
        }

        printf("[Rank %d] CSR build completed successfully\n", rank);
    }

    free(offset);
}



//====Local SpMV kernel

void local_spmv(const CSR* A, const double* x, double* y) {
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
    const size_t size = 16 * 1024 * 1024;
    static char *buffer = NULL;
    if (!buffer) buffer = (char *)malloc(size);
    for (size_t i = 0; i < size; i++) {
        buffer[i] = (char)(i);
    }
}

// MAIN

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
    
    // Setup 2D grid
    Grid2D grid;
    setup_2d_grid(&grid, size, rank);
    
    if (rank == 0) {
        printf("2D Process grid: %d x %d\n", grid.dims[0], grid.dims[1]);
    }
    
   
    
    double bcast_time = 0.0;
    double reduce_time = 0.0;
    double compute_time = 0.0;
    double local_comp_time = 0.0;



    int M = 0, N = 0, NNZ = 0;
    Entry* local_entries = NULL;
    int local_nnz = 0;
    
    //===Parallel I/O Read
    
    MPI_File fh;
    MPI_File_open(MPI_COMM_WORLD, argv[1], MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    
 
    MPI_Offset header_end = 0;
    if (rank == 0) {
        char line[256];
        FILE* f = fopen(argv[1], "r");
        if (!f) MPI_Abort(MPI_COMM_WORLD, 1);
        
        do {
            header_end = ftell(f);
            fgets(line, sizeof(line), f);
        } while (line[0] == '%');
        
        sscanf(line, "%d %d %d", &M, &N, &NNZ);
        header_end = ftell(f);
        fclose(f);
    }
    
    MPI_Bcast(&M, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&NNZ, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&header_end, 1, MPI_OFFSET, 0, MPI_COMM_WORLD);
    
    // compute the bounds
    int r0, r1, c0, c1, br, bc;
    get_block_bounds(&grid, M, N, &r0, &r1, &c0, &c1, &br, &bc);
    
    int local_rows = r1 - r0;
    int local_cols = c1 - c0;
    
    if (debug) {
        printf("Rank %d grid[%d,%d] block rows [%d,%d) cols [%d,%d) (local_rows=%d local_cols=%d)\n",
               rank, grid.proc_row, grid.proc_col, r0, r1, c0, c1, local_rows, local_cols);
    }
    
    
    
    MPI_Offset file_size;
    MPI_File_get_size(fh, &file_size);
    MPI_Offset data_size = file_size - header_end;
    
    // Each rank reads its part of the file
    MPI_Offset chunk_bytes = data_size / size;
    MPI_Offset my_start = header_end + rank * chunk_bytes;
    MPI_Offset my_end = (rank == size - 1) ? file_size : (header_end + (rank + 1) * chunk_bytes);
    
    int overlap = 256;
    if (rank > 0) my_start -= overlap;
    
    int my_size = (int)(my_end - my_start);
    char* buffer = malloc(my_size + 1);
    
    MPI_File_read_at_all(fh, my_start, buffer, my_size, MPI_CHAR, MPI_STATUS_IGNORE);
    buffer[my_size] = '\0';
    
    MPI_File_close(&fh);
    
    // Find first complete line
    char* parse_start = buffer;
    if (rank > 0) {
        char* first_newline = strchr(buffer, '\n');
        if (first_newline) {
            parse_start = first_newline + 1;
        }
    }
    
    // Find last complete line
    char* parse_end = buffer + my_size;
    if (rank < size - 1) {
        char* last_newline = parse_end;
        while (last_newline > parse_start && *last_newline != '\n') {
            last_newline--;
        }
        if (last_newline > parse_start) {
            parse_end = last_newline + 1;
        }
    }
    
    char saved_char = *parse_end;
    *parse_end = '\0';
    
    // Parse and distribute entries using 2D ownership
    int max_entries_per_rank = (NNZ / size) * 2;
    
    int* counts = calloc(size, sizeof(int));
    Entry** buffers = malloc(size * sizeof(Entry*));
    for (int p = 0; p < size; p++)
        buffers[p] = malloc(max_entries_per_rank * sizeof(Entry));
    
    char* line_ptr = parse_start;
    MPI_Offset first_line_pos = my_start + (parse_start - buffer);
    
    while (line_ptr < parse_end) {
        char* line_end = strchr(line_ptr, '\n');
        if (!line_end) break;
        
        *line_end = '\0';
        
        int i, j; 
        double v;
        if (sscanf(line_ptr, "%d %d %lf", &i, &j, &v) == 3) {
            i--; j--;  // Convert to 0-indexed
            
           
            int p = owner_2d_block(i, j, &grid, M, N);

            
            MPI_Offset current_pos = first_line_pos + (line_ptr - parse_start);
            if (rank == 0 || current_pos >= header_end + rank * chunk_bytes) {
                buffers[p][counts[p]++] = (Entry){i, j, v};
            }
        }
        
        line_ptr = line_end + 1;
    }
    
    *parse_end = saved_char;
    free(buffer);
    
    // All-to-all exchange
    int* send_counts = calloc(size, sizeof(int));
    int* recv_counts = calloc(size, sizeof(int));
    int* send_displs = calloc(size, sizeof(int));
    int* recv_displs = calloc(size, sizeof(int));
    
    memcpy(send_counts, counts, size * sizeof(int));
    
    MPI_Alltoall(send_counts, 1, MPI_INT, recv_counts, 1, MPI_INT, MPI_COMM_WORLD);
    
    for (int p = 1; p < size; p++) {
        send_displs[p] = send_displs[p-1] + send_counts[p-1];
        recv_displs[p] = recv_displs[p-1] + recv_counts[p-1];
    }
    
    int total_send = send_displs[size-1] + send_counts[size-1];
    int total_recv = recv_displs[size-1] + recv_counts[size-1];
    
    Entry* send_buffer = malloc(total_send * sizeof(Entry));
    for (int p = 0; p < size; p++) {
        memcpy(send_buffer + send_displs[p], buffers[p], send_counts[p] * sizeof(Entry));
    }
    
    local_entries = malloc(total_recv * sizeof(Entry));
    
    MPI_Datatype MPI_ENTRY;
    MPI_Type_contiguous(sizeof(Entry), MPI_BYTE, &MPI_ENTRY);
    MPI_Type_commit(&MPI_ENTRY);
    
    MPI_Alltoallv(send_buffer, send_counts, send_displs, MPI_ENTRY,
                  local_entries, recv_counts, recv_displs, MPI_ENTRY,
                  MPI_COMM_WORLD);
    
    MPI_Type_free(&MPI_ENTRY);
    
    local_nnz = total_recv;
    
    long long sum_local = (long long)local_nnz;
    long long sum_total = 0;
    MPI_Allreduce(&sum_local, &sum_total, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) printf("Total received entries across ranks = %lld (expected %d)\n", sum_total, NNZ);
    
    free(send_buffer);
    for (int p = 0; p < size; p++) free(buffers[p]);
    free(buffers); free(counts);
    free(send_counts); free(recv_counts);
    free(send_displs); free(recv_displs);
    
    // transforming COO to CSR
       
    CSR csr;

    build_csr(local_entries,local_nnz,r0, r1,c0, c1,rank,&csr,debug);
    
    free(local_entries);
    local_entries = NULL;

    
    
    if (debug) {
        printf("Rank %d (grid pos [%d,%d]) owns %d rows and %d nonzeros\n",
               rank, grid.proc_row, grid.proc_col, csr.nrows, csr.nnz);
    }
    
  //==== 2D SUMMA SpMV

    double* local_x = malloc(local_cols * sizeof(double));
    for (int i = 0; i < local_cols; i++)
        local_x[i] = 1.0;
    
    double* x_panel = malloc(bc * sizeof(double));
    double* local_y = malloc(local_rows * sizeof(double));
    double* y_final = malloc(local_rows * sizeof(double));
    
    const int iters = 15;
    double times[iters];
    
    double local_bcast_time = 0.0;
    double local_reduce_time = 0.0;
    double local_compute_time = 0.0;
    double comp_time=0.0;
    
    for (int it = 0; it < iters; it++) {
    
        memset(local_y, 0, local_rows * sizeof(double));
        flush_cache();
    
        MPI_Barrier(MPI_COMM_WORLD);
        double t0 = MPI_Wtime();
    
        // SUMMA
        for (int pc = 0; pc < grid.dims[1]; pc++) {

            int root = pc;
        
            if (grid.proc_col == pc) {
                memset(x_panel, 0, bc * sizeof(double));
                memcpy(x_panel, local_x, local_cols * sizeof(double));
            }
        
            double tb0 = MPI_Wtime();
            MPI_Bcast(x_panel, bc, MPI_DOUBLE, root, grid.row_comm);
            double tb1 = MPI_Wtime();
            local_bcast_time += tb1 - tb0;
        
            int j0 = pc * bc;
            int j1 = j0 + bc;
            if (j1 > N) j1 = N;
        
            double tc0 = MPI_Wtime();
            for (int lr = 0; lr < local_rows; lr++) {
                for (int kk = csr.rowptr[lr]; kk < csr.rowptr[lr + 1]; kk++) {
                    int j = csr.colind[kk];
                    if (j >= j0 && j < j1) {
                        local_y[lr] += csr.values[kk] * x_panel[j - j0];
                    }
                }
            }
            double tc1 = MPI_Wtime();
            local_compute_time += tc1 - tc0;
        }

    
        double tr0 = MPI_Wtime();
        MPI_Allreduce(local_y, y_final,
                      local_rows, MPI_DOUBLE, MPI_SUM,
                      grid.row_comm);
        double tr1 = MPI_Wtime();
        local_reduce_time += tr1 - tr0;

    
        MPI_Barrier(MPI_COMM_WORLD);
        double t1 = MPI_Wtime();
    
        times[it] = t1 - t0;
    }
    
    // ---- 90th percentile ----
    for (int i = 0; i < iters - 1; i++)
        for (int j = i + 1; j < iters; j++)
            if (times[j] < times[i]) {
                double tmp = times[i];
                times[i] = times[j];
                times[j] = tmp;
            }
    
    local_comp_time = times[(int)(0.9 * iters)];
    MPI_Reduce(&local_comp_time, &comp_time, 1,
           MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Reduce(&local_bcast_time, &bcast_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_reduce_time, &reduce_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_compute_time, &compute_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    
    // FLOPs
    long long local_flops = 2LL * csr.nnz;
    long long total_flops;
    MPI_Reduce(&local_flops, &total_flops, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    
    if (rank == 0) {
        double gflops = (double)total_flops / comp_time / 1e9;
        printf(COLOR_GREEN "Total FLOPs = %lld & Performance = %.3f GFLOP/s\n" COLOR_RESET,
               total_flops, gflops);
    }


    // Balance metrics
    int min_nnz, max_nnz, sum_nnz;
    
    MPI_Reduce(&local_nnz, &min_nnz, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_nnz, &max_nnz, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_nnz, &sum_nnz, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    
    if (rank == 0) {
        printf(COLOR_GREEN "NNZ per rank: min=%d avg=%.1f max=%d\n" COLOR_RESET,
               min_nnz, (double)sum_nnz/size, max_nnz);
    }

  // Gather y 

    int* y_counts = NULL;
    int* y_displs = NULL;
    double* y_packed = NULL;
    
    /* ---- gather local_rows sizes ---- */
    if (rank == 0) {
        y_counts = malloc(size * sizeof(int));
        y_displs = malloc(size * sizeof(int));
    }
    
    MPI_Gather(&local_rows, 1, MPI_INT,
               y_counts, 1, MPI_INT,
               0, MPI_COMM_WORLD);
    
    /* ---- allocate receive buffer on rank 0 ---- */
    if (rank == 0) {
        y_displs[0] = 0;
        for (int p = 1; p < size; p++)
            y_displs[p] = y_displs[p-1] + y_counts[p-1];
    
        y_packed = malloc((y_displs[size-1] + y_counts[size-1]) * sizeof(double));
    }
    
    /* ---- gather partial y vectors ---- */
    MPI_Gatherv(y_final, local_rows, MPI_DOUBLE,
                y_packed, y_counts, y_displs, MPI_DOUBLE,
                0, MPI_COMM_WORLD);
    
    /* ---- gather row bounds (ALL ranks must call) ---- */
    int my_bounds[2] = { r0, r1 };
    int* all_bounds = NULL;
    
    if (rank == 0) {
        all_bounds = malloc(2 * size * sizeof(int));
    }
    
    MPI_Gather(my_bounds, 2, MPI_INT,
               all_bounds, 2, MPI_INT,
               0, MPI_COMM_WORLD);
    
    /* ---- reorder on rank 0 only ---- */
    if (rank == 0) {
        double* y_global = calloc(M, sizeof(double));
    
        for (int p = 0; p < size; p++) {
            int pr0 = all_bounds[2*p];
            int pr1 = all_bounds[2*p + 1];
            for (int i = 0; i < pr1 - pr0; i++) {
                y_global[pr0 + i] = y_packed[y_displs[p] + i];
            }
        }
    
        free(y_global);
        free(all_bounds);
        free(y_packed);
        free(y_counts);
        free(y_displs);
    }



    // Print timing results
    if (rank == 0) {
        double bcast_avg   = bcast_time   / iters;
        double reduce_avg  = reduce_time  / iters;
        double compute_avg = compute_time / iters;
    
        double total = bcast_avg + reduce_avg + compute_avg;
    
        printf(COLOR_GREEN "2D SUMMA timing breakdown (per iteration avg):\n" COLOR_RESET);
        printf(COLOR_GREEN "  Computation:  %f s (%.1f%%)\n" COLOR_RESET,
               compute_avg, 100.0 * compute_avg / total);
        printf(COLOR_GREEN "  Row Bcast:    %f s (%.1f%%)\n" COLOR_RESET,
               bcast_avg, 100.0 * bcast_avg / total);
        printf(COLOR_GREEN "  Col Reduce:   %f s (%.1f%%)\n" COLOR_RESET,
               reduce_avg, 100.0 * reduce_avg / total);
    }



    //====Cleanup

    free(csr.rowptr);
    free(csr.colind);
    free(csr.values);
    free(local_y);
    free(local_x);
    free(x_panel);
    free(y_final);



    free_2d_grid(&grid);
    MPI_Finalize();
    return 0;
}