#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/*
  Generate a random sparse matrix in Matrix Market format
  Each row has exactly K nonzeros (including diagonal)
*/

int main(int argc, char** argv) {
    if (argc < 5) {
        fprintf(stderr,
            "Usage: %s <N> <K> <seed> <output.mtx>\n"
            "  N    : matrix size (NxN)\n"
            "  K    : nonzeros per row (K >= 1)\n"
            "  seed : random seed\n",
            argv[0]);
        return 1;
    }

    int N = atoi(argv[1]);
    int K = atoi(argv[2]);
    int seed = atoi(argv[3]);
    const char* out = argv[4];

    if (K < 1 || K > N) {
        fprintf(stderr, "Error: K must be in [1, N]\n");
        return 1;
    }

    srand(seed);

    long long NNZ = (long long)N * K;

    FILE* f = fopen(out, "w");
    if (!f) {
        perror("fopen");
        return 1;
    }

    /* Matrix Market header */
    fprintf(f, "%%%%MatrixMarket matrix coordinate real general\n");
    fprintf(f, "%% Random synthetic matrix: N=%d, K=%d\n", N, K);
    fprintf(f, "%d %d %lld\n", N, N, NNZ);

    for (int i = 0; i < N; i++) {
        /* Always include diagonal */
        fprintf(f, "%d %d %f\n", i + 1, i + 1, 1.0);

        /* Select remaining K-1 random columns */
        int written = 1;
        while (written < K) {
            int j = rand() % N;
            if (j == i) continue;  // avoid duplicate diagonal
            fprintf(f, "%d %d %f\n", i + 1, j + 1, 1.0);
            written++;
        }
    }

    fclose(f);
    return 0;
}
