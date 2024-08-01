#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

// Placeholder for file-reader.h functions
int read_num_of_temps(const char* filename) {
    FILE* file = fopen(filename, "r");
    int num;
    fscanf(file, "%d", &num);
    fclose(file);
    return num;
}

double* read_temps(const char* filename, int numOfTemps) {
    FILE* file = fopen(filename, "r");
    int num;
    fscanf(file, "%d", &num); // Skip the first line
    double* temps = (double*) malloc(numOfTemps * sizeof(double));
    for (int i = 0; i < numOfTemps; i++) {
        fscanf(file, "%lf", &temps[i]);
    }
    fclose(file);
    return temps;
}

void write_to_output_file(const char* filename, double* results, int numOfTemps) {
    FILE* file = fopen(filename, "w");
    fprintf(file, "%d\n", numOfTemps);
    for (int i = 0; i < numOfTemps; i++) {
        fprintf(file, "%f ", results[i]);
    }
    fclose(file);
}

// Function to calculate final temperatures using OpenMP (from previous code)
double get_final_temperatures(int N, int maxIter, double radTemp) {
    double** curr_t = (double**) malloc(N * sizeof(double*));
    double** prev_t = (double**) malloc(N * sizeof(double*));
    for (int i = 0; i < N; i++) {
        curr_t[i] = (double*) malloc(N * sizeof(double));
        prev_t[i] = (double*) malloc(N * sizeof(double));
    }

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (i == N - 1 && j >= floor((N-1) * 0.3) && j <= ceil((N-1) * 0.7)) {
                curr_t[i][j] = radTemp; 
            } else {
                curr_t[i][j] = 10.0;
            }
        }
    }

    for (int iter = 0; iter < maxIter; iter++) {
        double** temp = prev_t;
        prev_t = curr_t;
        curr_t = temp;

        #pragma omp parallel for collapse(2)
        for (int i = 1; i < N-1; i++) {
            for (int j = 1; j < N-1; j++) {
                curr_t[i][j] = (prev_t[i+1][j] + prev_t[i-1][j] +
                                prev_t[i][j+1] + prev_t[i][j-1]) / 4.0;
            }
        }
    }

    int pointX = floor((N-1) * 0.5);
    int pointY = floor((N-1) * 0.5);
    double result = curr_t[pointX][pointY];

    for (int i = 0; i < N; i++) {
        free(curr_t[i]);
        free(prev_t[i]);
    }
    free(curr_t);
    free(prev_t);

    return result;
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc != 5) {
        if (rank == 0) {
            printf("Usage: %s <N> <maxIter> <input_file> <output_file>\n", argv[0]);
        }
        MPI_Finalize();
        return 1;
    }

    int N = atoi(argv[1]);
    int maxIter = atoi(argv[2]);
    char* input_file = argv[3];
    char* output_file = argv[4];

    int numOfTemps, local_numOfTemps;
    double *temps = NULL, *local_temps = NULL, *local_results = NULL;

    if (rank == 0) {
        numOfTemps = read_num_of_temps(input_file);
        temps = read_temps(input_file, numOfTemps);
    }

    MPI_Bcast(&numOfTemps, 1, MPI_INT, 0, MPI_COMM_WORLD);

    local_numOfTemps = numOfTemps / size;
    int remainder = numOfTemps % size;
    if (rank < remainder) {
        local_numOfTemps++;
    }

    local_temps = (double*) malloc(local_numOfTemps * sizeof(double));
    local_results = (double*) malloc(local_numOfTemps * sizeof(double));

    int* sendcounts = (int*) malloc(size * sizeof(int));
    int* displs = (int*) malloc(size * sizeof(int));
    int offset = 0;
    for (int i = 0; i < size; i++) {
        sendcounts[i] = numOfTemps / size + (i < remainder ? 1 : 0);
        displs[i] = offset;
        offset += sendcounts[i];
    }

    MPI_Scatterv(temps, sendcounts, displs, MPI_DOUBLE, local_temps, local_numOfTemps, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    for (int i = 0; i < local_numOfTemps; i++) {
        local_results[i] = get_final_temperatures(N, maxIter, local_temps[i]);
    }

    double* results = NULL;
    if (rank == 0) {
        results = (double*) malloc(numOfTemps * sizeof(double));
    }
    MPI_Gatherv(local_results, local_numOfTemps, MPI_DOUBLE, results, sendcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        write_to_output_file(output_file, results, numOfTemps);
        free(results);
        free(temps);
    }

    free(local_temps);
    free(local_results);
    free(sendcounts);
    free(displs);

    MPI_Finalize();
    return 0;
}
