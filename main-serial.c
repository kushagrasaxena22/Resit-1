
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

// Function prototypes for file reading/writing
int read_num_of_temps(const char* filename);
double* read_temps(const char* filename, int numOfTemps);
void write_to_output_file(const char* filename, double* results, int numOfTemps);

// Function to calculate the final temperature at the center of the room
double get_final_temperatures(int N, int maxIter, double radTemp) {
    // Allocate memory for the temperature grid
    double** curr_t = (double**) malloc(N * sizeof(double*));
    double** prev_t = (double**) malloc(N * sizeof(double*));
    for (int i = 0; i < N; i++) {
        curr_t[i] = (double*) malloc(N * sizeof(double));
        prev_t[i] = (double*) malloc(N * sizeof(double));
    }

    // Initialize the grid
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (i == N - 1 && j >= floor((N-1) * 0.3) && j <= ceil((N-1) * 0.7)) {
                curr_t[i][j] = radTemp; // Radiator
            } else {
                curr_t[i][j] = 10.0; // Initial room temperature
            }
        }
    }

    // Iterative process using OpenMP
    for (int iter = 0; iter < maxIter; iter++) {
        // Swap current and previous grids
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

    // Get the temperature at the center of the room
    int pointX = floor((N-1) * 0.5);
    int pointY = floor((N-1) * 0.5);
    double result = curr_t[pointX][pointY];

    // Free allocated memory
    for (int i = 0; i < N; i++) {
        free(curr_t[i]);
        free(prev_t[i]);
    }
    free(curr_t);
    free(prev_t);

    return result;
}

// Function to read the number of temperatures from the input file
int read_num_of_temps(const char* filename) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        perror("Failed to open file");
        exit(EXIT_FAILURE);
    }
    int numOfTemps;
    fscanf(file, "%d", &numOfTemps);
    fclose(file);
    return numOfTemps;
}

// Function to read temperatures from the input file
double* read_temps(const char* filename, int numOfTemps) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        perror("Failed to open file");
        exit(EXIT_FAILURE);
    }
    fscanf(file, "%*d"); // Skip the first integer
    double* temps = (double*) malloc(numOfTemps * sizeof(double));
    for (int i = 0; i < numOfTemps; i++) {
        fscanf(file, "%lf", &temps[i]);
    }
    fclose(file);
    return temps;
}

// Function to write results to the output file
void write_to_output_file(const char* filename, double* results, int numOfTemps) {
    FILE* file = fopen(filename, "w");
    if (!file) {
        perror("Failed to open file");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < numOfTemps; i++) {
        fprintf(file, "%lf\n", results[i]);
    }
    fclose(file);
}

int main(int argc, char* argv[]) {
    if (argc != 5) {
        printf("Usage: %s <N> <maxIter> <input_file> <output_file>\n", argv[0]);
        return 1;
    }

    int N = atoi(argv[1]);
    int maxIter = atoi(argv[2]);
    char* input_file = argv[3];
    char* output_file = argv[4];

    // Read number of temperatures and the temperature values
    int numOfTemps = read_num_of_temps(input_file);
    double* temps = read_temps(input_file, numOfTemps);

    // Allocate memory for storing results
    double* results = (double*) malloc(numOfTemps * sizeof(double));

    // Compute final temperature for each radiator temperature
    for (int i = 0; i < numOfTemps; i++) {
        results[i] = get_final_temperatures(N, maxIter, temps[i]);
    }

    // Write results to the output file
    write_to_output_file(output_file, results, numOfTemps);

    // Free allocated memory
    free(temps);
    free(results);

    return 0;
}
