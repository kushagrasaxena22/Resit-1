#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

// Function prototype
double get_final_temperatures(int N, int maxIter, double radTemp);

int main() {
    int N = 100; // Grid size
    int maxIter = 10000; // Number of iterations
    double radTemp = 100.0; // Temperature of the radiator

    // Call the function and get the temperature at the center
    double center_temp = get_final_temperatures(N, maxIter, radTemp);

    // Print the result
    printf("Final temperature at the center: %f\n", center_temp);

    return 0;
}

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
