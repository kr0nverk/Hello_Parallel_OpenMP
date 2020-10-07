#include <iostream>
#include <iomanip>
#include <omp.h>
#include <ctime>

using namespace std;

double MaxMinInMatrix(unsigned num_size, double** Matrix) {
    double max_min_in_matrix = .0;
    
    double* Min = new double[num_size];
    for (int i = 0; i < num_size; ++i) {
        Min[i] = Matrix[i][0];
        for (int j = 0; j < num_size; ++j) {
            if (Matrix[i][j] < Min[i]) {
                Min[i] = Matrix[i][j];
            }
        }

        if (Min[i] > max_min_in_matrix) {
            max_min_in_matrix = Min[i];
        }
    }
    delete[]Min;

    return max_min_in_matrix;
}

double MaxMinInMatrixParallel(unsigned num_size, double** Matrix) {
    double max_min_in_matrix_parallel = .0;

    double* Min = new double[num_size];
    int i, j;

    #pragma omp parallel
    {
        #pragma omp for private(i,j)
        for (i = 0; i < num_size; ++i) {
            Min[i] = Matrix[i][0];
            for (j = 0; j < num_size; ++j) {
                if (Matrix[i][j] < Min[i]) {
                    Min[i] = Matrix[i][j];
                }
            }

            #pragma omp flush(max_min_in_matrix_parallel)
            if (Min[i] > max_min_in_matrix_parallel) {

                #pragma omp critical
                if (Min[i] > max_min_in_matrix_parallel) {
                    max_min_in_matrix_parallel = Min[i];
                }
            }
        }
    }
    delete[]Min;

    return max_min_in_matrix_parallel;
}

double MaxMinInTriangularMatrix(unsigned num_size, double** Matrix) {
    double max_min_in_triangular_matrix = .0;
    // lower

    double* Min = new double[num_size];
    for (int i = 0; i < num_size; ++i) {
        Min[i] = Matrix[i][0];
        for (int j = 0; j < i + 1; ++j) {
            if (Matrix[i][j] < Min[i]) {
                Min[i] = Matrix[i][j];
            }
        }
        if (Min[i] > max_min_in_triangular_matrix) {
            max_min_in_triangular_matrix = Min[i];
        }
    }
    delete[]Min;

    return max_min_in_triangular_matrix;
}

double MaxMinInTriangularMatrixParallel(unsigned num_size, double** Matrix) {
    double max_min_in_triangular_matrix_parallel = .0;
    // lower

    double* Min = new double[num_size];
    int i, j;

    #pragma omp parallel
    {
        #pragma omp for private(i,j)
        for (i = 0; i < num_size; ++i) {
            Min[i] = Matrix[i][0];
            for (j = 0; j < i + 1; ++j) {
                if (Matrix[i][j] < Min[i]) {
                    Min[i] = Matrix[i][j];
                }
            }

            #pragma omp flush(max_min_in_triangular_matrix_parallel)
            if (Min[i] > max_min_in_triangular_matrix_parallel) {
                
                #pragma omp critical
                if (Min[i] > max_min_in_triangular_matrix_parallel) {
                    max_min_in_triangular_matrix_parallel = Min[i];
                }
            }
        }
    }
    delete[]Min;

    return max_min_in_triangular_matrix_parallel;
}

double MaxMinInDiagonalMatrix(unsigned num_size, double** Matrix) {
    double max_min_in_diagonal_matrix = .0;

    return max_min_in_diagonal_matrix;
}

int main() {
    srand(time(0));
    const unsigned num_size = 10000;

    double** Matrix = new double* [num_size];
    for (int i = 0; i < num_size; ++i) {
        Matrix[i] = new double[num_size];
    }

    for (int i = 0; i < num_size; ++i) {
        for (int j = 0; j < num_size; ++j) {
            Matrix[i][j] = rand();
            //cout << setw(6) << Matrix[i][j];
        }
    }

    double start = omp_get_wtime();
    cout << "MaxMinInMatrix: " << MaxMinInMatrix(num_size, Matrix) << endl;
    double end = omp_get_wtime();
    cout << "Wtime: " << end - start << endl << endl;

    double start2 = omp_get_wtime();
    cout << "MaxMinInMatrixParallel: " << MaxMinInMatrixParallel(num_size, Matrix) << endl;
    double end2 = omp_get_wtime();
    cout << "Parallel Wtime: " << end2 - start2 << endl << endl;


    double start3 = omp_get_wtime();
    cout << "MaxMinInTriangularMatrix: " << MaxMinInTriangularMatrix(num_size, Matrix) << endl;
    double end3 = omp_get_wtime();
    cout << "Wtime: " << end3 - start3 << endl << endl;

    double start4 = omp_get_wtime();
    cout << "MaxMinInTriangularMatrixParallel: " << MaxMinInTriangularMatrixParallel(num_size, Matrix) << endl;
    double end4 = omp_get_wtime();
    cout << "Parallel Wtime: " << end4 - start4 << endl << endl;


    for (int i = 0; i < num_size; ++i) {
        delete[] Matrix[i];
    }
    delete[]Matrix;

	return 0;
}