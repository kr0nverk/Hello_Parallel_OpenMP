#include <iostream>
#include <iomanip>
#include <omp.h>
#include <ctime>

using namespace std;

/*MaxMinInMatrix*/
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

/*MaxMinInMatrixParallel*/
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

/*MaxMinInTriangularMatrix*/
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

/*MaxMinInTriangularMatrixParallel*/
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

/*MaxMinInDiagonalMatrix*/
double MaxMinInDiagonalMatrix(unsigned num_size, double** Matrix) {
    double max_min_in_diagonal_matrix = .0;

    return max_min_in_diagonal_matrix;
}

/*MaxMinInDiagonalMatrixParallel*/
double MaxMinInDiagonalMatrixParallel(unsigned num_size, double** Matrix) {
    double max_min_in_diagonal_matrix_parallel = .0;

    return max_min_in_diagonal_matrix_parallel;
}

int main() {
    srand(time(0));
    const unsigned num_size = 8192;

    double** Matrix = new double* [num_size];
    for (int i = 0; i < num_size; ++i) {
        Matrix[i] = new double[num_size];
    }

    for (int i = 0; i < num_size; ++i) {
        for (int j = 0; j < num_size; ++j) {
            Matrix[i][j] = rand();
            //cout << setw(6) << Matrix[i][j];
        }
        //cout << endl;
    }
    double start, end, result, time, start_parallel, end_parallel, result_parallel, time_parallel;

    cout << endl;
    cout << "                        Max Min In Matrix                      " << endl;
    cout << "  -------------------------------------------------------------" << endl;
    cout << "         i  |       times        |    a   |      results       " << endl;
    cout << "  -------------------------------------------------------------" << endl;
    for (int i = 1; i <= num_size; i *= 2) {
        start = omp_get_wtime();
        result = MaxMinInMatrix(i, Matrix);
        end = omp_get_wtime();
        time = end - start;

        start_parallel = omp_get_wtime();
        result_parallel = MaxMinInMatrixParallel(i, Matrix);
        end_parallel = omp_get_wtime();
        time_parallel = end_parallel - start_parallel;

        cout << setw(10) << right << i << "  |  " \
            << fixed << setprecision(4) << time << "    " \
            << fixed << setprecision(4) << time_parallel << "  |  " \
            << fixed << setprecision(2) << (time) / (time_parallel) << "  |  " \
            << setprecision(0) << result << "    " \
            << setprecision(0) << result_parallel << endl;
    }

    cout << endl;
    cout << "                   Max Min In Triangular Matrix                " << endl;
    cout << "  -------------------------------------------------------------" << endl;
    cout << "         i  |       times        |    a   |      results       " << endl;
    cout << "  -------------------------------------------------------------" << endl;
    for (int i = 1; i <= num_size; i *= 2) {
        start = omp_get_wtime();
        result = MaxMinInTriangularMatrix(i, Matrix);
        end = omp_get_wtime();
        time = end - start;

        start_parallel = omp_get_wtime();
        result_parallel = MaxMinInTriangularMatrixParallel(i, Matrix);
        end_parallel = omp_get_wtime();
        time_parallel = end_parallel - start_parallel;

        cout << setw(10) << right << i << "  |  " \
            << fixed << setprecision(4) << time << "    " \
            << fixed << setprecision(4) << time_parallel << "  |  " \
            << fixed << setprecision(2) << (time) / (time_parallel) << "  |  " \
            << fixed << setprecision(0) << result << "    " \
            << fixed << setprecision(0) << result_parallel << endl;
    }

    for (int i = 0; i < num_size; ++i) {
        delete[] Matrix[i];
    }
    delete[]Matrix;

	return 0;
}