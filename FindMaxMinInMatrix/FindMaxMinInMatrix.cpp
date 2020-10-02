#include <iostream>
#include <iomanip>

using namespace std;

double MaxMinInMatrix(unsigned N, double** a) {
    double max_min_in_matrix = .0;
    
    double* min = new double[N];
    for (int i = 0; i < N; ++i) {
        min[i] = a[i][0];
        for (int j = 0; j < N; ++j) {
            if (a[i][j] < min[i]) {
                min[i] = a[i][j];
            }
        }
        max_min_in_matrix = min[0];
        if (min[i] > max_min_in_matrix) {
            max_min_in_matrix = min[i];
        }
    }

    return max_min_in_matrix;
}

double MaxMinInTriangularMatrix(unsigned N, double** a) {
    double max_min_in_triangular_matrix = .0;
    // lower

    double* min = new double[N];
    for (int i = 0; i < N; ++i) {
        min[i] = a[i][0];
        for (int j = N - i - 1; j < N; ++j) {
            if (a[i][j] < min[i]) {
                min[i] = a[i][j];
            }
        }
        max_min_in_triangular_matrix = min[0];
        if (min[i] > max_min_in_triangular_matrix) {
            max_min_in_triangular_matrix = min[i];
        }
    }

    return max_min_in_triangular_matrix;
}

double MaxMinInDiagonalMatrix(unsigned N, double** a) {
    double max_min_in_diagonal_matrix = .0;

    return max_min_in_diagonal_matrix;
}

int main() {
    const unsigned N = 10;

    double** a = new double* [N];
    for (int i = 0; i < N; ++i) {
        a[i] = new double[N];
    }

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            a[i][j] = i + j;
            cout << setw(4) << a[i][j];
        }
        cout << endl;
    }

    cout << "MaxMinInMatrix " << MaxMinInMatrix(N, a) << endl;
    cout << "MaxMinInTriangularMatrix " << MaxMinInTriangularMatrix(N, a) << endl;


    for (int i = 0; i < N; ++i) {
        for (int j = N - i - 1; j < N; ++j) {
            a[i][j] = i;
            cout << setw(4) << a[i][j];
        }
        cout << endl;
    }


    for (int i = 0; i < N; ++i) {
        delete[] a[i];
    }
    delete[]a;

	return 0;
}