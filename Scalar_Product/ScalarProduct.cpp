#include <iostream>
#include <vector>
#include <ctime>
#include <omp.h>

using namespace std;

double ScalarProduct(vector<double> a, vector<double> b) {
    if (a.size() != b.size()) {
        cout << "different sizes";
        return -1;
    }

    double scalar_product = 0;

    for (int i = 0; i < a.size(); ++i) {
        scalar_product += (a[i]) * (b[i]);
    }

    return scalar_product;
}

double ScalarProductParallel(vector<double> a, vector<double> b) {
    if (a.size() != b.size()) {
        cout << "different sizes";
        return -1;
    }

    double scalar_product_parallel = 0;
    int i;

    #pragma omp parallel
    {
        #pragma omp for private(i)
        for (i = 0; i < a.size(); ++i) {
            scalar_product_parallel += a[i] * b[i];
        }
    }

    return scalar_product_parallel;
}

int main() {
    srand(time(0));
    vector<double> veca;
    vector<double> vecb;

    for (int i = 0; i < 1000; ++i) {
        veca.push_back(rand());
        vecb.push_back(rand());
    }

    double start = omp_get_wtime();
    cout << "Scalar product: " << ScalarProduct(veca, vecb) << endl;
    double end = omp_get_wtime();
    cout << "Wtime: " << end - start << endl;

    double start2 = omp_get_wtime();
    cout << "Parallel Scalar product: " << ScalarProductParallel(veca, vecb) << endl;
    double end2 = omp_get_wtime();
    cout << "Parallel Wtime: " << end2 - start2 << endl;

    return 0;
}