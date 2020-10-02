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
    for (int i = 0; i <= a.size() - 1; ++i) {
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

#pragma omp parallel
{
#pragma omp for private(i) shared(a, b, scalar_product_parallel)
        for (int i = 0; i < a.size(); ++i) {
            scalar_product_parallel += a[i] * b[i];
        }
}
    return scalar_product_parallel;
}

int main() {
    srand(time(0));
    vector<double> veca;
    vector<double> vecb;

    for (int i = 0; i < 1000000; ++i) {
        veca.push_back(rand());
        vecb.push_back(rand());
    }

    double start;
    double end;
    start = omp_get_wtime();
    cout << ScalarProduct(veca, vecb) << endl;
    end = omp_get_wtime();
    cout << end - start << endl;

    double start2;
    double end2;
    start2 = omp_get_wtime();
    cout << ScalarProductParallel(veca, vecb) << endl;
    end2 = omp_get_wtime();
    cout << end2 - start2 << endl;
    return 0;
}