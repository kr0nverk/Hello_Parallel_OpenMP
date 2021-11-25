#include <iostream>
#include <vector>
#include <ctime>
#include <omp.h>
#include <iomanip>

using namespace std;

/*ScalarProduct*/
double ScalarProduct(vector<double> a, vector<double> b, int n) {
    double scalar_product = 0;

    for (int i = 0; i < n; ++i) {
        scalar_product += (a[i]) * (b[i]);
    }

    return scalar_product;
}

/*ScalarProductParallel*/
double ScalarProductParallel(vector<double> a, vector<double> b, int n) {
    double scalar_product_parallel = 0;
    int i;
    double sum = 0;

    #pragma omp parallel shared(sum)
    {
        #pragma omp for private(i)
        for (i = 0; i < n; ++i) {
            sum += (double) a[i] * b[i];
        }
        #pragma omp critical
        {
            scalar_product_parallel += sum;
        }

    }

    return scalar_product_parallel;
}

int main() {
    srand(time(0));
    vector<double> veca;
    vector<double> vecb;

    int n = 10000000;

    for (int i = 0; i < n; ++i) {
        veca.push_back(rand()%9);
        vecb.push_back(rand()%9);
    }
    double start, end, time, result, start_parallel, end_parallel, time_parallel, result_parallel;

    cout << endl;
    cout << "                          ScalarProduct                        " << endl;
    cout << "  -------------------------------------------------------------" << endl;
    cout << "         i  |       times        |    a   |      results       " << endl;
    cout << "  -------------------------------------------------------------" << endl;
    for (int i = 1000000; i <= n; i += 1000000) {
        start = omp_get_wtime();
        result = ScalarProduct(veca, vecb, i);
        end = omp_get_wtime();
        time = end - start;

        start_parallel = omp_get_wtime();
        result_parallel = ScalarProductParallel(veca, vecb, i);
        end_parallel = omp_get_wtime();
        time_parallel = end_parallel - start_parallel;

        cout << setw(10) << right << i << "  |  " \
            << fixed << setprecision(4) << time << "    " \
            << fixed << setprecision(4) << time_parallel << "  |  " \
            << fixed << setprecision(2) << (time) / (time_parallel) << "  |  " \
            << fixed << setprecision(0) << result << "    " \
            << fixed << setprecision(0) << result_parallel << endl;
    }

    return 0;
}