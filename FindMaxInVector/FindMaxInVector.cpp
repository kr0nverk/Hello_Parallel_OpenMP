#include <iostream>
#include <vector>
#include <ctime>
#include <omp.h>

using namespace std;

double FindMaxInVector(vector<int> temp) {
    int find_max_in_vector = -1;

    for (int i = 0; i < temp.size(); ++i) {
        if (temp[i] > find_max_in_vector) {
            find_max_in_vector = temp[i];
        }
    }

    return find_max_in_vector;
}

double FindMaxInVectorParallel(vector<int> temp) {
    int find_max_in_vector_parallel = -1;
    int i;
#pragma omp parallel
{
    #pragma omp for private(i)
    for (i = 0; i < temp.size(); ++i) {

        #pragma omp flush(find_max_in_vector_parallel)
        if (temp[i] > find_max_in_vector_parallel) {

            #pragma omp critical
            if (temp[i] > find_max_in_vector_parallel) {
                    find_max_in_vector_parallel = temp[i];
            }
        }
    }
}
    return find_max_in_vector_parallel;
}

int main() {
    srand(time(0));
    vector<int> temp;
    for (int i = 0; i < 10000000; ++i) {
        temp.push_back(rand());
    }

    double start;
    double end;
    start = omp_get_wtime();
    cout << FindMaxInVector(temp) << endl;
    end = omp_get_wtime();
    cout << end - start << endl;
    
    double start2;
    double end2;
    start2 = omp_get_wtime();
    cout << FindMaxInVectorParallel(temp) << endl;
    end2 = omp_get_wtime();
    cout << end2 - start2 << endl;

    return 0;
}