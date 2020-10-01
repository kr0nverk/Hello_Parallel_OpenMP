#include <iostream>
#include <omp.h>
#include <vector>
#include <ctime>


using namespace std;


int main() {
	srand(time(0));
	vector<double> a;

	int num_subintervals = 10000000;

	for (int i = 0; i < num_subintervals; ++i) {
		a.push_back(rand()%100);
	}


	double start = omp_get_wtime();

	#pragma omp parallel
	{
		#pragma omp for
		for (int i = 0; i < num_subintervals; ++i) {
			a[i] += 1;
		}

	}

	double end = omp_get_wtime();
	cout << end - start << endl;

	double start2 = omp_get_wtime();

	for (int i = 0; i < num_subintervals; ++i) {
		a[i] += 1;
	}

	double end2 = omp_get_wtime();
	cout << end2 - start2 << endl;

	return 0;
}