#include <iostream>
#include <omp.h>

using namespace std;

int main() {
	int number_of_threads = 4;
	omp_set_dynamic(0);
	omp_set_num_threads(number_of_threads);
	#pragma omp parallel
	{
		#pragma omp critical
		cout << "Hello world\n";
	}

	return 0;
}