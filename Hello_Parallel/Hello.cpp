#include <iostream>
#include <omp.h>

using namespace std;

int main() {

	#pragma omp parallel
	{
		#pragma omp critical
		cout << "Hello world\n";
	}

	#pragma omp parallel
	{
		#pragma omp atomic
		cout << "Hello world\n";
	}

	return 0;
}