#include <iostream>
#include <omp.h>

int main() {

	#pragma omp parallel
	{
		#pragma omp critical
		std::cout << "Hello world\n";
	}
	return 0;
}