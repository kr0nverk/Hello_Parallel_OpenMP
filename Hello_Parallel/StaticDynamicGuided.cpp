#include <iostream>
#include <omp.h>
#include <iomanip>
using namespace std;

void NonParsllelMatrix(unsigned num_size, double* Matrix) {
	for (int i = 0; i < num_size; ++i) {
		if (i < num_size / 2) {
			Matrix[i] = (double)i;
		}
		else {
			Matrix[i] = rand();
			Matrix[i] += rand();
			Matrix[i] += rand();
			Matrix[i] += rand();
			Matrix[i] += rand();
		}
	}
}

void ParallelFor(unsigned num_size, double* Matrix) {
	int i;
#pragma omp parallel
#pragma omp for private(i)
	for (i = 0; i < num_size; ++i) {
		if (i < num_size / 2) {
			Matrix[i] = (double)i;
		}
		else {
			Matrix[i] = rand();
			Matrix[i] += rand();
			Matrix[i] += rand();
			Matrix[i] += rand();
			Matrix[i] += rand();
		}
	}
}

void ParallelForStatic(unsigned num_size, double* Matrix, int chunk) {
	int i;
#pragma omp parallel
#pragma omp for schedule(static, chunk) nowait private(i)
	for (i = 0; i < num_size; ++i) {
		if (i < num_size / 2) {
			Matrix[i] = (double)i;
		}
		else {
			Matrix[i] = rand();
			Matrix[i] += rand();
			Matrix[i] += rand();
			Matrix[i] += rand();
			Matrix[i] += rand();
		}
	}
}

void ParallelForDynamic(unsigned num_size, double* Matrix, int chunk) {
	int i;
#pragma omp parallel
#pragma omp for schedule(dynamic, chunk) nowait private(i)
	for (i = 0; i < num_size; ++i) {
		if (i < num_size / 2) {
			Matrix[i] = (double)i;
		}
		else {
			Matrix[i] = rand();
			Matrix[i] += rand();
			Matrix[i] += rand();
			Matrix[i] += rand();
			Matrix[i] += rand();
		}
	}
}

void ParallelForGuided(unsigned num_size, double* Matrix, int chunk) {
	int i;
#pragma omp parallel
#pragma omp for schedule(guided, chunk) nowait private(i)
	for (i = 0; i < num_size; ++i) {
		if (i < num_size / 2) {
			Matrix[i] = (double)i;
		}
		else {
			Matrix[i] = rand();
			Matrix[i] += rand();
			Matrix[i] += rand();
			Matrix[i] += rand();
			Matrix[i] += rand();
		}
	}
}

int main()
{
	srand(time(0));
	const unsigned num_size = 100000000;
	int sum = 0;
	double* Matrix = new double [num_size];

	sum = 0;
	cout << "N: " << num_size << endl;
	cout << "For" << endl;
	
	double start = omp_get_wtime();
	NonParsllelMatrix(num_size, Matrix);
	double end = omp_get_wtime();
	
	for (int i = 0; i < num_size; ++i) {
		sum = sum + Matrix[i];
	}
	cout << setprecision(5);
	cout << setw(10) << left << "Wtime: " << setw(11) << left << end - start <<" "<< sum << endl << endl;

	int i;
	int chunk;

	sum = 0;
	cout << "For parallel" << endl;
	
	start = omp_get_wtime();
	ParallelFor(num_size, Matrix);
	end = omp_get_wtime();
	
	for (int i = 0; i < num_size; ++i) {
		sum = sum + Matrix[i];
	}
	cout << setprecision(5);
	cout << setw(10) << left << "Wtime: " << setw(11) << left << end - start << " " << sum << endl << endl;


	sum = 0;
	cout << "For static" << endl;
	cout << setw(6) << right << " chunk" << setw(11) << left << " - Wtime:" << endl;
	for (int chunk = 1; chunk <= num_size/8; chunk *= 2) {
		
		start = omp_get_wtime();
		ParallelForStatic(num_size, Matrix, chunk);
		end = omp_get_wtime();
		
		for (int i = 0; i < num_size; ++i) {
			sum = sum + Matrix[i];
		}
		cout << setprecision(5);
		cout << setw(7) << right << chunk << " - " << setw(11) << left << end - start << " " << sum << endl;
	}
	cout << endl;


	sum = 0;
	cout << "For dynamic" << endl;
	cout << setw(6) << right << " chunk" << setw(11) << left << " - Wtime:" << endl;
	for (int chunk = 1; chunk <= num_size / 8; chunk *= 2) {
		
		start = omp_get_wtime();
		ParallelForDynamic(num_size, Matrix, chunk);
		end = omp_get_wtime();
		
		for (int i = 0; i < num_size; ++i) {
			sum = sum + Matrix[i];
		}
		cout << setprecision(5);
		cout << setw(7) << right << chunk << " - " << setw(11) << left << end - start << " " << sum << endl;
	}
	cout << endl;

	sum = 0;
	cout << "For guided" << endl;
	cout << setw(6) << right << " chunk" << setw(11) << left << " - Wtime:" << endl;
	for (int chunk = 1; chunk <= num_size / 8; chunk *= 2) {
		
		start = omp_get_wtime();
		ParallelForGuided(num_size, Matrix, chunk);
		end = omp_get_wtime();

		for (int i = 0; i < num_size; ++i) {
			sum = sum + Matrix[i];
		}
		cout << setprecision(5);
		cout << setw(7) << right << chunk << " - " << setw(11) << left << end - start << " " << sum << endl;
	}
	cout << endl;


	delete[]Matrix;

	return 0;
}