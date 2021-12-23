#include <iostream>
#include <omp.h>
#include <iomanip>
#include <Windows.h>

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
	HANDLE  hout = GetStdHandle(STD_OUTPUT_HANDLE);
	COORD  size;
	size.X = 128;
	size.Y = 1001;
	SetConsoleScreenBufferSize(hout, size);

	srand(time(0));
	const unsigned num_size = 1000000;
	int sum = 0;
	double* Matrix = new double [num_size];

	cout << endl;
	cout << "                                 Static Dynamic Guided                                 " << endl;
	cout << "  -------------------------------------------------------------------------------------" << endl;

	omp_set_dynamic(0);
	int max_number_of_threads = omp_get_max_threads();
	for (int number_of_threads = 2; number_of_threads <= max_number_of_threads; ++number_of_threads) {
		omp_set_num_threads(number_of_threads);

		cout << endl;
		cout << endl;
		cout << "  number of threads = " << number_of_threads << endl;

		double start = omp_get_wtime();
		NonParsllelMatrix(num_size, Matrix);
		double end = omp_get_wtime();
		double time = end - start;

		for (int i = 0; i < num_size; ++i) {
			sum = sum + Matrix[i];
		}

		start = omp_get_wtime();
		ParallelFor(num_size, Matrix);
		end = omp_get_wtime();
		double parallel_time = end - start;

		for (int i = 0; i < num_size; ++i) {
			sum = sum + Matrix[i];
		}


		double static_sum = .0;
		double dynamic_sum = .0;
		double guided_sum = .0;

		cout << "  time / parallel_time = " << fixed << setprecision(2) << (time) / (parallel_time) << endl;
		cout << "  -------------------------------------------------------------------------------------" << endl;
		cout << "   chunk |   s_s     s_d     s_g  |                       results                      " << endl;
		cout << "  -------------------------------------------------------------------------------------" << endl;

		for (int chunk = 1; chunk <= num_size / 8; chunk *= 2) {

			start = omp_get_wtime();
			ParallelForStatic(num_size, Matrix, chunk);
			end = omp_get_wtime();
			double static_time = end - start;

			for (int i = 0; i < num_size; ++i) {
				static_sum = static_sum + Matrix[i];
			}

			start = omp_get_wtime();
			ParallelForDynamic(num_size, Matrix, chunk);
			end = omp_get_wtime();
			double dynamic_time = end - start;

			for (int i = 0; i < num_size; ++i) {
				dynamic_sum = dynamic_sum + Matrix[i];
			}

			start = omp_get_wtime();
			ParallelForGuided(num_size, Matrix, chunk);
			end = omp_get_wtime();
			double guided_time = end - start;

			for (int i = 0; i < num_size; ++i) {
				guided_sum = guided_sum + Matrix[i];
			}

			cout << setprecision(5);
			cout << setw(7) << right << chunk << "  |  " \
				<< fixed << setprecision(2) << (time) / (static_time) << "    " \
				<< fixed << setprecision(2) << (time) / (dynamic_time) << "    " \
				<< fixed << setprecision(2) << (time) / (guided_time) << "  |  " \
				<< fixed << setprecision(0) << static_sum << " " \
				<< fixed << setprecision(0) << dynamic_sum << " " \
				<< fixed << setprecision(0) << guided_sum << endl;
		}
	}
	/*
	sum = 0;
	cout << endl;
	cout << "               For dynamic            " << endl;
	cout << "  ------------------------------------" << endl;
	cout << "  chunk - time           result       " << endl;
	cout << "  ------------------------------------" << endl;
	for (int chunk = 1; chunk <= num_size / 8; chunk *= 2) {
		
		start = omp_get_wtime();
		ParallelForDynamic(num_size, Matrix, chunk);
		end = omp_get_wtime();
		
		for (int i = 0; i < num_size; ++i) {
			sum = sum + Matrix[i];
		}
		cout << setprecision(5);
		cout << setw(7) << right << chunk << " - " << setw(11) << left << time / (end - start) << " " << sum << endl;
	}

	sum = 0;
	cout << endl;
	cout << "               For guided             " << endl;
	cout << "  ------------------------------------" << endl;
	cout << "  chunk - time           result       " << endl;
	cout << "  ------------------------------------" << endl;
	for (int chunk = 1; chunk <= num_size / 8; chunk *= 2) {
		
		start = omp_get_wtime();
		ParallelForGuided(num_size, Matrix, chunk);
		end = omp_get_wtime();

		for (int i = 0; i < num_size; ++i) {
			sum = sum + Matrix[i];
		}
		cout << setprecision(5);
		cout << setw(7) << right << chunk << " - " << setw(11) << left << time / (end - start) << " " << sum << endl;
	}
	*/

	delete[]Matrix;

	return 0;
}