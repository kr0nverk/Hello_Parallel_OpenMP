#include <iostream>
#include <omp.h>
#include <iomanip>

using namespace std;

/*NonParallel*/
double NonParallel(unsigned num_size, double* Matrix) {
	double sum = .0;

	for (unsigned i = 0; i < num_size; ++i) {
		sum += Matrix[i];
	}

	return sum;
}

/*Reduction*/
double Reduction(unsigned num_size, double* Matrix) {
	int i;
	double sum = .0;

	#pragma omp parallel for reduction(+:sum)
	for (i = 0; i < num_size; ++i) {
		sum += Matrix[i];
	}

	return sum;
}

/*Lock*/
double Lock(unsigned num_size, double* Matrix) {
	int i;
	double sum = .0;
	omp_lock_t myLock;
	omp_init_lock(&myLock);

	#pragma omp parallel for
	for (i = 0; i < num_size; ++i) {
		omp_set_lock(&myLock);
		sum += Matrix[i];
		omp_unset_lock(&myLock);
	}
	
	return sum;
}

/*Critical*/
double Critical(unsigned num_size, double* Matrix) {
	int i;
	double sum = .0;

	#pragma omp parallel for
	for (i = 0; i < num_size; ++i) {
		#pragma omp critical
		sum += Matrix[i];
	}

	return sum;
}

/*Atomic*/
double Atomic(unsigned num_size, double* Matrix) {
	int i;
	double sum = 0;

	#pragma omp parallel for
	for (i = 0; i < num_size; ++i) {
		#pragma omp atomic
		sum += Matrix[i];
	}

	return sum;
}

int main() {
	srand(time(0));
	const unsigned num_size = 10000000;

	double* Matrix = new double [num_size];
	for (int i = 0; i < num_size; ++i) {
		Matrix[i] = rand();
	}

	double start, end, result, result_reduction, result_lock, result_critical, result_atomic, time, time_reduction, time_lock, time_critical, time_atomic;

	cout << endl;
	cout << "                                                                             Reduction                                                               " << endl;
	cout << "  ---------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
	cout << "     i     |  time   time_r  time_l  time_c  time_a |   r     l     c     a   |                               results                             " << endl;
	cout << "  ---------------------------------------------------------------------------------------------------------------------------------------------------" << endl;

	for (int i = 500000; i <= num_size; i += 500000) {

		start = omp_get_wtime();
		result = NonParallel(i, Matrix);
		end = omp_get_wtime();
		time = end - start;

		start = omp_get_wtime();
		result_reduction = Reduction(i, Matrix);
		end = omp_get_wtime();
		time_reduction = end - start;

		start = omp_get_wtime();
		result_lock = Lock(i, Matrix);
		end = omp_get_wtime();
		time_lock = end - start;

		start = omp_get_wtime();
		result_critical = Critical(i, Matrix);
		end = omp_get_wtime();
		time_critical = end - start;

		start = omp_get_wtime();
		result_atomic = Atomic(i, Matrix);
		end = omp_get_wtime();
		time_atomic = end - start;

		cout << setw(10) << right << i << " | " \
			<< fixed << setprecision(4) << time << "  " \
			<< fixed << setprecision(4) << time_reduction << "  " \
			<< fixed << setprecision(4) << time_lock << "  " \
			<< fixed << setprecision(4) << time_critical << "  " \
			<< fixed << setprecision(4) << time_atomic << " | " \
			<< setw(5) << fixed << setprecision(2) << (time) / (time_reduction) << " " \
			<< setw(5) << fixed << setprecision(2) << (time) / (time_lock) << " " \
			<< setw(5) << fixed << setprecision(2) << (time) / (time_critical) << " " \
			<< setw(5) << fixed << setprecision(2) << (time) / (time_atomic) << " | " \
			<< setprecision(0) << result << " " \
			<< setprecision(0) << result_reduction << " " \
			<< setprecision(0) << result_lock << " " \
			<< setprecision(0) << result_critical << " " \
			<< setprecision(0) << result_atomic << " " << endl;

	}
	return 0;
}