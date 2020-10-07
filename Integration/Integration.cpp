#include <iostream>
#include <omp.h>

using namespace std;

double f(double x) {
	return x*x;
}

double RiemannIntegral(double a, double b, unsigned num_subintervals) {
	const double width = (b - a) / num_subintervals;

	double riemann_integral = 0.;
	for (int step = 0; step < num_subintervals; ++step) {
		const double x = a + width * step; //sample at left endpoint
		// const double x = a + width * step + width/2; //sample at midpoint
		// const double x = a + width * step + width; //sample at right endpoint

		riemann_integral += f(x) * width;
	}
	
	return riemann_integral;
}

double RiemannIntegralParallel(double a, double b, unsigned num_subintervals) {
	const double width = (b - a) / num_subintervals;

	double riemann_integral_parallel = 0.;
	int step;

	#pragma omp parallel
	{
		#pragma omp for private(step) reduction(+:riemann_integral_parallel)
		for (step = 0; step < num_subintervals; ++step) {
			const double x = a + width * step; //sample at left endpoint
			// const double x = a + width * step + width/2; //sample at midpoint
			// const double x = a + width * step + width; //sample at right endpoint

			riemann_integral_parallel += f(x) * width;
		}
	}

	return riemann_integral_parallel;
}

double TrapezoidalIntegral(double a, double b, unsigned num_subintervals) {
	const double width = (b - a) / num_subintervals;

	double trapezoidal_integral = 0.;
	for (int step = 0; step < num_subintervals; ++step) {
		const double x1 = a + step * width;
		const double x2 = a + (step + 1) * width;

		trapezoidal_integral += 0.5 * (x2 - x1) * (f(x1) + f(x2));
	}

	return trapezoidal_integral;
}

double TrapezoidalIntegralParallel(double a, double b, unsigned num_subintervals) {
	const double width = (b - a) / num_subintervals;

	double trapezoidal_integral_parallel = 0.;
	int step;
	#pragma omp parallel
	{
		#pragma omp for private(step) reduction(+:trapezoidal_integral_parallel)
		for (step = 0; step < num_subintervals; ++step) {
			const double x1 = a + step * width;
			const double x2 = a + (step + 1) * width;

			trapezoidal_integral_parallel += 0.5 * (x2 - x1) * (f(x1) + f(x2));
		}
	}
	return trapezoidal_integral_parallel;
}

double SimpsonIntegral(double a, double b, unsigned num_subintervals) {
	const double width = (b - a) / num_subintervals;

	double simpson_integral = 0.;
	for (int step = 0; step < num_subintervals; ++step) {
		const double x1 = a + step * width;
		const double x2 = a + (step + 1) * width;

		simpson_integral += (x2 - x1) / 6.0 * (f(x1) + 4.0 * f(0.5 * (x1 + x2)) + f(x2));
	}

	return simpson_integral;
}

double SimpsonIntegralParallel(double a, double b, unsigned num_subintervals) {
	const double width = (b - a) / num_subintervals;

	double simpson_integral_parallel = 0.;
	int step;

	#pragma omp parallel
	{
		#pragma omp for private(step) reduction(+:simpson_integral_parallel)
		for (int step = 0; step < num_subintervals; ++step) {
			const double x1 = a + step * width;
			const double x2 = a + (step + 1) * width;

			simpson_integral_parallel += (x2 - x1) / 6.0 * (f(x1) + 4.0 * f(0.5 * (x1 + x2)) + f(x2));
		}
	}
	return simpson_integral_parallel;
}


int main() {
	double a = -1.;
	double b = 1.;
	const unsigned num_subintervals = 10000000;

	double start = omp_get_wtime();
	cout << "RiemannIntegral: " << RiemannIntegral(a, b, num_subintervals) << endl;
	double end = omp_get_wtime();
	cout << "Wtime: " << end - start << endl << endl;

	double start2 = omp_get_wtime();
	cout << "RiemannIntegralParallel: " << RiemannIntegralParallel(a, b, num_subintervals) << endl;
	double end2 = omp_get_wtime();
	cout << "Parallel Wtime: " << end2 - start2 << endl << endl;

	double start3 = omp_get_wtime();
	cout << "TrapezoidalIntegral: " << TrapezoidalIntegral(a, b, num_subintervals) << endl;
	double end3 = omp_get_wtime();
	cout << "Wtime: " << end3 - start3 << endl << endl;

	double start4 = omp_get_wtime();
	cout << "TrapezoidalIntegralParallel: " << TrapezoidalIntegralParallel(a, b, num_subintervals) << endl;
	double end4 = omp_get_wtime();
	cout << "Parallel Wtime: " << end4 - start4 << endl << endl;

	double start5 = omp_get_wtime();
	cout << "SimpsonIntegral: " << SimpsonIntegral(a, b, num_subintervals) << endl;
	double end5 = omp_get_wtime();
	cout << "Wtime: " << end5 - start5 << endl << endl;

	double start6 = omp_get_wtime();
	cout << "SimpsonIntegralParallel: " << SimpsonIntegralParallel(a, b, num_subintervals) << endl;
	double end6 = omp_get_wtime();
	cout << "Parallel Wtime: " << end6 - start6 << endl << endl;
	
	return 0;
}