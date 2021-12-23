#include <iostream>
#include <omp.h>
#include <iomanip>
#include <Windows.h>

using namespace std;

/*f*/
double f(double x) {
	return x*x;
}

/*RiemannIntegral*/
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

/*RiemannIntegralParallel*/
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

/*TrapezoidalIntegral*/
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

/*TrapezoidalIntegralParallel*/
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

/*SimpsonIntegral*/
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

/*SimpsonIntegralParallel*/
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

	HANDLE  hout = GetStdHandle(STD_OUTPUT_HANDLE);
	COORD  size;
	size.X = 128;
	size.Y = 1001;
	SetConsoleScreenBufferSize(hout, size);

	double a = -1.;
	double b = 1.;
	const unsigned num_subintervals = 67108864;
	double start, end, result, time, start_parallel, end_parallel, result_parallel, time_parallel;

	omp_set_dynamic(0);
	int max_number_of_threads = omp_get_max_threads();
	for (int number_of_threads = 2; number_of_threads <= max_number_of_threads; ++number_of_threads) {
		omp_set_num_threads(number_of_threads);


		cout << endl;
		cout << endl;
		cout << "  number of threads = " << number_of_threads << endl;
		cout << "  -------------------------------------------------------------" << endl;
		cout << "                        Riemann Integral                       " << endl;
		cout << "  -------------------------------------------------------------" << endl;
		cout << "         i  |       times        |    s   |      results       " << endl;
		cout << "  -------------------------------------------------------------" << endl;
		for (int i = 2; i <= num_subintervals; i *= 2) {
			start = omp_get_wtime();
			result = RiemannIntegral(a, b, i);
			end = omp_get_wtime();
			time = end - start;

			start_parallel = omp_get_wtime();
			result_parallel = RiemannIntegralParallel(a, b, i);
			end_parallel = omp_get_wtime();
			time_parallel = end_parallel - start_parallel;

			cout << setw(10) << right << i << "  |  " \
				<< fixed << setprecision(4) << time << "    " \
				<< fixed << setprecision(4) << time_parallel << "  | " \
				<< setw(5) << fixed << setprecision(2) << (time) / (time_parallel) << "  |  " \
				<< fixed << setprecision(4) << result << "    " \
				<< fixed << setprecision(4) << result_parallel << endl;
		}

		cout << endl;
		cout << "                      Trapezoidal Integral                     " << endl;
		cout << "  -------------------------------------------------------------" << endl;
		cout << "         i  |       times        |    s   |      results       " << endl;
		cout << "  -------------------------------------------------------------" << endl;
		for (int i = 2; i <= num_subintervals; i *= 2) {
			start = omp_get_wtime();
			result = TrapezoidalIntegral(a, b, i);
			end = omp_get_wtime();
			time = end - start;

			start_parallel = omp_get_wtime();
			result_parallel = TrapezoidalIntegralParallel(a, b, i);
			end_parallel = omp_get_wtime();
			time_parallel = end_parallel - start_parallel;

			cout << setw(10) << right << i << "  |  " \
				<< fixed << setprecision(4) << time << "    " \
				<< fixed << setprecision(4) << time_parallel << "  | " \
				<< setw(5) << fixed << setprecision(2) << (time) / (time_parallel) << "  |  " \
				<< fixed << setprecision(4) << result << "    " \
				<< fixed << setprecision(4) << result_parallel << endl;
		}

		cout << endl;
		cout << "                        Simpson Integral                       " << endl;
		cout << "  -------------------------------------------------------------" << endl;
		cout << "         i  |       times        |    s   |      results       " << endl;
		cout << "  -------------------------------------------------------------" << endl;
		for (int i = 2; i <= num_subintervals; i *= 2) {
			start = omp_get_wtime();
			result = SimpsonIntegral(a, b, i);
			end = omp_get_wtime();
			time = end - start;

			start_parallel = omp_get_wtime();
			result_parallel = SimpsonIntegralParallel(a, b, i);
			end_parallel = omp_get_wtime();
			time_parallel = end_parallel - start_parallel;

			cout << setw(10) << right << i << "  |  " \
				<< fixed << setprecision(4) << time << "    " \
				<< fixed << setprecision(4) << time_parallel << "  | " \
				<< setw(5) << fixed << setprecision(2) << (time) / (time_parallel) << "  |  " \
				<< fixed << setprecision(4) << result << "    " \
				<< fixed << setprecision(4) << result_parallel << endl;
		}
	}
	return 0;
}