#include <iostream>

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

int main() {
	double a = -1.;
	double b = 1.;
	unsigned num_subintervals = 100;

	cout << "RiemannIntegral " << RiemannIntegral(a, b, num_subintervals) << endl;
	cout << "TrapezoidalIntegral " << TrapezoidalIntegral(a, b, num_subintervals) << endl;
	cout << "SimpsonIntegral " << SimpsonIntegral(a, b, num_subintervals) << endl;
	
	
	return 0;
}