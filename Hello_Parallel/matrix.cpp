
#include <iostream>
using namespace std;
int main()
{
	int i, j, k, N, M, L;
	double** a, ** b, ** c;

	N = M = L = 1100;

	a = new double* [N];
	for (i = 0; i < N; i++) {
		a[i] = new double[M];
	}
	b = new double* [M];
	for (i = 0; i < M; i++) {
		b[i] = new double[L];
	}
	c = new double* [N];
	for (i = 0; i < N; i++) {
		c[i] = new double[L];
	}

	for (i = 0; i < N; i++) {
		for (j = 0; j < M; j++) {
			a[i][j] = rand() % 11;
		}
	}

	for (i = 0; i < M; i++) {
		for (j = 0; j < L; j++) {
			b[i][j] = rand() % 11;
		}
	}
#pragma omp parallel
{
#pragma omp parallel for shared(a, b, c) private(i, j, k)
	for (i = 0; i < N; i++) {
		for (j = 0; j < L; j++) {
			for (c[i][j] = 0, k = 0; k < M; k++) {
				c[i][j] += a[i][k] * b[k][j];
			}
		}
	}
}
	for (i = 0; i < N; i++) {
		delete[] a[i];
	}
	delete[]a;

	for (i = 0; i < M; i++) {
		delete[] b[i];
	}
	delete[] b;
	for (i = 0; i < N; i++) {
		delete[] c[i];
	}
	delete[]c;

	return 0;
}