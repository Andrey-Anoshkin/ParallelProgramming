#include <iostream>
#include <omp.h>
#include <time.h>
#include <cmath>
#define PI 3.1415926535897932384626433832795

using namespace std;

double f1(const double x, const double y, const double a1, const double b1, const double a2, const double b2) {
	return (1 + exp(sin(PI * x) * cos(PI * y))) / (b1 - a1) / (b2 - a2);
}

double q(int i, int j, int n, int m) {
	double q = 1;
	if (i == 0 || i == n)
		q *= 0.5;
	if (j == 0 || j == m)
		q *= 0.5;
	return q;
}

void integral(const double a1, const double b1, const double a2, const double b2, const double h1, const double h2, double* res) {
	double sum = 0;
	int n = (int)((b1 - a1) / h1);
	int m = (int)((b2 - a2) / h2);

	#pragma omp parallel for reduction(+: sum)
	for (int i = 0; i <= n; ++i) {
		double x = a1 + i * h1 + h1 / 2;
		for (int j = 0; j <= m; ++j) {
			double y = a2 + j * h2 + h2 / 2;
			sum += f1(x, y, a1, b1, a2, b2) * h1 * h2 * q(i, j, n, m);
		}
	}

	*res = sum;
}

double experiment(double* res, void f(const double a1, const double b1, const double a2, const double b2, const double h1, const double h2, double* res)) {
	double stime = 0, ftime = 0;
	double a1 = 0, b1 = 16;
	double a2 = 0, b2 = 16;
	double h1 = 0.005, h2 = 0.005;
	stime = clock();
	f(a1, b1, a2, b2, h1, h2, res);
	ftime = clock();

	return (ftime - stime) / CLOCKS_PER_SEC;
}

void calculate(void f(const double a1, const double b1, const double a2, const double b2, const double h1, const double h2, double* res)) {
	double avg_time = 0;
	double min_time = 0;
	double max_time = 0;
	double res = 0;
	int numbExp = 10;

	min_time = max_time = experiment(&res, f);
	avg_time = min_time / numbExp;
	for (int i = 0; i < numbExp - 1; ++i) {
		double time = experiment(&res, f);
		if (time > max_time)
			max_time = time;
		if (time < min_time)
			min_time = time;

		avg_time += time / numbExp;
	}

	cout << "Execution time: " << avg_time << "; " << min_time << "; " << max_time << "\n";
	cout.precision(8);
	cout << "Integral value: " << res << "\n";
}

int main() {

	calculate(integral); 
	cout << "\n";

	system("pause");
	return 0;
}