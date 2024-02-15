#include <iostream>
#include <omp.h>
#include <time.h>
#define PI 3.1415926535897932384626433832795

using namespace std;

double f1(double x) {
	return 1.0 / (1 + x * x);
}

void integral_posl(const double a, const double b, const double h, double* res) {
	double sum = 0;
	int n = (int)((b - a) / h);

	for (int i = 0; i < n; ++i) {
		double x = a + i * h + h / 2;
		sum += f1(x) * h;
	}

	*res = sum;
}

void integral_paral(const double a, const double b, const double h, double* res) {
	double sum = 0;
	int n = (int)((b - a) / h);

	#pragma omp parallel for reduction(+: sum)
	for (int i = 0; i < n; ++i) {
		double x = a + i * h + h / 2;
		sum += f1(x) * h;
	}

	*res = sum;
}

void integral_Simpson(const double a, const double b, const double h, double* res) {
	double sum = f1(a) + f1(b);
	int n = (int)((b - a) / 2 / h);

	#pragma omp parallel for reduction(+: sum)
	for (int i = 1; i <= n; ++i) {
		double x = a + h * (2 * i - 1);
		sum += f1(x) * 4;
	}

	#pragma omp parallel for reduction(+: sum)
	for (int i = 1; i < n; ++i) {
		double x = a + 2 * i * h;
		sum += 2 * f1(x);
	}

	*res = h / 3 * sum;
}

double experiment(double* res, void f(const double a, const double b, const double h, double* res)) {
	double stime = 0, ftime = 0;
	double a = 0;
	double b = 1e+6;
	double h = 0.01;
	stime = clock();
	f(a, b, h, res);
	ftime = clock();

	return (ftime - stime) / CLOCKS_PER_SEC;
}

void calculate(void f(const double a, const double b, const double h, double* res)) {
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
	
	bool calcPosl = true, calcParal = true, calcSimpson = true;

	// Последовательное вычисление

	if (calcPosl) {
		cout << "Posl:\n";
		calculate(integral_posl);
		cout << "\n";
	}

	// Параллельное вычисление
	
	if (calcParal) {
		cout << "Paral:\n";
		calculate(integral_paral);
		cout << "\n";
	}

	// Параллельное вычисление по формуле Симпсона

	if (calcSimpson) {
		cout << "Simpson:\n";
		calculate(integral_Simpson);
		cout << "\n";
	}

	system("pause");
	return 0;
}