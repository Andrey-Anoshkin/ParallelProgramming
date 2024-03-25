//MPI_Simpliest_c_Bind.cpp
#include "mpi.h"
#include <iostream>
#include <time.h>
#define PI 3.1415926535897932384626433832795

using namespace std;

int NProc, ProcId;
MPI_Status st;

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
	MPI_Barrier(MPI_COMM_WORLD);
	double sum = 0;
	int n = (int)((b - a) / h);

	for (int i = ProcId; i < n; i += NProc) {
		double x = a + i * h + h / 2;
		sum += f1(x) * h;
	}

	if (ProcId != 0) 
		MPI_Send(&sum, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	else {
		for (int i = 1; i < NProc; ++i) {
			double s;
			MPI_Recv(&s, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &st);
			sum += s;
		}

		*res = sum;
	}
}

void integral_Simpson(const double a, const double b, const double h, double* res) {
	MPI_Barrier(MPI_COMM_WORLD);
	double sum = f1(a) + f1(b);
	int n = (int)((b - a) / 2 / h);

	for (int i = ProcId + 1; i <= n; i += NProc) {
		double x = a + h * (2 * i - 1);
		sum += f1(x) * 4;
	}

	for (int i = ProcId + 1; i < n; i += NProc) {
		double x = a + 2 * i * h;
		sum += 2 * f1(x);
	}

	if (ProcId != 0)
		MPI_Send(&sum, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	else {
		for (int i = 1; i < NProc; ++i) {
			double s;
			MPI_Recv(&s, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &st);
			sum += s;
		}

		*res = h / 3 * sum;
	}
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
	int numbExp = 1;

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

	if (ProcId == 0) {
		cout << "Execution time: " << avg_time << "; " << min_time << "; " << max_time << "\n";
		cout.precision(8);
		cout << "Integral value: " << res << "\n";
	}
}

int main() {

	bool calcPosl = true, calcParal = true, calcSimpson = true;
	
	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &NProc);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcId);

	// Последовательное вычисление

	if (calcPosl && ProcId == 0) {
		cout << "Posl:\n";
		calculate(integral_posl);
		cout << "\n";
	}

	// Параллельное вычисление

	if (calcParal) {
		if (ProcId == 0)
			cout << "Paral:\n";
		calculate(integral_paral);
		if (ProcId == 0)
			cout << "\n";
	}

	// Параллельное вычисление по формуле Симпсона

	if (calcSimpson) {
		if (ProcId == 0)
			cout << "Simpson:\n";
		calculate(integral_Simpson);
		if (ProcId == 0)
			cout << "\n";
	}

	MPI_Finalize();

	return 0;
}
