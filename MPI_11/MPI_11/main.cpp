#include <iostream>
#include "mpi.h"
#include <time.h>
#include <cmath>
#define PI 3.1415926535897932384626433832795

using namespace std;

int NProc, ProcId;
MPI_Status st;

double f1(const double x, const double y, const double a1, const double b1, const double a2, const double b2) {
	return (1 + exp(sin(PI * x) * cos(PI * y))) / (b1 - a1) / (b2 - a2);
}

void integral_posl(const double a1, const double b1, const double a2, const double b2, const double h1, const double h2, double* res) {
	double sum = 0;
	int n = (int)((b1 - a1) / h1);
	int m = (int)((b2 - a2) / h2);

	for (int i = 0; i < n; ++i) {
		double x = a1 + i * h1 + h1 / 2;
		for (int j = 0; j < m; ++j) {
			double y = a2 + j * h2 + h2 / 2;
			sum += f1(x, y, a1, b1, a2, b2) * h1 * h2;
		}
	}

	*res = sum;
}

void integral_paral_1(const double a1, const double b1, const double a2, const double b2, const double h1, const double h2, double* res) {
	double sum = 0;
	int n = (int)((b1 - a1) / h1);
	int m = (int)((b2 - a2) / h2);

	for (int i = ProcId; i < n; i += NProc) {
		double x = a1 + i * h1 + h1 / 2;
		for (int j = 0; j < m; ++j) {
			double y = a2 + j * h2 + h2 / 2;
			sum += f1(x, y, a1, b1, a2, b2) * h1 * h2;
		}
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

void integral_paral_2(const double a1, const double b1, const double a2, const double b2, const double h1, const double h2, double* res) {
	double sum = 0;
	int n = (int)((b1 - a1) / h1);
	int m = (int)((b2 - a2) / h2);

	for (int i = ProcId / 2; i < n; i += NProc / 2) {
		double x = a1 + i * h1 + h1 / 2;
		for (int j = ProcId % 2; j < m; j += NProc / 2) {
			double y = a2 + j * h2 + h2 / 2;
			sum += f1(x, y, a1, b1, a2, b2) * h1 * h2;
		}
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

	if (ProcId == 0) {
		cout << "Execution time: " << avg_time << "; " << min_time << "; " << max_time << "\n";
		cout.precision(8);
		cout << "Integral value: " << res << "\n";
	}
}

int main() {

	bool calcPosl = true, calcParal_1 = true, calcParal_2 = true;

	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &NProc);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcId);

	// ���������������� ����������

	if (calcPosl && ProcId == 0) {
		cout << "Posl:\n";
		calculate(integral_posl);
		cout << "\n";
	}

	// ������������ ���������� (1)

	if (calcParal_1) {
		if (ProcId == 0)
			cout << "Paral_1:\n";
		calculate(integral_paral_1);
		if (ProcId == 0)
			cout << "\n";
	}

	// ������������ ���������� (2)

	if (calcParal_2) {
		if (ProcId == 0)
			cout << "Paral_2:\n";
		calculate(integral_paral_2);
		if (ProcId == 0)
			cout << "\n";
	}

	MPI_Finalize();

	return 0;
}