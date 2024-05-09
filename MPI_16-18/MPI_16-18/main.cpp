#include <iomanip>
#include <iostream>
#include <cmath>
#include <complex>
#include <time.h>
#include <vector>
#include "mpi.h"

#define PI 3.14159265358979323846

using namespace std;

int ProcId, NProc;
MPI_Status st;

//Function for simple initialization of input signal elements
void DummyDataInitialization(complex<double>* input, int size) {
	for (int i = 0; i < size; ++i)
		input[i] = 0;
	input[size - size / 4] = 1;
}

// Function for random initialization of objects' elements
void RandomDataInitialization(complex<double>* input, int size) {
	srand(unsigned(clock()));
	for (int i = 0; i < size; ++i)
		input[i] = complex<double>(rand() / 1000.0, rand() / 1000.0);
}

void MyDataInitialization(complex<double>* input, int size) {
	int task = 2;
	int partialSize = size / NProc;
	complex<double>* partialInput = new complex<double>[partialSize];

	for (int i = 0; i < partialSize; ++i) {
		switch (task) {
		case 1: 
			if (!ProcId && !i)
				partialInput[0] = sin(2 * PI * 10 * 1.0 / size * 0); 
			else
				partialInput[i] = sin(2 * PI * 10 * 1.0 / size * (i + ProcId * partialSize)); 
			break;
		case 2: 
			if (!ProcId && !i)
				partialInput[0] = 1; 
			else
				partialInput[i] = -log(2 * sin(PI * 1.0 / size * (i + ProcId * partialSize)));
			break;
		}
	}

	MPI_Gather(partialInput, partialSize, MPI_C_DOUBLE_COMPLEX, input, partialSize, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
	MPI_Bcast(input, size, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
}

//Function for memory allocation and data initialization
void ProcessInitialization(complex<double>*& inputSignal, complex<double>*& outputSignal, int& size) {

	size = 32;

	cout << "Input signal length = " << size << endl;
	inputSignal = new complex<double>[size];
	outputSignal = new complex<double>[size];

	//Initialization of input signal elements - tests
	//DummyDataInitialization(inputSignal, size);

	//Computational experiments
	//RandomDataInitialization(inputSignal, size);

	MyDataInitialization(inputSignal, size);
}

void BitReversing(complex<double>* inputSignal, complex<double>* outputSignal, int size) {
	int bitsCount = 0;
	//bitsCount = log2(size)
	for (int tmp_size = size; tmp_size > 1; tmp_size /= 2, bitsCount++);

	//ind - index in input array
	//revInd - correspondent to ind index in output array

	for (int ind = 0; ind < size; ind++) {
		int mask = 1 << (bitsCount - 1);
		int revInd = 0;
		for (int i = 0; i < bitsCount; i++) {
			bool val = ind & mask;
			revInd |= val << i;
			mask = mask >> 1;
		}
		outputSignal[revInd] = inputSignal[ind];
	}
}

void ParallelBitReversing(complex<double>* inputSignal, complex<double>* outputSignal, int size) {
	int bitsCount = 0;
	int partialSize = size / NProc;
	int counter = 0;
	complex<double>* partialOutputSignal = new complex<double>[partialSize];
	//bitsCount = log2(size)
	for (int tmp_size = size; tmp_size > 1; tmp_size /= 2, ++bitsCount);

	//ind - index in input array
	//revInd - correspondent to ind index in output array

	for (int ind = ProcId * partialSize; counter < partialSize; ++ind, ++counter) {
		int mask = 1 << (bitsCount - 1);
		int revInd = 0;
		for (int i = 0; i < bitsCount; ++i) {
			bool val = ind & mask;
			revInd |= val << i;
			mask = mask >> 1;
		}
		partialOutputSignal[counter] = inputSignal[revInd];
	}

	MPI_Gather(partialOutputSignal, partialSize, MPI_C_DOUBLE_COMPLEX, outputSignal, partialSize, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
	MPI_Bcast(outputSignal, size, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
}

__inline void Butterfly(complex<double>* signal, complex<double> u, int offset, int butterflySize, int i, complex<double>* partialSignal, int size, int j) {
	complex<double> tem = signal[offset + butterflySize] * u;
	if (size / butterflySize / 2 < NProc) {
		signal[offset + butterflySize] = signal[offset] - tem;
		signal[offset] += tem;
	}
	else {
		partialSignal[i] = signal[offset] - tem;
		partialSignal[i] = signal[offset] + tem;
	}
}

void SerialFFTCalculation(complex<double>* signal, int size) {
	int m = 0;
	for (int tmp_size = size; tmp_size > 1; tmp_size /= 2, m++);
	for (int p = 0; p < m; p++) {
		int butterflyOffset = 1 << (p + 1);
		int butterflySize = butterflyOffset >> 1;
		double coeff = PI / butterflySize;
		for (int i = 0; i < size / butterflyOffset; i++)
			for (int j = 0; j < butterflySize; j++)
				Butterfly(signal, complex<double>(cos(-j * coeff),
					sin(-j * coeff)), j + i * butterflyOffset, butterflySize);
	}
}

// FFT computation
void SerialFFT(complex<double>* inputSignal, complex<double>* outputSignal, int size) {
	BitReversing(inputSignal, outputSignal, size);
	SerialFFTCalculation(outputSignal, size);
}

void ParallelFFTCalculation(complex<double>* signal, int size) {
	int m = 0; // m = log2(size)
	for (int tmp_size = size; tmp_size > 1; tmp_size /= 2, ++m);

	for (int p = 0; p < m; ++p) {
		int butterflyOffset = 1 << (p + 1); // butterflyOffset = 2^(p + 1)
		int butterflySize = butterflyOffset >> 1; // butterflySize = 2^p = butterflyOffset / 2
		double coeff = PI / butterflySize;
		complex<double>* partialSignal = new complex<double>[size / NProc];

		for (int i = 0; i < ((size / butterflyOffset < NProc) ? (size / butterflyOffset) : (size / butterflyOffset / NProc)); ++i) {
			int ind = ((size / butterflyOffset < NProc) ? i : (i + size / butterflyOffset / NProc * ProcId));
			for (int j = 0; j < butterflySize; ++j)
				Butterfly(signal, complex<double>(cos(-j * coeff),
					sin(-j * coeff)), j + ind * butterflyOffset, butterflySize, i, partialSignal, size, j);
		}

		if (size / butterflyOffset >= NProc) {
			MPI_Gather(signal + size / NProc * ProcId, size / NProc, MPI_C_DOUBLE_COMPLEX, signal, size / NProc, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
			MPI_Bcast(signal, size, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
		}
	}
}

void ParallelFFT(complex<double>* inputSignal, complex<double>* outputSignal, int size) {
	BitReversing(inputSignal, outputSignal, size);
	//ParallelBitReversing(inputSignal, outputSignal, size);
	ParallelFFTCalculation(outputSignal, size);
}

void PrintSignal(complex<double>* signal, int size) {
	cout << "Result signal" << endl;
	for (int i = 0; i < size; i++)
		cout << signal[i] << endl;
}

void TestResult(complex<double>* inputSignal, complex<double>* outputSignal, int size) {
	// Buffer for storing the result of serial FFT
	complex<double>* testSerialSignal;
	double Accuracy = 1.e-6; // Comparison accuracy

	// Flag, that shows whether the vectors are identical
	int equal = 0;
	testSerialSignal = new complex<double>[size];
	SerialFFT(inputSignal, testSerialSignal, size);
	for (int i = 0; i < size; i++) {
		if (abs(outputSignal[i] - testSerialSignal[i]) >= Accuracy)
			equal = 1;
	}

	if (equal == 1)
		printf("The results of serial and parallel algorithms are NOT identical. Check your code.\n");
	else
		printf("The results of serial and parallel algorithms are identical.\n");

}

double calcFun(double t) {
	double sum = 0;
	for (int i = 1; i < 1000; ++i) {
		sum += cos(i * 2 * PI * t) / i;
	}

	return sum;
}

int main() {

	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &NProc);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcId);

	complex<double>* inputSignal = NULL;
	complex<double>* outputSignal = NULL;
	int size = 0;
	const int repeatCount = 16;
	double startTime;
	double duration;
	double minDuration = DBL_MAX;

	bool parallel = false;
	bool test = true;

	cout << "Fast Fourier Transform" << endl;

	if (parallel) {
		cout << "Parallel:\n";
		ProcessInitialization(inputSignal, outputSignal, size);
		for (int i = 0; i < repeatCount; i++) {
			startTime = clock();

			// FFT computation
			ParallelFFT(inputSignal, outputSignal, size);
			duration = (clock() - startTime) / CLOCKS_PER_SEC;
			if (duration < minDuration)
				minDuration = duration;
		}
		cout << setprecision(6);
		cout << "Execution time is " << minDuration << " s. " << endl;

		// Result signal output
		//PrintSignal(outputSignal, size);

		// ========================================== Task-1 ===================================================

		/*vector<double> A(size);
		double norm = 0;
		for (int i = 0; i < size; ++i) {
			A[i] = (outputSignal[i].real() * outputSignal[i].real() * 4) + (outputSignal[i].imag() * outputSignal[i].imag() * 4);
			norm += A[i];
			A[i] = sqrt(A[i]);
		}

		norm = sqrt(norm);
		for (int i = 0; i < size; ++i) {
			A[i] /= norm;
			cout << A[i] << "\n";
		}

		cout << "\n";*/

		// ========================================== Task-2 ===================================================

		for (double t = 1.0 / size; t < 1; t += 100.0 / size) {
			cout << "t = " << t << "\n";
			cout << "Accurate value: " << -log(2 * sin(PI * t)) << "\n";
			double sum = outputSignal[0].real() / size;
			for (int i = 1; i < size / 2; ++i)
				sum += outputSignal[i].real() * 2 * cos(i * 2 * PI * t) / size + outputSignal[i].imag() * -2 * sin(i * 2 * PI * t) / size;

			cout << "Fourier value: " << calcFun(t) << "\n";
			cout << "Calculated value: " << sum << "\n";
		}
	}

	if (test) {
		ProcessInitialization(inputSignal, outputSignal, size);
		ParallelFFT(inputSignal, outputSignal, size);
		//PrintSignal(outputSignal, size);
		TestResult(inputSignal, outputSignal, size);
	}

	MPI_Finalize();

	return 0;
}