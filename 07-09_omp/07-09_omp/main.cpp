#include <iomanip>
#include <iostream>
#include <cmath>
#include <complex>
#include <time.h>
#include <omp.h>
#include <vector>

#define PI 3.14159265358979323846

using namespace std;

//Function for simple initialization of input signal elements
void DummyDataInitialization(complex<double>*mas, int size) {
	for (int i = 0; i < size; ++i)
		mas[i] = 0;
	mas[size - size / 4] = 1;
}

// Function for random initialization of objects' elements
void RandomDataInitialization(complex<double>* mas, int size) {
	srand(unsigned(clock()));
	for (int i = 0; i < size; ++i)
		mas[i] = complex<double>(rand() / 1000.0, rand() / 1000.0);
}

void MyDataInitialization(complex<double>* mas, int size) {
	int task = 2;
	for (int i = 1; i < size; ++i) {
		switch (task) {
			case 1: mas[0] = sin(2 * PI * 10 * 1.0 / size * 0); mas[i] = sin(2 * PI * 10 * 1.0 / size * i); break;
			case 2: mas[0] = 1; mas[i] = -log(2 * sin(PI * 1.0 / size * i)); break;
		}
	}
}

//Function for memory allocation and data initialization
void ProcessInitialization(complex<double>*& inputSignal, complex<double>*& outputSignal, int& size) {
	// Setting the size of signals
	do {
		cout << "Enter the input signal length: ";
		cin >> size;
		if (size < 4)
			cout << "Input signal length should be >= 4" << endl;
		else {
			int tmpSize = size;
			while (tmpSize != 1)
			{
					if (tmpSize % 2 != 0)
					{
						cout << "Input signal length should be powers of two" << endl;
						size = -1;
						break;
					}
				tmpSize /= 2;
			}
		}
	} while (size < 4);

	cout << "Input signal length = " << size << endl;
	inputSignal = new complex<double>[size];
	outputSignal = new complex<double>[size];

	//Initialization of input signal elements - tests
	//DummyDataInitialization(inputSignal, size);

	//Computational experiments
	//RandomDataInitialization(inputSignal, size);

	MyDataInitialization(inputSignal, size);
}

//Function for computational process temination
void ProcessTermination(complex<double>*& inputSignal, complex<double>*& outputSignal) {
	delete[] inputSignal;
	inputSignal = NULL;
	delete[] outputSignal;
	outputSignal = NULL;
}

void BitReversing(complex<double>* inputSignal, complex<double>* outputSignal, int size) {
	/*int j = 0, i = 0;
	while (i < size) {
		if (j > i) {
			outputSignal[i] = inputSignal[j];
			outputSignal[j] = inputSignal[i];
		}
		else
			if (j == i)
				outputSignal[i] = inputSignal[i];

		int m = size >> 1;
		while ((m >= 1) && (j >= m)) {
			j -= m;
			m = m >> 1;
		}

		j += m;
		i++;
	}*/
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
	//bitsCount = log2(size)
	for (int tmp_size = size; tmp_size > 1; tmp_size /= 2, bitsCount++);

	//ind - index in input array
	//revInd - correspondent to ind index in output array
	#pragma omp parallel for
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

__inline void Butterfly(complex<double>* signal, complex<double> u, int offset, int butterflySize) {
	complex<double> tem = signal[offset + butterflySize] * u;
	signal[offset + butterflySize] = signal[offset] - tem;
	signal[offset] += tem;
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
	int m = 0;
	for (int tmp_size = size; tmp_size > 1; tmp_size /= 2, m++);

	for (int p = 0; p < m; ++p) {
		int butterflyOffset = 1 << (p + 1);
		int butterflySize = butterflyOffset >> 1;
		double coeff = PI / butterflySize;

		#pragma omp parallel for schedule(dynamic)
		for (int i = 0; i < size / butterflyOffset; i++) 
			for (int j = 0; j < butterflySize; j++)
				Butterfly(signal, complex<double>(cos(-j * coeff),
					sin(-j * coeff)), j + i * butterflyOffset, butterflySize);
	}
}

void ParallelFFT(complex<double>* inputSignal, complex<double>* outputSignal, int size) {
	//BitReversing(inputSignal, outputSignal, size);
	ParallelBitReversing(inputSignal, outputSignal, size);
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
		printf("The results of serial and parallel algorithms are NOT identical. Check your code.");
	else
		printf("The results of serial and parallel algorithms are identical.");

	delete[] testSerialSignal;
}

double calcFun(double t) {
	double sum = 0;
	for (int i = 1; i < 1000; ++i) {
		sum += cos(i * 2 * PI * t) / i;
	}

	return sum;
}

int main() {
	complex<double>* inputSignal = NULL;
	complex<double>* outputSignal = NULL;
	int size = 0;
	const int repeatCount = 16;
	double startTime;
	double duration;
	double minDuration = DBL_MAX;

	bool serial = false;
	bool parallel = true;
	bool test = false;
	
	cout << "Fast Fourier Transform" << endl;

	if (serial) {
		cout << "Serial:\n";
		// Memory allocation and data initialization
		ProcessInitialization(inputSignal, outputSignal, size);
		for (int i = 0; i < repeatCount; i++) {
			startTime = clock();

			// FFT computation
			SerialFFT(inputSignal, outputSignal, size);
			duration = (clock() - startTime) / CLOCKS_PER_SEC;
			if (duration < minDuration)
				minDuration = duration;
		}
		cout << setprecision(6);
		cout << "Execution time is " << minDuration << " s. " << endl;

		// Result signal output
		//PrintSignal(outputSignal, size);

		// Computational process termination
		ProcessTermination(inputSignal, outputSignal);
		cout << "\n";
	}
	
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

		// Computational process termination
		ProcessTermination(inputSignal, outputSignal);
	}

	if (test) {
		ProcessInitialization(inputSignal, outputSignal, size);
		ParallelFFT(inputSignal, outputSignal, size);
		PrintSignal(outputSignal, size);
		TestResult(inputSignal, outputSignal, size);
		ProcessTermination(inputSignal, outputSignal);
	}

	return 0;
}