#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <time.h>
#include <math.h>
#include <omp.h>
#include <algorithm>

// Function for simple initialization of the matrix
// and the vector elements
void DummyDataInitialization(double* pMatrix, double* pVector, int Size) {
	for (int i = 0; i < Size; ++i) {
		pVector[i] = i + 1.0;
		for (int j = 0; j < Size; ++j)
			if (j <= i)
				pMatrix[i * Size + j] = 1;
			else
				pMatrix[i * Size + j] = 0;
	}
}

// Function for random initialization of the matrix
// and the vector elements
void RandomDataInitialization(double* pMatrix, double* pVector, int Size) {
	srand(unsigned(clock()));
	for (int i = 0; i < Size; ++i) {
		pVector[i] = rand() / double(1000);
		for (int j = 0; j < Size; ++j)
			if (j <= i)
				pMatrix[i * Size + j] = rand() / double(1000);
			else
				pMatrix[i * Size + j] = 0;
	}
}

// Function for memory allocation and definition of the objectselements 
void ProcessInitialization(double*& pMatrix, double*& pVector, double*& pResult, int& Size) {
	// Setting the size of the matrix and the vector
	do {
		printf("\nEnter size of the matrix and the vector: ");
		scanf_s("%d", &Size);
		printf("\nChosen size = %d \n", Size);
		if (Size <= 0)
			printf("\nSize of objects must be greater than 0!\n");
	} while (Size <= 0);


	// Memory allocation
	pMatrix = new double[Size * Size];
	pVector = new double[Size];
	pResult = new double[Size];

	// Initialization of the matrix and the vector elements
	DummyDataInitialization(pMatrix, pVector, Size);
	//RandomDataInitialization(pMatrix, pVector, Size);
}

// Function for formatted matrix output
void PrintMatrix(double* pMatrix, int RowCount, int ColCount) {
	for (int i = 0; i < RowCount; ++i) {
		for (int j = 0; j < ColCount; ++j)
			printf("%7.4f ", pMatrix[i * RowCount + j]);
		printf("\n");
	}
}

// Function for formatted vector output
void PrintVector(double* pVector, int Size) {
	for (int i = 0; i < Size; ++i)
		printf("%7.4f ", pVector[i]);
}

// Function for the execution of Gauss algorithm
void ParallelResultCalculation(double* pMatrix, double* pVector, double* pResult, int Size) {

	for (int i = 0; i < Size; pResult[i] = pVector[i] / pMatrix[i * Size + i], ++i);
	double maxDif = 1e+6;
	double w = 1.5;

	for (; maxDif > 1e-9;) {
		maxDif = 0;
		double* pNew = new double[Size];
		#pragma omp parallel for
		for (int i = 0; i < Size; ++i) {
			double sum = 0;
			for (int j = i + 1; j < Size; ++j)
				sum += pMatrix[i * Size + j] * pResult[j];
			
			pNew[i] = (1 - w) * pResult[i] + w * (pVector[i] - sum) / pMatrix[i * Size + i];
		}

		for (int i = 0; i < Size; ++i) {
			double sum = 0;
			for (int j = 0; j < i; ++j)
				sum += pMatrix[i * Size + j] * pNew[j];
			pNew[i] += w * -sum / pMatrix[i * Size + i];
		}

		for (int i = 0; i < Size; ++i) {
			maxDif = std::max(maxDif, fabs(pNew[i] - pResult[i]));
			pResult[i] = pNew[i];
		}

		delete[] pNew;
	}
}

// Function for computational process termination
void ProcessTermination(double* pMatrix, double* pVector, double* pResult) {
	delete[] pMatrix;
	delete[] pVector;
	delete[] pResult;
}

// Function for testing the result
void TestResult(double* pMatrix, double* pVector, double* pResult, int Size) {
	/* Buffer for storing the vector, that is a result of multiplication
	of the linear system matrix by the vector of unknowns */
	double* pRightPartVector;

	// Flag, that shows wheather the right parts vectors are identical or not
	int equal = 0;
	double Accuracy = 1e-3; // Comparison accuracy
	pRightPartVector = new double[Size];
	for (int i = 0; i < Size; ++i) {
		pRightPartVector[i] = 0;
		for (int j = 0; j < Size; ++j)
			pRightPartVector[i] += pMatrix[i * Size + j] * pResult[j];
	}

	for (int i = 0; i < Size; i++)
		if (fabs(pRightPartVector[i] - pVector[i]) > Accuracy)
			equal = 1;

	if (equal == 1)
		printf("\nThe result of the parallel Gauss algorithm is NOT correct. Check your code.");
	else
		printf("\nThe result of the parallel Gauss algorithm is correct.");

	delete[] pRightPartVector;
}
int main() {
	double* pMatrix; // The matrix of the linear system
	double* pVector; // The right parts of the linear system
	double* pResult; // The result vector
	int Size; // The sizes of the initial matrix and the vector
	double start, finish, duration;
	printf("Parallel upper relaxation algorithm for solving linear systems\n");

	// Memory allocation and definition of objects' elements
	ProcessInitialization(pMatrix, pVector, pResult, Size);

	// The matrix and the vector output
	/*printf("Initial Matrix \n");
	PrintMatrix(pMatrix, Size, Size);
	printf("Initial Vector \n");
	PrintVector(pVector, Size);*/

	// Execution of Gauss algorithm
	start = clock();
	ParallelResultCalculation(pMatrix, pVector, pResult, Size);
	finish = clock();
	duration = (finish - start) / CLOCKS_PER_SEC;

	// Testing the result
	TestResult(pMatrix, pVector, pResult, Size);

	// Printing the result vector
	/*printf("\nResult Vector: \n");
	PrintVector(pResult, Size);*/

	// Printing the execution time of Gauss method
	printf("\nTime of execution: %f\n", duration);

	// Computational process termination
	ProcessTermination(pMatrix, pVector, pResult);

	return 0;
}