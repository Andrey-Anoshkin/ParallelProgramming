#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <time.h>
#include <math.h>
#include <algorithm>
#include "mpi.h"

int NProc, ProcId;
MPI_Status st;

void PrintMatrix(double* pMatrix, int RowCount, int ColCount) {
	for (int i = 0; i < RowCount; ++i) {
		for (int j = 0; j < ColCount; ++j)
			printf("%7.4f ", pMatrix[i * RowCount + j]);
		printf("\n");
	}
}

void PrintVector(double* pVector, int Size) {
	for (int i = 0; i < Size; ++i)
		printf("%7.4f ", pVector[i]);
}

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

void ProcessInitialization(double*& pMatrix, double*& pVector, double*& pResult,
	double*& pPartialMatrix, double*& pPartialVector, int& Size, int& PartialSize,
	int& AdditionalSize, double*& pAdditionalMatrix, double*& pAdditionalVector) {
	MPI_Barrier(MPI_COMM_WORLD);

	Size = 3000;

	// Memory allocation
	pMatrix = new double[Size * Size];
	pVector = new double[Size];
	pResult = new double[Size];

	PartialSize = Size / NProc;

	// Partial arrays memory allocation
	pPartialMatrix = new double[Size * PartialSize];
	pPartialVector = new double[PartialSize];

	AdditionalSize = ((ProcId > 0) && (ProcId <= (Size % NProc)) && ((Size % NProc) > 0));

	// Additional arrays memory allocation
	pAdditionalMatrix = new double[Size * AdditionalSize];
	pAdditionalVector = new double[AdditionalSize];

	if (ProcId == 0) {
		// Initialization of the matrix and the vector elements
		//DummyDataInitialization(pMatrix, pVector, Size);
		RandomDataInitialization(pMatrix, pVector, Size);
		//PrintMatrix(pMatrix, Size, Size); printf("\n");
		for (int i = 0; i < Size; pResult[i] = pVector[i] / pMatrix[i * Size + i], ++i);
	}

	MPI_Bcast(pResult, Size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	MPI_Scatter(pMatrix, Size * PartialSize, MPI_DOUBLE, pPartialMatrix, Size * PartialSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Scatter(pVector, PartialSize, MPI_DOUBLE, pPartialVector, PartialSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if (Size % NProc) {
		if (ProcId == 0)
			for (int i = 0; i < Size % NProc; ++i) {
				MPI_Send(pMatrix + Size * (PartialSize * NProc + i), Size, MPI_DOUBLE, i + 1, 0, MPI_COMM_WORLD);
				MPI_Send(pVector + (PartialSize * NProc + i), 1, MPI_DOUBLE, i + 1, 0, MPI_COMM_WORLD);
			}
		else if (AdditionalSize) {
			MPI_Recv(pAdditionalMatrix, Size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &st);
			MPI_Recv(pAdditionalVector, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &st);
		}
	}
}

void ParallelResultCalculation(double* pPartialMatrix, double* pPartialVector, double* pResult, int PartialSize, 
	double* pAdditionalMatrix, double* pAdditionalVector, int AdditionalSize, int Size) {
	double maxDif = 1e+6;

	for (; maxDif > 1e-9;) {
		MPI_Barrier(MPI_COMM_WORLD);
		maxDif = 1;
		double mDif = 0;
		double* pNew = new double[PartialSize];
		double pNewAdditional;

		for (int i = 0; i < PartialSize; ++i) {
			pNew[i] = pPartialVector[i];
			for (int j = 0; j < Size; ++j)
				if (PartialSize * ProcId + i != j)
					pNew[i] -= pPartialMatrix[i * Size + j] * pResult[j];

			pNew[i] /= pPartialMatrix[i * Size + PartialSize * ProcId + i];
		}

		if (AdditionalSize) {
			pNewAdditional = pAdditionalVector[0];
			for (int j = 0; j < Size; ++j)
				if (j != PartialSize * NProc + ProcId - 1)
					pNewAdditional -= pAdditionalMatrix[j] * pResult[j];

			pNewAdditional /= pAdditionalMatrix[PartialSize * NProc + ProcId - 1];
		}

		for (int i = 0; i < PartialSize; ++i)
			mDif = std::max(mDif, fabs(pNew[i] - pResult[ProcId * PartialSize + i]));
			
		if (AdditionalSize)
			mDif = std::max(mDif, fabs(pNewAdditional - pResult[PartialSize * NProc + ProcId - 1]));

		/*if (ProcId == 0) {
			PrintVector(pResult, Size); 
			printf("\n");
		}*/

		MPI_Gather(pNew, PartialSize, MPI_DOUBLE, pResult, PartialSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		if (ProcId == 0) 
			for (int i = 0; i < Size % NProc; ++i)
				MPI_Recv(pResult + PartialSize * NProc + i, 1, MPI_DOUBLE, i + 1, 0, MPI_COMM_WORLD, &st);
		else if (AdditionalSize) 
			MPI_Send(&pNewAdditional, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		
		MPI_Bcast(pResult, Size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		MPI_Reduce(&mDif, &maxDif, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
		MPI_Bcast(&maxDif, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}
}

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
		printf("\nThe result of the parallel Gauss algorithm is NOT correct. Check your code.\n");
	else
		printf("\nThe result of the parallel Gauss algorithm is correct.\n");

	delete[] pRightPartVector;
}
int main() {
	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &NProc);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcId);

	double* pMatrix; // The matrix of the linear system
	double* pVector; // The right parts of the linear system
	double* pResult; // The result vector
	int Size; // The sizes of the initial matrix and the vector
	double start, finish, duration;

	if (ProcId == 0)
		printf("Parallel Yakobi algorithm for solving linear systems\n");

	int PartialSize;
	double* pPartialMatrix;
	double* pPartialVector;

	int AdditionalSize;
	double* pAdditionalMatrix;
	double* pAdditionalVector;

	// Memory allocation and definition of objects' elements
	ProcessInitialization(pMatrix, pVector, pResult,
		pPartialMatrix, pPartialVector, Size, PartialSize,
		AdditionalSize, pAdditionalMatrix, pAdditionalVector);

	

	// Execution of Gauss algorithm
	start = clock();
	ParallelResultCalculation(pPartialMatrix, pPartialVector, pResult, PartialSize,
		pAdditionalMatrix, pAdditionalVector, AdditionalSize, Size);
	finish = clock();
	duration = (finish - start) / CLOCKS_PER_SEC;

	if (ProcId == 0)
		TestResult(pMatrix, pVector, pResult, Size);

	/*if (ProcId == 0)
		PrintVector(pResult, Size);*/

	// Printing the result vector
	/*printf("\nResult Vector: \n");
	PrintVector(pResult, Size);*/

	if (ProcId == 0)
		printf("\nTime of execution: %f\n", duration);

	MPI_Finalize();

	return 0;
}