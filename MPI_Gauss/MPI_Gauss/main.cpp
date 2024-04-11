// SerialGauss.cpp
#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <time.h>
#include <math.h>
#include "mpi.h"
#include <iostream>

int NProc, ProcId;
MPI_Status st;

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
			pMatrix[i * Size + j] = rand() / double(1000);
	}
}

void MyDataInitialization(double* pMatrix, double* pVector, int Size) {
	pMatrix[0] = 1; pMatrix[1] = 1; pMatrix[2] = 4; pMatrix[3] = 4; pMatrix[4] = 9; pVector[0] = -9;
	pMatrix[5] = 2; pMatrix[6] = 2; pMatrix[7] = 17; pMatrix[8] = 17; pMatrix[9] = 82; pVector[1] = -146;
	pMatrix[10] = 2; pMatrix[11] = 0; pMatrix[12] = 3; pMatrix[13] = -1; pMatrix[14] = 4; pVector[2] = -10;
	pMatrix[15] = 0; pMatrix[16] = 1; pMatrix[17] = 4; pMatrix[18] = 12; pMatrix[19] = 27; pVector[3] = -26;
	pMatrix[20] = 1; pMatrix[21] = 2; pMatrix[22] = 2; pMatrix[23] = 10; pMatrix[24] = 0; pVector[4] = 37;
}

// Function for memory allocation and definition of the objectselements 
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
		//MyDataInitialization(pMatrix, pVector, Size);
	}
	
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

// Finding the pivot row
int ParallelFindPivotRow(double* pPartialMatrix, int PartialSize, double* pAdditionalMatrix, int AdditionalSize, int Iter, int* pPivotIter, int Size) {
	MPI_Barrier(MPI_COMM_WORLD);
	int PivotRow = -1; // The index of the pivot row
	int MaxValue = 0; // The value of the pivot element

	for (int i = 0; i < PartialSize; ++i) {
		if ((pPivotIter[i + PartialSize * ProcId] == -1) && (fabs(pPartialMatrix[i * Size + Iter]) >= MaxValue)) {
			PivotRow = i + PartialSize * ProcId;
			MaxValue = fabs(pPartialMatrix[i * Size + Iter]);
		}
	}

	if (AdditionalSize) {
		if (pPivotIter[PartialSize * NProc + ProcId - 1] == -1 && fabs(pAdditionalMatrix[Iter]) >= MaxValue) {
			PivotRow = PartialSize * NProc + ProcId - 1;
			MaxValue = fabs(pAdditionalMatrix[Iter]);
		}
	}

	if (ProcId != 0) {
		MPI_Send(&MaxValue, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
		MPI_Send(&PivotRow, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
	}
	else {
		for (int i = 1; i < NProc; ++i) {
			int TMaxValue, TPivotRow;
			MPI_Recv(&TMaxValue, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &st);
			MPI_Recv(&TPivotRow, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &st);
			if (TMaxValue > MaxValue) {
				MaxValue = TMaxValue;
				PivotRow = TPivotRow;
			}
		}
	}

	MPI_Bcast(&PivotRow, 1, MPI_INT, 0, MPI_COMM_WORLD);

	return PivotRow;
}

// Column elimination
void ParallelColumnElimination(double* pPartialMatrix, double* pPartialVector, double* pAdditionalMatrix, double* pAdditionalVector, int Pivot, int Iter, int Size, int PartialSize, int AdditionalSize, int* pPivotIter) {
	MPI_Barrier(MPI_COMM_WORLD);
	double PivotValue, PivotFactor;
	double* PivotRow = new double[Size - Iter];

	if (Pivot <= NProc * PartialSize - 1) {
		if (Pivot >= ProcId * PartialSize && Pivot < (ProcId + 1) * PartialSize) {
			PivotValue = pPartialVector[(Pivot - ProcId * PartialSize)];
			for (int i = 0; i < Size - Iter; ++i)
				PivotRow[i] = pPartialMatrix[(Pivot - ProcId * PartialSize) * Size + Iter + i];
			for (int i = 0; i < NProc; ++i)
				if (i != ProcId) {
					MPI_Send(PivotRow, Size - Iter, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
					MPI_Send(&PivotValue, 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
				}
		}
		else {
			MPI_Recv(PivotRow, Size - Iter, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &st);
			MPI_Recv(&PivotValue, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &st);
		}
	}
	else
		if (Pivot == NProc * PartialSize + ProcId - 1) {
			PivotValue = pAdditionalVector[0];
			for (int i = 0; i < Size - Iter; ++i)
				PivotRow[i] = pAdditionalMatrix[Iter + i];
			for (int i = 0; i < NProc; ++i)
				if (i != ProcId) {
					MPI_Send(PivotRow, Size - Iter, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
					MPI_Send(&PivotValue, 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
				}
		}
		else {
			MPI_Recv(PivotRow, Size - Iter, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &st);
			MPI_Recv(&PivotValue, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &st);
		}

	for (int i = 0; i < PartialSize; ++i) {
		if (pPivotIter[i + PartialSize * ProcId] == -1) {
			PivotFactor = pPartialMatrix[i * Size + Iter] / PivotRow[0];
			for (int j = Iter; j < Size; ++j)
				pPartialMatrix[i * Size + j] -= PivotFactor * PivotRow[j - Iter];
			pPartialVector[i] -= PivotFactor * PivotValue;
		}
	}

	if (AdditionalSize) {
		if (pPivotIter[NProc * PartialSize + ProcId - 1] == -1) {
			PivotFactor = pAdditionalMatrix[Iter] / PivotRow[0];
			for (int j = Iter; j < Size; ++j)
				pAdditionalMatrix[j] -= PivotFactor * PivotRow[j - Iter];
			pAdditionalVector[0] -= PivotFactor * PivotValue;
		}
	}
}

// Gaussian elimination
void ParallelGaussianElimination(double* pPartialMatrix, double* pPartialVector, int PartialSize,
	double* pAdditionalMatrix, double* pAdditionalVector, int AdditionalSize, int Size, int* pPivotPos, int* pPivotIter) {
	//MPI_Barrier(MPI_COMM_WORLD);
	int PivotRow; // The number of the current pivot row
	for (int Iter = 0; Iter < Size; ++Iter) {
		// Finding the pivot row
		PivotRow = ParallelFindPivotRow(pPartialMatrix, PartialSize, pAdditionalMatrix, AdditionalSize, Iter, pPivotIter, Size);
		pPivotPos[Iter] = PivotRow;
		pPivotIter[PivotRow] = Iter;

		/*if (ProcId == 0) {
			for (int i = 0; i < Size; ++i)
				std::cout << pPivotPos[i] << " ";
			std::cout << "\n";
		}*/

		ParallelColumnElimination(pPartialMatrix, pPartialVector, pAdditionalMatrix, pAdditionalVector, PivotRow, Iter, Size, PartialSize, AdditionalSize, pPivotIter);
	}
}

// Back substution
void ParallelBackSubstitution(double* pMatrix, double* pVector, double* pResult, int Size, int* pPivotPos) {
	int RowIndex, Row;
	for (int i = Size - 1; i >= 0; --i) {
		RowIndex = pPivotPos[i];
		pResult[i] = pVector[RowIndex] / pMatrix[Size * RowIndex + i];
		for (int j = 0; j < i; ++j) {
			Row = pPivotPos[j];
			pVector[Row] -= pMatrix[Row * Size + i] * pResult[i];
			pMatrix[Row * Size + i] = 0;
		}
	}
}

// Function for the execution of Gauss algorithm
void ParallelResultCalculation(double* pMatrix, double* pVector, double* pResult, int Size,
	double* pPartialMatrix, double* pPartialVector, int PartialSize,
	double* pAdditionalMatrix, double* pAdditionalVector, int AdditionalSize, int*& pPivotPos, int*& pPivotIter) {
	MPI_Barrier(MPI_COMM_WORLD);
	// Memory allocation
	pPivotPos = new int[Size];
	pPivotIter = new int[Size];

	for (int i = 0; i < Size; pPivotIter[i] = -1, ++i);

	// Gaussian elimination
	ParallelGaussianElimination(pPartialMatrix, pPartialVector, PartialSize, pAdditionalMatrix, pAdditionalVector, AdditionalSize, Size, pPivotPos, pPivotIter);

	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Gather(pPartialMatrix, Size * PartialSize, MPI_DOUBLE, pMatrix, Size * PartialSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gather(pPartialVector, PartialSize, MPI_DOUBLE, pVector, PartialSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if (ProcId == 0) {
		for (int i = 0; i < Size % NProc; ++i) {
			MPI_Recv(pMatrix + PartialSize * NProc * Size + i * Size, Size, MPI_DOUBLE, i + 1, 0, MPI_COMM_WORLD, &st);
			MPI_Recv(pVector + PartialSize * NProc + i, 1, MPI_DOUBLE, i + 1, 0, MPI_COMM_WORLD, &st);
		}
	}
	else if (AdditionalSize) {
		MPI_Send(pAdditionalMatrix, Size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		MPI_Send(pAdditionalVector, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	}

	/*if (ProcId == 0)
		PrintMatrix(pMatrix, Size, Size);*/

	// Back substitution
	if (ProcId == 0)
		ParallelBackSubstitution(pMatrix, pVector, pResult, Size, pPivotPos);
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

	int PartialSize;
	double* pPartialMatrix;
	double* pPartialVector;

	int AdditionalSize;
	double* pAdditionalMatrix;
	double* pAdditionalVector;

	int* pPivotPos; // The Number of pivot rows selected at the iterations
	int* pPivotIter; // The Iterations, at which the rows were pivots

	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &NProc);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcId);

	if (ProcId == 0)
		printf("Parallel Gauss algorithm for solving linear systems\n");

	MPI_Barrier(MPI_COMM_WORLD);
	// Memory allocation and definition of objects' elements
	ProcessInitialization(pMatrix, pVector, pResult, 
		pPartialMatrix, pPartialVector, Size, PartialSize, 
		AdditionalSize, pAdditionalMatrix, pAdditionalVector);

	//if (ProcId == 0) {
	//	// The matrix and the vector output
	//	printf("Initial Matrix \n");
	//	PrintMatrix(pMatrix, Size, Size);
	//	printf("Initial Vector \n");
	//	PrintVector(pVector, Size);
	//}

	// Execution of Gauss algorithm
	start = clock();
	ParallelResultCalculation(pMatrix, pVector, pResult, Size,
		pPartialMatrix, pPartialVector, PartialSize,
		pAdditionalMatrix, pAdditionalVector, AdditionalSize, pPivotPos, pPivotIter);
	finish = clock();
	duration = (finish - start) / CLOCKS_PER_SEC;

	MPI_Barrier(MPI_COMM_WORLD);
	if (ProcId == 0) {
		// Testing the result
		TestResult(pMatrix, pVector, pResult, Size);

		//// Printing the result vector
		//printf("\nResult Vector: \n");
		//PrintVector(pResult, Size);

		// Printing the execution time of Gauss method
		printf("\nTime of execution: %f\n", duration);
	}

	MPI_Finalize();

	return 0;
}
