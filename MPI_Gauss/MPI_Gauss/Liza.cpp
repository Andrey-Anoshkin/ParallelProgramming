//
//#include <iostream>
//#include <stdio.h>
//#include <stdlib.h>
//#include <conio.h>
//#include <time.h>
//#include <math.h>
//#include "mpi.h"
//
//using namespace std;
//
//MPI_Status st;
//int NProc, ProcId; int size1, rank1;
//
//typedef struct {
//	int PivotRow;
//	double MaxValue;
//} TThreadPivotRow;
//
//void PrintMatrix(double* pMatrix, int RowCount, int ColCount) {
//	int i, j; // Loop variables
//	for (i = 0; i < RowCount; i++) {
//		for (j = 0; j < ColCount; j++)
//			cout << pMatrix[i * RowCount + j] << " ";
//		cout << "\n";
//	}
//}
//
//template <typename T>
//void PrintVector(T* pVector, int Size) {
//	int i;
//	for (i = 0; i < Size; i++)
//		cout << pVector[i] << " ";
//}
//
//void DummyDataInitialization(double* pMatrix, double* pVector, int* pPivotIter, int Size) {
//	int i, j; // Loop variables
//	for (i = 0; i < Size; i++) {
//		pPivotIter[i] = -1;
//		pVector[i] = i + 1.0;
//		for (j = 0; j < Size; j++) {
//			if (j <= i)
//				pMatrix[i * Size + j] = 1;
//			else
//				pMatrix[i * Size + j] = 0;
//		}
//	}
//}
//
//void RandomDataInitialization(double* pM, double* pV, int* pPivotIter, int Size) {
//	srand(unsigned(clock()));
//	for (int i = 0; i < Size; i++) {
//		pPivotIter[i] = -1;
//		pV[i] = rand() / double(1000);
//		for (int j = 0; j < Size; j++)
//			pM[i * Size + j] = rand() / double(1000);
//	}
//}
//
//void TaskDataInitialization(double* pMatrix, double* pVector, int* pPivotIter, int& Size) {
//	Size = 5;
//	double left[5][5]{
//		{1, 2, 3, 4, 5},
//		{2, 3, 7, 10, 13},
//		{3, 5, 11, 16, 21},
//		{2, -7, 7, 7, 2},
//		{1, 4, 5, 3, 10}
//	};
//	double right[5]{ 2, 12, 17, 57, 7 };
//	for (int i = 0; i < Size; i++) {
//		pVector[i] = right[i];
//		pPivotIter[i] = -1;
//		for (int j = 0; j < Size; j++)
//			pMatrix[i * Size + j] = left[i][j];
//	}
//}
//
//void ProcessInitialization(double*& pMatrix, double*& pVector, double*& pResult, int*& pPivotPos, int*& pPivotIter,
//	int*& PivotRowR, double*& MaxValueR, int& Size) {
//	Size = 1000;
//
//	pMatrix = new double[Size * Size];
//	pVector = new double[Size];
//	pResult = new double[Size];
//
//	pPivotPos = new int[Size];
//	pPivotIter = new int[Size];
//
//	PivotRowR = new int[Size];
//	MaxValueR = new double[Size];
//
//
//	if (rank1 == 0) {
//		// Memory allocation
//		// Initialization of the matrix and the vector elements
//		//DummyDataInitialization(pMatrix, pVector, pPivotIter, Size);
//		RandomDataInitialization(pMatrix, pVector, pPivotIter, Size);
//		//TaskDataInitialization(pMatrix, pVector, pPivotIter, Size);
//	}
//	MPI_Barrier(MPI_COMM_WORLD);
//	MPI_Bcast(pMatrix, Size * Size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//	MPI_Bcast(pVector, Size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//	MPI_Bcast(pPivotIter, Size, MPI_INT, 0, MPI_COMM_WORLD);
//
//
//	MPI_Barrier(MPI_COMM_WORLD);
//}
//
//int find_max(double* maxValues, int* pivotRows, int size) {
//	TThreadPivotRow t;
//	t.MaxValue = 0;
//	t.PivotRow = -1;
//	for (int i = 0; i < NProc; i++) {
//		if (fabs(maxValues[i]) >= t.MaxValue) {
//			t.MaxValue = fabs(maxValues[i]);
//			t.PivotRow = pivotRows[i];
//		}
//	}
//	return t.PivotRow;
//}
//
//void ParallelFindPivotRow(double* pMatrix, int* pPivotIter,
//	int* PivotRowR, double* MaxValueR, int Size, int Iter, int& res) {
//	MPI_Barrier(MPI_COMM_WORLD);
//	int t = -1;
//	// Choose the row, that stores the maximum element
//	int PivotRowC = -1;
//	double MaxValueC = 0;
//	MPI_Barrier(MPI_COMM_WORLD);
//
//	for (int i = rank1; i < Size; i += NProc) {
//		if (fabs(pMatrix[i * Size + Iter]) >= MaxValueC && pPivotIter[i] == -1) {
//			PivotRowC = i;
//			MaxValueC = fabs(pMatrix[i * Size + Iter]);
//		}
//	}
//	MPI_Gather(&MaxValueC, 1, MPI_DOUBLE, MaxValueR, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//	MPI_Gather(&PivotRowC, 1, MPI_INT, PivotRowR, 1, MPI_INT, 0, MPI_COMM_WORLD);
//	if (rank1 == 0) {
//		t = find_max(MaxValueR, PivotRowR, Size);
//		//cout << "max: " << t << "\n";
//	}
//	MPI_Barrier(MPI_COMM_WORLD);
//	res = t;
//}
//
//void ParallelColumnElimination(double* pMatrix, double* pVector, int* pPivotIter, int Pivot, int Iter, int Size) {
//	MPI_Barrier(MPI_COMM_WORLD);
//	double PivotValue, PivotFactor;
//	PivotValue = pMatrix[Pivot * Size + Iter];
//	int num = Size / NProc;
//	double* r = new double[num * Size];
//	MPI_Barrier(MPI_COMM_WORLD);
//
//	int left = Size - num * NProc;
//	double* pV_i = new double[num];
//
//	for (int i = 0; i < num; ++i) {
//		MPI_Scatter(pMatrix + i * NProc * Size, Size, MPI_DOUBLE, r + i * Size, Size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//		MPI_Scatter(pVector + i * NProc, 1, MPI_DOUBLE, pV_i + i, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//	}
//
//	/*MPI_Barrier(MPI_COMM_WORLD);
//	if (rank1 == 1) {
//		cout << "new pv_i\n"; PrintVector(pV_i, num); cout << "\n";
//	}
//	MPI_Barrier(MPI_COMM_WORLD);*/
//
//	double* r1 = new double[Size];
//
//	int it = Size - left - 1;
//	double pt;
//	MPI_Barrier(MPI_COMM_WORLD);
//
//	if (left) {
//		if (rank1 == 0) {
//			for (int i = 0; i < left; i++) {
//				MPI_Send(pMatrix + num * Size * NProc + i * Size, Size, MPI_DOUBLE, i + 1, 0, MPI_COMM_WORLD);
//				MPI_Send(pVector + num * NProc + i, 1, MPI_DOUBLE, i + 1, 0, MPI_COMM_WORLD);
//			}
//		}
//		else if (rank1 <= left) {
//			MPI_Recv(r1, Size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &st);
//			MPI_Recv(&pt, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &st);
//
//			it += rank1;
//			if (pPivotIter[it] == -1) {
//				PivotFactor = r1[Iter] / PivotValue;
//				for (int j = Iter; j < Size; j++) {
//					r1[j] -= PivotFactor * pMatrix[Pivot * Size + j];
//				}
//				pt -= PivotFactor * pVector[Pivot];
//			}
//		}
//		/*MPI_Barrier(MPI_COMM_WORLD);
//		if (rank1 == 1) {
//			cout << "new left matrix\n"; PrintVector(r1, Size); cout << "\n";
//		}
//		MPI_Barrier(MPI_COMM_WORLD);*/
//	}
//	MPI_Barrier(MPI_COMM_WORLD);
//
//
//	for (int i = 0; i < num; i++) {
//		if (pPivotIter[i * NProc + rank1] == -1) {
//			PivotFactor = r[i * Size + Iter] / PivotValue;
//			for (int j = Iter; j < Size; j++) {
//				r[i * Size + j] -= PivotFactor * pMatrix[Pivot * Size + j];
//			}
//			pV_i[i] -= PivotFactor * pVector[Pivot];
//		}
//	}
//	MPI_Barrier(MPI_COMM_WORLD);
//
//	/*MPI_Barrier(MPI_COMM_WORLD);
//	if (rank1 == 1) { cout << "new pV_i\n"; PrintVector(pV_i, num); cout << "\n"; }
//	MPI_Barrier(MPI_COMM_WORLD);*/
//
//
//	for (int i = 0; i < num; i++) {
//		MPI_Gather(pV_i + i, 1, MPI_DOUBLE, pVector + i * NProc, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//		MPI_Gather(r + i * Size, Size, MPI_DOUBLE, pMatrix + i * Size * NProc, Size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//		//MPI_Barrier(MPI_COMM_WORLD);
//	}
//	if (left) {
//		if (rank1 == 0) {
//			for (int i = 0; i < left; ++i) {
//				MPI_Recv(pVector + num * Size + i, 1, MPI_DOUBLE, i + 1, 0, MPI_COMM_WORLD, &st);
//				MPI_Recv(pMatrix + num * Size * NProc + i * Size, Size, MPI_DOUBLE, i + 1, 0, MPI_COMM_WORLD, &st);
//			}
//		}
//		else if (rank1 <= left) {
//			MPI_Send(&pt, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
//			MPI_Send(r1, Size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
//		}
//	}
//
//	MPI_Barrier(MPI_COMM_WORLD);
//
//}
//
//void ParallelGaussianElimination(double* pMatrix, double* pVector, int* pPivotPos, int* pPivotIter,
//	int* PivotRowR, double* MaxValueR, int Size) {
//	int PivotRow; // The number of the current pivot row
//
//	//if (rank1 == 0) PrintVector(pPivotIter, Size); cout << "\n";
//	MPI_Barrier(MPI_COMM_WORLD);
//
//	for (int Iter = 0; Iter < Size; Iter++) {
//		// Finding the pivot row
//		ParallelFindPivotRow(pMatrix, pPivotIter, PivotRowR, MaxValueR, Size, Iter, PivotRow);
//		MPI_Bcast(&PivotRow, 1, MPI_INT, 0, MPI_COMM_WORLD);
//		pPivotPos[Iter] = PivotRow;
//		pPivotIter[PivotRow] = Iter;
//
//		/*MPI_Barrier(MPI_COMM_WORLD);
//		if (rank1 == 0) {
//			cout << "new pivot\n"; PrintVector(pPivotIter, Size); cout << "\n";
//		}
//		MPI_Barrier(MPI_COMM_WORLD);*/
//
//		ParallelColumnElimination(pMatrix, pVector, pPivotIter, PivotRow, Iter, Size);
//	}
//	MPI_Barrier(MPI_COMM_WORLD);
//}
//
//void ParallelBackSubstitution(double* pMatrix, double* pVector, int* pPivotPos, double* pResult, int Size) {
//	int RowIndex, Row;
//	for (int i = Size - 1; i >= 0; i--) {
//		RowIndex = pPivotPos[i];
//		pResult[i] = pVector[RowIndex] / pMatrix[Size * RowIndex + i];
//		for (int j = 0; j < i; j++) {
//			Row = pPivotPos[j];
//			pVector[Row] -= pMatrix[Row * Size + i] * pResult[i];
//			pMatrix[Row * Size + i] = 0;
//		}
//	}
//}
//
//void ParallelResultCalculation(double* pMatrix, double* pVector, double* pResult, int* pPivotPos, int* pPivotIter,
//	int* PivotRowR, double* MaxValueR, int Size) {
//
//	ParallelGaussianElimination(pMatrix, pVector, pPivotPos, pPivotIter, PivotRowR, MaxValueR, Size);
//
//	if (rank1 == 0)
//		ParallelBackSubstitution(pMatrix, pVector, pPivotPos, pResult, Size);
//}
//
//void TestResult(double* pMatrix, double* pVector, double* pResult, int Size) {
//	/* Buffer for storing the vector, that is a result of multiplication
//	of the linear system matrix by the vector of unknowns */
//	double* pRightPartVector;
//	// Flag, that shows wheather the right parts
//	// vectors are identical or not
//	int equal = 0;
//	double Accuracy = 1.e-3; // Comparison accuracy
//	pRightPartVector = new double[Size];
//	for (int i = 0; i < Size; i++) {
//		pRightPartVector[i] = 0;
//		for (int j = 0; j < Size; j++) {
//			pRightPartVector[i] +=
//				pMatrix[i * Size + j] * pResult[j];
//		}
//	}
//	for (int i = 0; i < Size; i++) {
//		if (fabs(pRightPartVector[i] - pVector[i]) > Accuracy)
//			equal = 1;
//	}
//	if (equal == 1)
//		printf("The result of the parallel Gauss algorithm is NOT correct."
//			"Check your code.");
//	else
//		printf("The result of the parallel Gauss algorithm is correct.");
//	//delete[] pRightPartVector;
//}
//
//int main(int argc, char* argv[]) {
//	double* pMatrix; double* pVector; double* pResult; double* MaxValueR;
//	int* pPivotPos; int* pPivotIter; int* PivotRowR;
//	int Size; // The size of the matrix and the vectors
//
//	MPI_Init(NULL, NULL);
//	MPI_Comm_size(MPI_COMM_WORLD, &NProc);
//	MPI_Comm_rank(MPI_COMM_WORLD, &rank1);
//	double tm = MPI_Wtime();
//	MPI_Barrier(MPI_COMM_WORLD);
//
//	// Data initialization
//	ProcessInitialization(pMatrix, pVector, pResult, pPivotPos, pPivotIter, PivotRowR, MaxValueR, Size);
//
//	MPI_Barrier(MPI_COMM_WORLD);
//
//
//
//	double start = clock();
//	ParallelResultCalculation(pMatrix, pVector, pResult, pPivotPos, pPivotIter, PivotRowR, MaxValueR, Size);
//	double finish = clock();
//	double duration = (finish - start) / CLOCKS_PER_SEC;
//	MPI_Barrier(MPI_COMM_WORLD);
//	if (rank1 == 0) cout << " Elapsed time = " << duration << endl;
//	MPI_Barrier(MPI_COMM_WORLD);
//
//
//	//if (rank1 == 0) PrintVector(pResult, Size);
//	// Testing the result
//	if (rank1 == 0) TestResult(pMatrix, pVector, pResult, Size);
//	MPI_Finalize();
//
//
//	return 0;
//}
