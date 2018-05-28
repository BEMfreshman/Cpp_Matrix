#include <iostream>
//#include "vld.h"
#include "Matrix.h"
#include "lu.h"
#include "transform.h"


/*************************************** *
 *             2017.07.17                *
 *****************************************/

using namespace std;

int main()
{
	int Row = 3;
	int Col = 3;

	double **arr = new double* [Row];
	for (int i = 0; i < Row; i++)
	{
		arr[i] = new double[Col];
	}
	double tmp[3][3] = { 1,4,7,2,5,8,3,6,10 };
	for (int i = 0; i < Row; i++)
	{
		for (int j = 0; j < Col; j++)
		{
			arr[i][j] = tmp[i][j];
		}
	}

	Matrix<double> A(Row, Col, arr);

	cout << A << endl;

	cout << endl;

	LU<double> lu(A);

	Matrix<double> L;
	Matrix<double> U;

	
	
	/*
	lu.DefaultFact();
	if (lu.GetL(L) == 0)
	{
		return 0;
	}

	if (lu.GetU(U) == 0)
	{
		return 0;
	}


	cout << "L" << endl;
	cout << L << endl;

	cout << "U" << endl;
	cout << U << endl;

	cout << "L*U" << endl;
	cout << L*U << endl;
	*/


	Matrix<double> P;
	Matrix<double> Q;
	lu.FullPivotFact();

	lu.GetP(P);
	lu.GetQ(Q);
	lu.GetL(L);
	lu.GetU(U);

	cout << "P" << endl;
	cout << P << endl;

	cout << "Q" << endl;
	cout << Q << endl;

	cout << "A" << endl;
	cout << A << endl;


	cout << "PAQ" << endl;
	cout << P*A*Q << endl;

	cout << "L" << endl;
	cout << L << endl;

	cout << "U" << endl;
	cout << U << endl;

	cout << "LU" << endl;
	cout << L*U << endl;

	double tmp1[3][3] = { 10,3,6,7,1,4,8,2,5 };
	for (int i = 0; i < Row; i++)
	{
		for (int j = 0; j < Col; j++)
		{
			arr[i][j] = tmp1[i][j];
		}
	}

	for (int i = 0; i < Row; i++)
	{
		delete[] arr[i];
	}
	delete[] arr;

    return 0;
}
