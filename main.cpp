#include <iostream>
//#include "vld.h"
#include "Matrix.h"
#include "lu.h"
#include "Cholesky.h"
#include "transform.h"


/*************************************** *
 *             2017.07.17                *
 *****************************************/

using namespace std;

int main()
{
	int Row = 4;
	int Col = 4;

	double **arr = new double* [Row];
	for (int i = 0; i < Row; i++)
	{
		arr[i] = new double[Col];
	}
	double tmp[4][4] = { 4,-2,4,2,-2,10,-2,-7,4,-2,8,4,2,-7,4,7 };
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

	/*LU<double> lu(A);

	Matrix<double> L;
	Matrix<double> U;*/

	
	
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


	/*Matrix<double> P;
	lu.ColPivotFact();

	lu.GetP(P);
	lu.GetL(L);
	lu.GetU(U);

	cout << "P" << endl;
	cout << P << endl;

	

	cout << "A" << endl;
	cout << A << endl;


	cout << "PA" << endl;
	cout << P*A << endl;

	cout << "L" << endl;
	cout << L << endl;

	cout << "U" << endl;
	cout << U << endl;

	cout << "LU" << endl;
	cout << L*U << endl;*/


	
	Cholesky<double> CH(A);
	CH.DefaultFact();


	Matrix<double> L;
	L = CH.GetL();

	cout << "L" << endl;
	cout << L << endl;

	cout << "L * L'" << endl;
	cout << L*L.TransPose() << endl;



	for (int i = 0; i < Row; i++)
	{
		delete[] arr[i];
	}
	delete[] arr;

    return 0;
}
