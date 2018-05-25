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

	Matrix<double> mat(Row, Col, arr);

	cout << mat << endl;

	cout << endl;

	//double detRes = Det(mat);

	//cout << detRes << endl;

	LU<double> lu(mat);
	lu.DefaultFact();

	Matrix<double> L;
	Matrix<double> U;
	double detres;

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


	for (int i = 0; i < Row; i++)
	{
		delete[] arr[i];
	}
	delete[] arr;

    return 0;
}
