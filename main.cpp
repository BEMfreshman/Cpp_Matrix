#include <iostream>
#include "vld.h"
#include "Matrix.h"
#include "decomposition.h"


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
	double tmp[3][3] = { 2, 1, 3, 4, 2, 1, 6, -3, 4 };
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

	double detRes = Det(mat);

	cout << detRes << endl;



	for (int i = 0; i < Row; i++)
	{
		delete[] arr[i];
	}
	delete[] arr;

    return 0;
}
