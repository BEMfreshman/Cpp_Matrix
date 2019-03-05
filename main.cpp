#include <iostream>
#include <thread>
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
    size_t Row = 4;
    size_t Col = 4;

	double **arr = new double* [Row];
    for (size_t i = 0; i < Row; i++)
	{
		arr[i] = new double[Col];
	}
	double tmp[4][4] = { 4,-2,4,2,-2,10,-2,-7,4,-2,8,4,2,-7,4,7 };
    for (size_t i = 0; i < Row; i++)
	{
        for (size_t j = 0; j < Col; j++)
		{
			arr[i][j] = tmp[i][j];
		}
	}

	Matrix<double> A(Row, Col, arr);

	cout << A << endl;
	cout << endl;


	Cholesky<double> CH(A);
    Matrix<double> L = CH.LDeCompose();


	cout << "L" << endl;
	cout << L << endl;

    cout << "L * L'" << endl;
    cout << L * L.TransPose() << endl;



    for (size_t i = 0; i < Row; i++)
	{
		delete[] arr[i];
	}
	delete[] arr;

    return 0;
}
