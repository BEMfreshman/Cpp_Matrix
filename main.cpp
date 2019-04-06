#include <iostream>
#include <thread>
//#include "vld.h"
#include "Matrix.h"
#include "lu.h"
#include "Cholesky.h"
#include "QR.h"
#include "utility.h"
#include "GaussSolver.h"



/*************************************** *
 *             2017.07.17                *
 *****************************************/

using namespace std;

int main()
{
    size_t Row = 3;
    size_t Col = 3;

    double **arr = new double* [Row];
    for (size_t i = 0; i < Row; i++)
    {
        arr[i] = new double[Col];
    }
    double tmp[3][3] = { 0,3,1,0,4,-2,2,1,1};
    for (size_t i = 0; i < Row; i++)
    {
        for (size_t j = 0; j < Col; j++)
        {
            arr[i][j] = tmp[i][j];
        }
    }

    Matrix<double> A(Row, Col, arr);

    cout << "A" << endl;
    cout << A << endl;


    for (size_t i = 0; i < Row; i++)
    {
        delete[] arr[i];
    }
    delete[] arr;


    QR<double> qr(A);
    vector<Matrix<double>> vecMat= qr.QRHouseHolder();

    cout << "Q" << endl;
    cout << vecMat[0] << endl;

    cout << "R" << endl;
    cout << vecMat[1] << endl;


    cout << " Product" << endl;
    cout << vecMat[0]* vecMat[1] << endl;

    return 0;
}


