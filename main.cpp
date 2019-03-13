#include <iostream>
#include <thread>
//#include "vld.h"
#include "Matrix.h"
#include "lu.h"
#include "Cholesky.h"
#include "transform.h"
#include "GaussSolver.h"



/*************************************** *
 *             2017.07.17                *
 *****************************************/

using namespace std;

int main()
{
    size_t Row = 5;
    size_t Col = 5;

    double **arr = new double* [Row];
    for (size_t i = 0; i < Row; i++)
    {
        arr[i] = new double[Col];
    }
    double tmp[5][5] = { 1,2,3,4,5,0,1,2,3,4,0,0,1,2,3,0,0,0,1,2,0,0,0,0,1};
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

//
//
//    Cholesky<double> CH(A);
//    Matrix<double> L = CH.LDeCompose();
//
//
//    cout << "L" << endl;
//    cout << L << endl;
//
//    cout << "L * L'" << endl;
//    cout << L * L.TransPose() << endl;
//
//
//

    Matrix<double> B(Row,1);
    B.SetZeros();

    B(0,0) = 1;
    B(1,0) = 2;
    B(2,0) = 3;
    B(3,0) = 4;
    B(4,0) = 5;

    cout << "B" << endl;
    cout << B << endl;


    GaussSolver<double> GS(A,B);

    GS.UptriSolve();

    Matrix<double> ans = GS.getAns();


    cout << "ans" << endl;
    cout << ans << endl;




    for (size_t i = 0; i < Row; i++)
    {
        delete[] arr[i];
    }
    delete[] arr;



    return 0;
}


