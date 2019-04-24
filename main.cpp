#include <iostream>
#include <thread>
//#include "vld.h"
#include "Matrix.h"
#include "Vector.h"
#include "LU.h"
#include "Cholesky.h"
#include "QR.h"
#include "Utility.h"
#include "GaussSolver.h"



/*************************************** *
 *             2017.07.17                *
 *****************************************/

using namespace std;

int main()
{
    Matrix<double> A(3,3);
    A(0,0) = 1.0;
    A(0,1) = 2.0;
    A(0,2) = 3.0;
    A(1,0) = 4.0;
    A(1,1) = 5.0;
    A(1,2) = 6.0;
    A(2,0) = 7.0;
    A(2,1) = 8.0;
    A(2,2) = 0.0;

    cout << "A Matrix" << endl;
    cout << A << endl;

    cout << "Condition Number" << endl;
    cout << A.Getcond() << endl;


//    Vector<double> A_Vector(3);
//    A_Vector(0) = 1;
//    A_Vector(1) = 2;
//    A_Vector(2) = 3;
//
//    Matrix<double> A_mat(3,1);
//    A_mat(0,0) = 1;
//    A_mat(1,0) = 2;
//    A_mat(2,0) = 3;
//
//    cout << "Vector Norm 1:" << endl;
//    cout << A_Vector.norm_1() << endl;
//
//    cout << "Matrix Norm 1:" << endl;
//    cout << A_mat.norm_1() << endl;


    return 0;
}


