#include <iostream>
#include <thread>
//#include "vld.h"
#include "Matrix.h"
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

    cout << "Condition Number" << endl;
    cout << A.Getcond() << endl;

    return 0;
}


