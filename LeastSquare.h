//
// Created by 杨文露 on 2019-03-07.
//

#ifndef CPP_MATRIX_LEASTSQUARE_H
#define CPP_MATRIX_LEASTSQUARE_H

#include "Matrix.h"
#include "Vector.h"
#include "Cholesky.h"


template <typename T>
class LeastSquare
{
public:
    LeastSquare(const Matrix<T>& A,const Vector<T>& b);
    const Vector<T> Solve_Regularized();


private:
    Matrix<T> A;
    Vector<T> b;

};



template <typename T>
LeastSquare<T>::LeastSquare(const Matrix<T>& A_Input,const Vector<T>& b_Input)
    :A(A_Input),b(b_Input)
{

}

template <typename T>
const Vector<T> LeastSquare<T>::Solve_Regularized()
{
    //正则化求解最小二乘法问题

    //具体理论参见《数值线性代数》 第二版 P80 P81

    Matrix<T> LeftMat(A.GetNumCol(),A.GetNumCol());
    Vector<T> RightVec(A.GetNumCol(),1);

    LeftMat = A.TransPose() * A;
    RightVec = A.TransPose() * b;

    Cholesky<T> chol(LeftMat);
    Matrix<T> L = chol.LDeCompose();

}


#endif //CPP_MATRIX_LEASTSQUARE_H
