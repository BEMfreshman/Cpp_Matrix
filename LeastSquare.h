//
// Created by 杨文露 on 2019-03-07.
//

#ifndef CPP_MATRIX_LEASTSQUARE_H
#define CPP_MATRIX_LEASTSQUARE_H

#include "Matrix.h"
#include "Vector.h"
#include "Cholesky.h"
#include "QR.h"

//理论 《数值线性代数》 第二版 p92

template <typename T>
class LeastSquare
{
public:
    LeastSquare(const Matrix<T>& A,const Vector<T>& b);
    const Vector<T> Solve_Regularized();
    const Vector<T> Solve_QR();


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

    //ToDo: 完成最小二乘法中正则化方法

    Matrix<T> y(A.GetNumCol(),1);
    Matrix<T> ans(L.GetNumRow(),1);

    if(L.isLowTri())
    {
        ans = Utility<T>::UptriSolve(L.TransPose(),Utility<T>::LowtriSolve(L,RightVec));
        return ans;
    }
    else if(L.isUpTri())
    {
        ans = Utility<T>::LowtriSolve(L.TransPose(),Utility<T>::UptriSolve(L,RightVec));
        return ans;
    }
    else
    {
        throw runtime_error("Chol fact failed");
    }

}

template <typename T>
const Vector<T> LeastSquare<T>::Solve_QR()
{
    //QR法求解最小二乘问题

    QR<T> qr(A);
    vector<Matrix<T>> QRMatrix = qr.QRHouseHolder();

    Matrix<T> Q = QRMatrix[0];
    Matrix<T> R = QRMatrix[1];

    Matrix<T> c = Q.TransPose() * b;
    Matrix<T> c1 = c.ExtractBlock(0,0,R.GetNumRow(),1);

    if(R.isUpTri())
    {
        Matrix<T> b = UptriSolve(R,c1);
        return b;
    }
    else
    {
        throw runtime_error("QR Fact failed");
    }

}


#endif //CPP_MATRIX_LEASTSQUARE_H
