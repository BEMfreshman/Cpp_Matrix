//
// Created by YangWenlu on 2019/4/6.
//

#ifndef CPP_MATRIX_QR_H
#define CPP_MATRIX_QR_H

#include "Matrix.h"
#include "Vector.h"
#include "utility.h"
#include <vector>

using namespace std;

template <typename T>
class QR
{
public:
    QR(const Matrix<T>& A);

    vector<Matrix<T>> QRHouseHolder();

private:
    Matrix<T> A;

};




template <typename T>
QR<T>::QR(const Matrix<T>& AInput)
:A(AInput)
{

}

template <typename T>
vector<Matrix<T>> QR<T>::QRHouseHolder()
{
    /*
     * 按照《数值线性代数》 P94~P95页的理论和存储方案完成该算法，
     * LAPACK dgeqrf()
     */

    size_t NumRow = A.GetNumRow();
    size_t NumCol = A.GetNumCol();

    if(NumRow < NumCol)
    {
        throw runtime_error("NumRow is less than NumCol");
    }

    auto* d = new double[NumCol];


    size_t iterNum = NumRow == NumCol ? NumCol -1 :NumCol;

    for(size_t i = 0 ; i < iterNum;i++)
    {
        Vector<T> ColVec(A.ExtractBlock(i,i,NumCol-i,1));
        vector<Matrix<T>> MuBeta = ColVec.House();
        Matrix<T> MuMat = MuBeta[0];
        Matrix<T> BetaMat = MuBeta[1];


        cout << "Beta" << endl;
        cout << BetaMat << endl;

        cout << "MuMat" << endl;
        cout << MuMat << endl;

        Matrix<T> H(MuMat.GetNumRow(),MuMat.GetNumRow());
        H.IdentityMatrix();
        H -=BetaMat(0,0) * MuMat * MuMat.TransPose();

        cout << "H" << endl;
        cout << H << endl;


        Matrix<T> NewABlock = A.ExtractBlock(i,i,NumRow - i,NumCol - i);


        NewABlock  =  H * NewABlock;
        A.SetBlock(i,i,NumRow - i,NumCol - i,NewABlock);

        cout << "A after NewBlock" << endl;
        cout << A << endl;

        A.SetBlock(i+1,i,NumCol - (i + 1),1,MuMat.ExtractBlock(1,0,NumCol - i - 1,1));

        cout  << "A Last" << endl;
        cout << A << endl;

        d[i] = static_cast<double> (BetaMat(0,0));
    }

    Matrix<T> Q(NumRow,NumRow);
    Matrix<T> R(NumRow,NumCol);


    //组装R
    R.SetZeros();
    for(size_t i = 0; i < NumRow;i++)
    {
        for(size_t j = 0 ; j < NumCol;j++)
        {
            if(i <= j)
            {
                R(i,j) = A(i,j);
            }
        }
    }


    //计算Q

    Q.IdentityMatrix();
    for(size_t i = 0 ; i < iterNum;i++)
    {
        Matrix<T> MuMat(NumRow,1);
        double beta = d[i];

        MuMat.SetZeros();
        MuMat(i,0) = 1;

        Matrix<T> MuStorageInA = A.ExtractBlock(i + 1,i,NumCol - (i + 1),1);

        MuMat.SetBlock(i+1,0,NumRow - (i + 1),1,MuStorageInA);

        Matrix<T> Hi(NumRow,NumRow);
        Hi.IdentityMatrix();
        Hi -=beta * MuMat * MuMat.TransPose();

        cout << "Hi" << endl;
        cout << Hi << endl;

        Q *= Hi;
    }

    delete [] d;
    vector<Matrix<T>> RC;
    RC.push_back(Q);
    RC.push_back(R);
    return RC;
}



#endif //CPP_MATRIX_QR_H
