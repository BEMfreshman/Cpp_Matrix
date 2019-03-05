#ifndef CHOLESKY_H
#define CHOLESKY_H

#include <iostream>
#include <vector>
#include "Matrix.h"

using namespace std;

template<typename T>
class Cholesky
{
public:
	Cholesky(const Matrix<T>& A);
	~Cholesky();

    Matrix<T> LDeCompose();               //A = L * L';
    vector<Matrix<T>> LDDeCompose();      //A = L * D * L';



private:
	Matrix<T> A;
};


template<typename T>
Cholesky<T>::Cholesky(const Matrix<T>& A)
{
	this->A = A;
}

template<typename T>
Cholesky<T>::~Cholesky()
{

}

template<typename T>
Matrix<T> Cholesky<T>::LDeCompose()
{
	//A = L * L';
	int row = A.GetNumRow();

	Matrix<T> matColBelowPivot;
	Matrix<T> matRowLeftPivot;
	Matrix<T> matLeftBelowPivot;


	Matrix<T> mattmp;

    for (size_t i = 0; i < A.GetNumCol(); i++)
	{
		
		if (i != A.GetNumCol() - 1)
		{
			if (i != 0)
			{
				matRowLeftPivot = A.ExtractBlock(i, 0, 1, i);
				matColBelowPivot = A.ExtractBlock(i+1, i, row - i - 1, 1);
				matLeftBelowPivot = A.ExtractBlock(i+1, 0, row - i - 1, i);

				mattmp = matRowLeftPivot * matRowLeftPivot.TransPose();

				double lk_1p = mattmp(0, 0);  //sum(lkp*lkp)   p = 1 To k - 1

				A(i, i) -= lk_1p;

				matColBelowPivot -= matLeftBelowPivot * matRowLeftPivot.TransPose();
				A.SetBlock(i+1, i, row - i - 1, 1, matColBelowPivot);
			}

			A(i, i) = sqrt(A(i, i));

			matColBelowPivot = A.ExtractBlock(i + 1, i, row - i - 1, 1);
			matColBelowPivot /= A(i, i);

			A.SetBlock(i + 1, i, row - i - 1, 1, matColBelowPivot);


		}
		else
		{
			matRowLeftPivot = A.ExtractBlock(i, 0, 1, i);
			mattmp = matRowLeftPivot * matRowLeftPivot.TransPose();
			
			double lk_1p = mattmp(0, 0);  //sum(lkp*lkp)   p = 1 To k - 1
			A(i, i) -= lk_1p;

			A(i, i) = sqrt(A(i, i));
		}
	}

	

    Matrix<T> L(A.GetNumRow(), A.GetNumCol());
	L.SetZeros();

    for (size_t i = 0; i < A.GetNumRow(); i++)
	{
        for (size_t j = 0; j < A.GetNumCol(); j++)
		{
			if (i >= j)
			{
				L(i, j) = A(i, j);
			}
		}
	}
    return L;
}

template<typename T>
vector<Matrix<T>> Cholesky<T>::LDDeCompose()
{
    //算法参见 《数值线性代数》 第二版  P31

    int row = A.GetNumRow();
    int col = A.GetNumCol();

    Matrix<T> v;
    Matrix<T> matRowLeftPivot;
    Matrix<T> tmp;
    double tmpValue;

    for (size_t j = 0; j < col; j++)
    {
        if (j != 0)
        {
            v.Resize(j, 1);
            for (size_t i = 0; i < j; i++)
            {
                v(i, 0) = A(j, i)*A(i, i);
            }

            matRowLeftPivot = A.ExtractBlock(j, 0, 1, j);
            tmp = matRowLeftPivot * v;
            A(j,j) -= tmp(0, 0);


            v.SetZeros();
            for (size_t i = j + 1; i < row; i++)
            {
                for (size_t k = 0; k < j; k++)
                {
                    v(k, 0) = A(i, k)*A(k, k);
                }

                tmp = matRowLeftPivot * v;
                tmpValue = tmp(0, 0);
                A(i, j) = (A(i, j) - tmpValue) / A(j, j);

            }

        }
        else
        {
            for (size_t i = 1; i < row; i++)
            {
                A(i, j) = A(i, j) / A(j, j);
            }
        }
    }


    Matrix<T> L(A.GetNumRow(), A.GetNumCol());
    Matrix<T> D(A.GetNumRow(), A.GetNumCol());
    L.SetZeros();
    D.SetZeros();

    for (size_t i = 0; i < A.GetNumRow(); i++)
    {
        for (size_t j = 0; j < A.GetNumCol(); j++)
        {
            if (i > j)
            {
                L(i, j) = A(i, j);
            }
            else if (i == j)
            {
                L(i, i) = 1.0;
                D(i, j) = A(i, j);
            }
        }
    }

    vector<Matrix<T>> RC;
    RC.push_back(L);
    RC.push_back(D);
    return RC;
}


#endif
