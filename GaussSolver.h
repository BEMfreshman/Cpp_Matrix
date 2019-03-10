#ifndef GAUSSSOLVER_H
#define GAUSSSOLVER_H

#include <iostream>
#include <vector>
#include "Matrix.h"

using namespace std;

template<typename T>
class GaussSolver
{
public:
	GaussSolver(const Matrix<T>& A, const Matrix<T>& b);

	int FullPivotSolve();          //全选主元消去
	int ColPivotSolve();           //列主元消去


	int LowtriSolve();
	int UptriSolve();


	const Matrix<T> getAns() const;


private:
	Matrix<T> A;
	Matrix<T> b;
	Matrix<T> ans;
};

template<typename T>
GaussSolver<T>::GaussSolver(const Matrix<T>& A, const Matrix<T>& b)
{
	this->A = A;
	this->b = b;
	ans.Resize(b.GetNumRow(), b.GetNumCol());
}




template<typename T>
int GaussSolver<T>::LowtriSolve()
{
	//解下三角形矩阵
	//采用前带法

	//理论参见 《数值线性代数》 第二版 P12

    size_t row = A.GetNumRow();

	ans.Resize(b.GetNumRow(), 1);
    for (size_t i = 0; i < A.GetNumRow() - 1; i++)
	{
		Matrix<T> bn(b.GetNumRow() - i - 1, 1);
		Matrix<T> An(b.GetNumRow() - i - 1, 1);
		bn = b.ExtractBlock(i + 1, 0, bn.GetNumRow(), bn.GetNumCol());
		An = A.ExtractBlock(i + 1, i, An.GetNumRow(), An.GetNumCol());

		ans(i,0) = b(i,0) / A(i, i);

		bn -= An*ans(i, 0);

		b.SetBlock(i + 1, 0, bn.GetNumRow(), bn.GetNumCol(), bn);
	}

	ans(row - 1, 0) = b(row - 1, 0) / A(row - 1, row - 1);

    return EXIT_SUCCESS;

}

template <typename T>
int GaussSolver<T>::UptriSolve()
{
    //解上三角矩阵
    //采用回代法

    //理论参见 《数值线性代数》 第二版 P13

    ans.Resize(b.GetNumRow(),1);


    for(size_t i = A.GetNumRow() - 1; i >= 1; i--)
    {
        Matrix<T> bn(i,1);
        Matrix<T> An(i,1);

        bn = b.ExtractBlock(0,0,bn.GetNumRow(),bn.GetNumCol());
        An = A.ExtractBlock(0,i,An.GetNumRow(),An.GetNumCol());


        ans(i,0) = b(i,0) / A(i,i);
        bn -= An * ans(i,0);

        b.SetBlock(0,0,bn.GetNumRow(),bn.GetNumCol(),bn);
    }

    ans(0,0) = b(0,0) / A(0,0);
    return EXIT_SUCCESS;

}

template<typename T>
const Matrix<T> GaussSolver<T>::getAns() const
{
	return ans;
}

template<typename T>
int GaussSolver<T>::FullPivotSolve()          //全选主元消去
{



}

template<typename T>
int GaussSolver<T>::ColPivotSolve()           //列主元消去
{



}



#endif
