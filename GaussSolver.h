#ifndef _GAUSSSOLVER_H__
#define _GAUSSSOLVER_H_

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

	int row = A.GetNumRow();

	ans.Resize(b.GetNumRow(), 1);
	for (int i = 0; i < A.GetNumRow() - 1; i++)
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

	return 1;

}

template<typename T>
const Matrix<T> GaussSolver<T>::getAns() const
{
	return ans;
}



#endif