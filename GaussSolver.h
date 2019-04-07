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
