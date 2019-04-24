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
    GaussSolver() {};
	GaussSolver(const Matrix<T>& A, const Matrix<T>& b);
    GaussSolver<T> &operator=(const GaussSolver<T>& GS);


	const Matrix<T> Solve() const;              //一般的LU解法
	const Matrix<T> FullPivotSolve() const;          //全选主元消去
	const Matrix<T> ColPivotSolve() const;           //列主元消去


private:
	Matrix<T> A;
	Matrix<T> b;
};

template<typename T>
GaussSolver<T>::GaussSolver(const Matrix<T>& A, const Matrix<T>& b)
{

	if(!A.isSquare())
	{
		throw runtime_error("A is not square matrix");
	}

	if(A.GetNumRow() != b.GetNumRow())
	{
		throw runtime_error("The Row Number of A is not equal to b");
	}

	if(b.GetNumCol() != 1)
	{
		throw runtime_error("b is not Col Vector");
	}


	this->A = A;
	this->b = b;
}

template <typename T>
GaussSolver<T>& GaussSolver<T>::operator=(const GaussSolver<T>& GS)
{
    this->A = GS.A;
    this->b = GS.b;
}


template<typename T>
const Matrix<T> GaussSolver<T>::FullPivotSolve()  const        //全选主元消去
{


}

template<typename T>
const Matrix<T> GaussSolver<T>::ColPivotSolve()  const         //列主元消去
{



}

template <typename T>
const Matrix<T> GaussSolver<T>::Solve() const
{
	LU<T> lu(A);
	Matrix<T> L(A.GetNumRow(),A.GetNumCol());
	Matrix<T> U(A.GetNumRow(),A.GetNumCol());

	try
	{
		vector<Matrix<T>> LUMatrix = lu.LUDeCompose();
		L = LUMatrix[0];
		U = LUMatrix[1];
	}
	catch (runtime_error& e)
	{
		throw runtime_error("LU failed");
	}

	Matrix<T> ans(A.GetNumRow(),1);


//	cout << "L" << endl;
//	cout << L << endl;
//
//
//	cout << "U" << endl;
//	cout << U << endl;
//
//	cout << "L * U" << endl;
//	cout << L * U << endl;
//
//	cout << "b" << endl;
//	cout << b << endl;

	try
	{
		ans = Utility<T>::UptriSolve(U,Utility<T>::LowtriSolve(L,b));
	}
	catch (runtime_error& e)
	{
		throw runtime_error("LU solve failed");
	}

	return ans;
}



#endif
