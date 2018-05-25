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

	int FullPivotSolve();          //ȫѡ��Ԫ��ȥ
	int ColPivotSolve();           //����Ԫ��ȥ


private:
	Matrix<T> A;
	Matrix<T> b;
	Matrix<T> ans;
};



#endif