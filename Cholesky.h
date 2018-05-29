#ifndef __CHOLESKY_H__
#define __CHOLESKY_H__

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

	int DefaultFact();      //A = L * L';
	int AdvanceFact();      //A = L * D * L';

	const Matrix<T> GetL() const;
	const Matrix<T> GetD() const;

private:
	Matrix<T> A;

	Matrix<T> L;
	Matrix<T> D;


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
int Cholesky<T>::DefaultFact()
{
	//A = L * L';
	int row = A.GetNumRow();

	Matrix<T> matColBelowPivot;
	Matrix<T> matRowLeftPivot;
	Matrix<T> matLeftBelowPivot;


	Matrix<T> mattmp;

	for (int i = 0; i < A.GetNumCol(); i++)
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

			cout << "A" << endl;
			cout << A << endl;

			A(i, i) = sqrt(A(i, i));

			matColBelowPivot = A.ExtractBlock(i + 1, i, row - i - 1, 1);
			matColBelowPivot /= A(i, i);

			A.SetBlock(i + 1, i, row - i - 1, 1, matColBelowPivot);

			cout << "A" << endl;
			cout << A << endl;
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

	

	L.Resize(A.GetNumRow(), A.GetNumCol());
	L.SetZeros();

	for (int i = 0; i < A.GetNumRow(); i++)
	{
		for (int j = 0; j < A.GetNumCol(); j++)
		{
			if (i >= j)
			{
				L(i, j) = A(i, j);
			}
		}
	}

	return 1;

}

template<typename T>
const Matrix<T> Cholesky<T>::GetL() const
{
	return L;
}




#endif