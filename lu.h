#ifndef __LU_H__
#define __LU_H__

//������࣬�������㷨��ĸ���

#include <iostream>
#include <vector>
#include "Matrix.h"
#include "transform.h"
using namespace std;

template <typename T>
class LU
{
public:
	LU();
	LU(const Matrix<T>& A);
	LU(const Matrix<T>& A, const Matrix<T>& b);

	~LU();

	//����
	int Decompose();   //Ĭ�϶�A������зֽ�
	int Solve();         //��ⷽ����   Ax  =  b

	int Det(T* Val);      //����A������ʽֵ

	int GetL(Matrix<T>& L) const;
	int GetU(Matrix<T>& U) const;

private:

	int LUDecomposeFlag;
	int FirstTranFormTimes;

	Matrix<T> A;
	Matrix<T> b;

	Matrix<T> L;
	Matrix<T> U;

	Matrix<T> x;

};


template <typename T>
LU<T>::LU() :FirstTranFormTimes(0), LUDecomposeFlag(0)
{

}


template<typename T>
LU<T>::LU(const Matrix<T>& A_) :A(A_), FirstTranFormTimes(0), LUDecomposeFlag(0)
{
	
}

template<typename T>
LU<T>::LU(const Matrix<T>& A_, const Matrix<T>& b_) :A(A_), b(b_), FirstTranFormTimes(0), LUDecomposeFlag(0)
{

}

template<typename T>
LU<T>::~LU()
{

}

template <typename T>
int LU<T>::Decompose()
{
	std::vector<int> RowReturn;
	std::vector<int> ColReturn;
	std::vector<double> NumReturn;

	int RowNum = A.GetNumRow();
	int ColNum = A.GetNumCol();

	U.Resize(RowNum, ColNum);
	L.Resize(RowNum, ColNum);
	U = A;

	ToRowEchelonForm<T>(U, RowReturn, ColReturn, NumReturn, FirstTranFormTimes);

	L.IdentityMatrix();
	// ��λ����

	for (int i = 0; i < RowReturn.size(); i++)
	{
		int rowreturn = RowReturn[i];
		int colreturn = ColReturn[i];
		double numreturn = NumReturn[i];
		L(rowreturn, colreturn) = -numreturn;
	}

	LUDecomposeFlag = 1;
	return 1;
}

template <typename T>
int LU<T>::Det(T *Val)
{
	//�������ʽ
	int RowNum = A.GetNumRow();
	int ColNum = A.GetNumCol();

	if (RowNum != ColNum)
	{
		cout << "���󣬷Ƿ��󣬲��ɼ�������ʽ" << endl;
		exit(1);
	}
	else if (RowNum == 1)
	{
		//ֻ��һ��Ԫ��
		(*Val) = A(0, 0);
		return 1;
	}

	if (Decompose() == 0)
	{
		printf("LU�ֽ�ʧ��\n");
		return 0;
	}

	(*Val) = U(0,0);

	for (int i = 1; i < RowNum; i++)
	{
		(*Val) *= U(i, i);
	}

	if (FirstTranFormTimes % 2 != 0)
	{
		//���������ε�һ��任
		(*Val) = -(*Val);
	}
	return 1;
}

template <typename T>
int LU<T>::GetL(Matrix<T>& LReturn) const
{
	if (LUDecomposeFlag == 0)
	{
		printf("δ����LU�ֽ�\n");
		return 0;
	}
	else
	{
		LReturn = L;
		return 1;
	}
}

template <typename T>
int LU<T>::GetU(Matrix<T>& UReturn) const
{
	if (LUDecomposeFlag == 0)
	{
		printf("δ����LU�ֽ�\n");
			return 0;
	}
	else
	{
		UReturn = U;
		return 1;
	}
}


#endif