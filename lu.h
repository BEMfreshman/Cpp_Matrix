#ifndef __LU_H__
#define __LU_H__

//抽象基类，是所有算法类的父类

#include <iostream>
#include <vector>
#include "Matrix.h"
#include "transform.h"
#include "GaussSolver.h"

using namespace std;

template <typename T>
class LU
{
public:
	LU();
	LU(const Matrix<T>& A);
	LU(const Matrix<T>& A, const Matrix<T>& b);

	~LU();

	//操作
	//int Decompose();   //默认对A矩阵进行分解

	int DefaultFact();        //默认分解方法
	int FullPivotFact();      //全主元分解方法
	int ColPivotFact();       //列主元分解方法

	int Solve();         //求解方程组   Ax  =  b

	int Det(T* Val);      //计算A的行列式值

	int GetL(Matrix<T>& L) const;
	int GetU(Matrix<T>& U) const;

	int GetP(Matrix<T>& P) const;
	int GetQ(Matrix<T>& Q) const;


private:

	int LUDecomposeFlag;
	int FirstTranFormTimes;

	Matrix<T> A;
	Matrix<T> b;

	Matrix<T> L;
	Matrix<T> U;

	Matrix<T> P;
	Matrix<T> Q;

private:
	const Matrix<T> ProducePorQMatrix(int p, int q);
	//交换第p行（列）和第q行（列）
	//理论参见 P21

	int ithGaussFact(int i);
	int ithGaussFact(int i, Matrix<T>& Li);
	//选取第i行第i列作为主元进行高斯消去


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

template<typename T>
const Matrix<T> LU<T>::ProducePorQMatrix(int p, int q)
{
	Matrix<T> mat(A.GetNumRow(), A.GetNumCol());
	mat.IdentityMatrix();     //单位矩阵

	//交换mat的第p列和第q列

	Matrix<T> pmat(A.GetNumRow(), 1);
	Matrix<T> qmat(A.GetNumRow(), 1);

	pmat = mat.ExtractBlock(0, p, pmat.GetNumRow(), pmat.GetNumCol());
	qmat = mat.ExtractBlock(0, q, qmat.GetNumRow(), qmat.GetNumCol());

	mat.SetBlock(0, p, pmat.GetNumRow(), pmat.GetNumCol(), qmat);
	mat.SetBlock(0, q, qmat.GetNumRow(), qmat.GetNumCol(), pmat);

	return mat;

}

template<typename T>
int LU<T>::ithGaussFact(int i)
{
	int row = A.GetNumRow();

	Matrix<T> matColBelowPivot(row - i - 1, 1);
	Matrix<T> matRowRightPivot(1, row - i - 1);
	Matrix<T> matBelowRightPivot(row - i - 1, row - i - 1);

	double pivot = A(i, i);
	if (abs(pivot) < sqrt(EPS))
	{
		//主元为0
		return 0;
	}

	matColBelowPivot = A.ExtractBlock(i + 1, i, matColBelowPivot.GetNumRow(), matColBelowPivot.GetNumCol());
	matRowRightPivot = A.ExtractBlock(i, i + 1, matRowRightPivot.GetNumRow(), matRowRightPivot.GetNumCol());
	matBelowRightPivot = A.ExtractBlock(i + 1, i + 1, matBelowRightPivot.GetNumRow(), matBelowRightPivot.GetNumCol());

	matColBelowPivot /= pivot;
	matBelowRightPivot -= matColBelowPivot * matRowRightPivot;


	A.SetBlock(i + 1, i, matColBelowPivot.GetNumRow(), matColBelowPivot.GetNumCol(), matColBelowPivot);
	A.SetBlock(i + 1, i + 1, matBelowRightPivot.GetNumRow(), matBelowRightPivot.GetNumCol(), matBelowRightPivot);

	
	return 1;
}

template<typename T>
int LU<T>::ithGaussFact(int i,Matrix<T>& InvLi)
{
	int row = A.GetNumRow();

	InvLi.Resize(A.GetNumRow(), A.GetNumCol());
	InvLi.IdentityMatrix();

	Matrix<T> matColBelowPivot(row - i - 1, 1);
	Matrix<T> matRowRightPivot(1, row - i - 1);
	Matrix<T> matBelowRightPivot(row - i - 1, row - i - 1);

	double pivot = A(i, i);
	if (abs(pivot) < sqrt(EPS))
	{
		//主元为0
		return 0;
	}

	matColBelowPivot = A.ExtractBlock(i + 1, i, matColBelowPivot.GetNumRow(), matColBelowPivot.GetNumCol());
	matRowRightPivot = A.ExtractBlock(i, i + 1, matRowRightPivot.GetNumRow(), matRowRightPivot.GetNumCol());
	matBelowRightPivot = A.ExtractBlock(i + 1, i + 1, matBelowRightPivot.GetNumRow(), matBelowRightPivot.GetNumCol());

	matColBelowPivot /= pivot;
	matBelowRightPivot -= matColBelowPivot * matRowRightPivot;

	InvLi.SetBlock(i + 1, i, matColBelowPivot.GetNumRow(), matColBelowPivot.GetNumCol(), matColBelowPivot);
	
	/*
	Li本身是单位矩阵经过初等行变换得来的矩阵，所以InvLi（即它的逆矩阵）为其非对角线元素皆取负数即可
	
	Li = [                                InvLi = [
	        1   0   0                               1    0    0
	        2   1   0                              -2    1    0    
			3   0   1                              -3    0    1
	     ]                                        ]
	
	
	*/

	

	A.SetBlock(i + 1, i, matColBelowPivot.GetNumRow(), matColBelowPivot.GetNumCol(), matColBelowPivot);
	A.SetBlock(i + 1, i + 1, matBelowRightPivot.GetNumRow(), matBelowRightPivot.GetNumCol(), matBelowRightPivot);

	
	return 1;
}

template<typename T>
int LU<T>::DefaultFact()
{
	//当主元为0时会失败
	//理论《数值线性代数》 P18

	int row = A.GetNumRow();
	int col = A.GetNumCol();

	L.Resize(row, col);
	U.Resize(row, col);
	L.SetZeros();
	U.SetZeros();


	for (int i = 0; i < row - 1; i++)
	{
		int reFlag = ithGaussFact(i);
		if (reFlag == 0)
		{
			return 0;
		}
	}

	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < col; j++)
		{
			if (i > j)
			{
				L(i, j) = A(i, j);
			}
			else if (i == j)
			{
				L(i, j) = 1;
				U(i, j) = A(i, j);
			}
			else 
			{
				U(i, j) = A(i, j);
			}
		}
	}
	LUDecomposeFlag = 1;
	return 1;

}

template <typename T>
int LU<T>::FullPivotFact()
{
	//全选主元三角分解
	//具体理论参见《数值线性代数》 P21——P25
	//对A矩阵的分解结果为
	//PAQ = LU

	int row = A.GetNumRow();
	int col = A.GetNumCol();

	
	P.Resize(row, col);
	P.SetZeros();
	Q.Resize(row, col);
	Q.SetZeros();
	L.Resize(row, col);
	L.SetZeros();
	U.Resize(row, col);
	U.SetZeros();
	

	vector<Matrix<T>> PVec;
	vector<Matrix<T>> QVec;
	//vector<Matrix<T>> InvLiVec;

	T Pivot;
	int PivotRow;
	int PivotCol;
	for (int i = 0; i < row - 1; i++)
	{
		Matrix<T> PEach;
		Matrix<T> QEach;
		Matrix<T> InvLiEach;
		A.FindMax(i, i, row - i, col - i, &Pivot, &PivotRow, &PivotCol);

		if (abs(Pivot) < sqrt(EPS))
		{
			printf("矩阵奇异");
			return 0;
		}

		PEach = ProducePorQMatrix(PivotRow, i);
		QEach = ProducePorQMatrix(PivotCol, i);

		cout << "Before" << endl;
		cout << A << endl;
		A = PEach * A * QEach;

		cout << "After" << endl;
		cout << A << endl;

		PVec.push_back(PEach);
		QVec.push_back(QEach);

		int reFlag = ithGaussFact(i);
		if (reFlag == 0)
		{
			return 0;
		}

		cout << "After i Fact" << endl;
		cout << A << endl;

		//InvLiVec.push_back(InvLiEach);
	}


	Matrix<T> InvLiPi;
	for (int i = 0; i < QVec.size(); i++)
	{
		if (i == 0)
		{
			Q = QVec[i];
			//L = InvLiVec[i];
		}
		else
		{
			Q *= QVec[i];
			//L = PVec[i] * L*PVec[i] * InvLiVec[i];
		}
	}

	for (int i = PVec.size() - 1; i >= 0; i--)
	{
		if (i == PVec.size() - 1)
		{
			P = PVec[i];
		}
		else
		{
			P *= PVec[i];
		}
	}

	/*Matrix<T> tmpI(LiPi.GetNumRow(), LiPi.GetNumCol());
	tmpI.IdentityMatrix();

	GaussSolver<T> GS(LiPi, tmpI);
	GS.LowtriSolve();

	InvLiPi = GS.getAns();*/

	//L = P*InvLiPi;
	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < col; j++)
		{
			if (i > j)
			{
				L(i, j) = A(i, j);
			}
			else if (i == j)
			{
				L(i, j) = 1;
				U(i, j) = A(i, j);
			}
			else
			{
				U(i, j) = A(i, j);
			}
		}
	}

	LUDecomposeFlag = 1;
	return 1;
}

//template <typename T>
//int LU<T>::Decompose()
//{
//	std::vector<int> RowReturn;
//	std::vector<int> ColReturn;
//	std::vector<double> NumReturn;
//
//	int RowNum = A.GetNumRow();
//	int ColNum = A.GetNumCol();
//
//	U.Resize(RowNum, ColNum);
//	L.Resize(RowNum, ColNum);
//	U = A;
//
//	ToRowEchelonForm<T>(U, RowReturn, ColReturn, NumReturn, FirstTranFormTimes);
//
//	L.IdentityMatrix();
//	// 单位矩阵
//
//	for (int i = 0; i < RowReturn.size(); i++)
//	{
//		int rowreturn = RowReturn[i];
//		int colreturn = ColReturn[i];
//		double numreturn = NumReturn[i];
//		L(rowreturn, colreturn) = -numreturn;
//	}
//
//	LUDecomposeFlag = 1;
//	return 1;
//}

template <typename T>
int LU<T>::Det(T *Val)
{
	//求解行列式
	int RowNum = A.GetNumRow();
	int ColNum = A.GetNumCol();

	if (RowNum != ColNum)
	{
		cout << "错误，非方阵，不可计算行列式" << endl;
		exit(1);
	}
	else if (RowNum == 1)
	{
		//只有一个元素
		(*Val) = A(0, 0);
		return 1;
	}

	if (Decompose() == 0)
	{
		printf("LU分解失败\n");
		return 0;
	}

	(*Val) = U(0,0);

	for (int i = 1; i < RowNum; i++)
	{
		(*Val) *= U(i, i);
	}

	if (FirstTranFormTimes % 2 != 0)
	{
		//做了奇数次第一类变换
		(*Val) = -(*Val);
	}
	return 1;
}

template <typename T>
int LU<T>::GetL(Matrix<T>& LReturn) const
{
	if (LUDecomposeFlag == 0)
	{
		printf("未进行LU分解\n");
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
		printf("未进行LU分解\n");
			return 0;
	}
	else
	{
		UReturn = U;
		return 1;
	}
}

template<typename T>
int LU<T>::GetP(Matrix<T>& Pmat) const
{
	if (LUDecomposeFlag == 0)
	{
		printf("未进行LU分解\n");
		return 0;
	}
	else
	{
		Pmat = P;
		return 1;
	}
}

template<typename T>
int LU<T>::GetQ(Matrix<T>& Qmat) const
{
	if (LUDecomposeFlag == 0)
	{
		printf("未进行LU分解\n");
		return 0;
	}
	else
	{
		Qmat = Q;
		return 1;
	}
}


#endif