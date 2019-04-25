#ifndef LU_H
#define LU_H

#include <iostream>
#include <vector>
#include "Matrix.h"
#include "Utility.h"
//#include "GaussSolver.h"

using namespace std;

template <typename T>
class LU
{
public:
	LU();
    explicit LU(const Matrix<T>& A);

	~LU();

    vector<Matrix<T>> LUDeCompose();
    vector<Matrix<T>> PQLUDeCompose();
    vector<Matrix<T>> PLUDeCompose();

    int getFirstTransFormTimes() {return FirstTranFormTimes;};


private:

	int LUDecomposeFlag;
	int FirstTranFormTimes;

	Matrix<T> A;

private:
	const Matrix<T> ProducePorQMatrix(size_t p, size_t q);
	//交换第p行（列）和第q行（列）
	//理论参见 P21

	int ithGaussFact(size_t i);
	int ithGaussFact(size_t i, Matrix<T>& Li);
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
LU<T>::~LU()
{

}

template<typename T>
const Matrix<T> LU<T>::ProducePorQMatrix(size_t p, size_t q)
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
int LU<T>::ithGaussFact(size_t i)
{
	size_t row = A.GetNumRow();

	Matrix<T> matColBelowPivot(row - i - 1, 1);
	Matrix<T> matRowRightPivot(1, row - i - 1);
	Matrix<T> matBelowRightPivot(row - i - 1, row - i - 1);

	double pivot = A(i, i);
	if (abs(pivot) < sqrt(EPS))
	{
		//主元为0
		return EXIT_FAILURE;
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
int LU<T>::ithGaussFact(size_t i,Matrix<T>& InvLi)
{
	size_t row = A.GetNumRow();

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


template <typename T>
vector<Matrix<T>> LU<T>::LUDeCompose()
{
    std::vector<int> RowReturn;
    std::vector<int> ColReturn;
    std::vector<double> NumReturn;

    size_t RowNum = A.GetNumRow();
    size_t ColNum = A.GetNumCol();

    Matrix<T> U(RowNum, ColNum);
    Matrix<T> L(RowNum, ColNum);
    //当主元为0时会失败
    //理论《数值线性代数》 P18

    size_t row = A.GetNumRow();
    size_t col = A.GetNumCol();

    L.Resize(row, col);
    U.Resize(row, col);
    L.SetZeros();
    U.SetZeros();


    for (size_t i = 0; i < row - 1; i++)
    {
        int reFlag = ithGaussFact(i);
        if (reFlag == 0)
        {
            throw runtime_error("There are zero items in diag, LU failed");
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

    vector<Matrix<T>> RC;
    RC.push_back(L);
    RC.push_back(U);

    return RC;
}

template <typename T>
vector<Matrix<T>> LU<T>::PQLUDeCompose()
{

    //全选主元三角分解
    //具体理论参见《数值线性代数》 P21――P25
    //对A矩阵的分解结果为
    //PAQ = LU

    size_t row = A.GetNumRow();
    size_t col = A.GetNumCol();

    Matrix<T> P(row, col);
    P.SetZeros();
    Matrix<T> Q(row, col);
    Q.SetZeros();
    Matrix<T>L(row, col);
    L.SetZeros();
    Matrix<T>U(row, col);
    U.SetZeros();


    vector<Matrix<T>> PVec;
    vector<Matrix<T>> QVec;
    //vector<Matrix<T>> InvLiVec;

    T Pivot;
    int PivotRow;
    int PivotCol;
    for (size_t i = 0; i < row - 1; i++)
    {
        Matrix<T> PEach;
        Matrix<T> QEach;
        Matrix<T> InvLiEach;
        A.FindMax(i, i, row - i, col - i, &Pivot, &PivotRow, &PivotCol);

        if (abs(Pivot) < sqrt(EPS))
        {
            printf("矩阵奇异");
            return EXIT_FAILURE;
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
            return EXIT_FAILURE;
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

    for (size_t i = PVec.size() - 1; i >= 0; i--)
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


    for (size_t i = 0; i < row; i++)
    {
        for (size_t j = 0; j < col; j++)
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

    vector<Matrix<T>> RC;
    RC.push_back(P);
    RC.push_back(Q);
    RC.push_back(L);
    RC.push_back(U);
    return RC;
}

template <typename T>
vector<Matrix<T>> LU<T>::PLUDeCompose()
{

    //列主元三角分解
    //具体理论参见《数值线性代数》 P26
    //对A矩阵的分解结果为
    //PA = LU

    size_t row = A.GetNumRow();
    size_t col = A.GetNumCol();


    Matrix<T> P(row, col);
    P.SetZeros();
    Matrix<T>L(row, col);
    L.SetZeros();
    Matrix<T>U(row, col);
    U.SetZeros();


    vector<Matrix<T>> PVec;

    T Pivot;
    size_t PivotRow;
    size_t PivotCol;
    for (size_t i = 0; i < row - 1; i++)
    {
        Matrix<T> PEach;
        Matrix<T> InvLiEach;
        A.FindMax(i, i, row - i, 1, &Pivot, &PivotRow, &PivotCol);

        if (abs(Pivot) < sqrt(EPS))
        {
            printf("矩阵奇异");
            throw runtime_error("Matrix Singularity");
        }

        PEach = ProducePorQMatrix(PivotRow, i);

//        cout << "Before" << endl;
//        cout << A << endl;
        A = PEach * A;

//        cout << "After" << endl;
//        cout << A << endl;

        PVec.push_back(PEach);

        int reFlag = ithGaussFact(i);
        if (reFlag == 0)
        {
            throw runtime_error("Gauss Fact failed");
        }

//        cout << "After i Fact" << endl;
//        cout << A << endl;

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

    vector<Matrix<T>> RC;
    RC.push_back(P);
    RC.push_back(L);
    RC.push_back(U);

    return RC;
}

#endif
