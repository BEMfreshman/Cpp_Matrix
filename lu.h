#ifndef LU_H
#define LU_H

#include <iostream>
#include <vector>
#include "Matrix.h"
#include "transform.h"
//#include "GaussSolver.h"

using namespace std;

template <typename T>
class LU
{
public:
	LU();
	LU(const Matrix<T>& A);
	LU(const Matrix<T>& A, const Matrix<T>& b);

	~LU();

	int Det(T* Val);      //����A������ʽֵ

//	int GetL(Matrix<T>& L) const;
//	int GetU(Matrix<T>& U) const;

//	int GetP(Matrix<T>& P) const;
//	int GetQ(Matrix<T>& Q) const;

    vector<Matrix<T>> LUDeCompose();
    vector<Matrix<T>> PQLUDeCompose();
    vector<Matrix<T>> PLUDeCompose();


private:

	int LUDecomposeFlag;
	int FirstTranFormTimes;

	Matrix<T> A;
	Matrix<T> b;

	Matrix<T> P;
	Matrix<T> Q;

private:
	const Matrix<T> ProducePorQMatrix(int p, int q);
	//������p�У��У��͵�q�У��У�
	//���۲μ� P21

	int ithGaussFact(int i);
	int ithGaussFact(int i, Matrix<T>& Li);
	//ѡȡ��i�е�i����Ϊ��Ԫ���и�˹��ȥ


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
	mat.IdentityMatrix();     //��λ����

	//����mat�ĵ�p�к͵�q��

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
		//��ԪΪ0
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
		//��ԪΪ0
		return 0;
	}

	matColBelowPivot = A.ExtractBlock(i + 1, i, matColBelowPivot.GetNumRow(), matColBelowPivot.GetNumCol());
	matRowRightPivot = A.ExtractBlock(i, i + 1, matRowRightPivot.GetNumRow(), matRowRightPivot.GetNumCol());
	matBelowRightPivot = A.ExtractBlock(i + 1, i + 1, matBelowRightPivot.GetNumRow(), matBelowRightPivot.GetNumCol());

	matColBelowPivot /= pivot;
	matBelowRightPivot -= matColBelowPivot * matRowRightPivot;

	InvLi.SetBlock(i + 1, i, matColBelowPivot.GetNumRow(), matColBelowPivot.GetNumCol(), matColBelowPivot);
	
	/*
	Li�����ǵ�λ���󾭹������б任�����ľ�������InvLi�������������Ϊ��ǶԽ���Ԫ�ؽ�ȡ��������
	
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

    int RowNum = A.GetNumRow();
    int ColNum = A.GetNumCol();

    Matrix<T> U(RowNum, ColNum);
    Matrix<T> L(RowNum, ColNum);
    //����ԪΪ0ʱ��ʧ��
    //���ۡ���ֵ���Դ����� P18

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

    vector<Matrix<T>> RC;
    RC.push_back(L,U);

    return RC;
}

template <typename T>
vector<Matrix<T>> LU<T>::PQLUDeCompose()
{

    //ȫѡ��Ԫ���Ƿֽ�
    //�������۲μ�����ֵ���Դ����� P21����P25
    //��A����ķֽ���Ϊ
    //PAQ = LU

    int row = A.GetNumRow();
    int col = A.GetNumCol();

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
    for (int i = 0; i < row - 1; i++)
    {
        Matrix<T> PEach;
        Matrix<T> QEach;
        Matrix<T> InvLiEach;
        A.FindMax(i, i, row - i, col - i, &Pivot, &PivotRow, &PivotCol);

        if (abs(Pivot) < sqrt(EPS))
        {
            printf("��������");
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
    RC.push_back(P,Q,L,U);
    return RC;
}

template <typename T>
vector<Matrix<T>> LU<T>::PLUDeCompose()
{

    //����Ԫ���Ƿֽ�
    //�������۲μ�����ֵ���Դ����� P26
    //��A����ķֽ���Ϊ
    //PA = LU

    int row = A.GetNumRow();
    int col = A.GetNumCol();


    Matrix<T> P(row, col);
    P.SetZeros();
    Matrix<T>L(row, col);
    L.SetZeros();
    Matrix<T>U(row, col);
    U.SetZeros();


    vector<Matrix<T>> PVec;

    T Pivot;
    int PivotRow;
    int PivotCol;
    for (int i = 0; i < row - 1; i++)
    {
        Matrix<T> PEach;
        Matrix<T> InvLiEach;
        A.FindMax(i, i, row - i, 1, &Pivot, &PivotRow, &PivotCol);

        if (abs(Pivot) < sqrt(EPS))
        {
            printf("��������");
            return 0;
        }

        PEach = ProducePorQMatrix(PivotRow, i);

        cout << "Before" << endl;
        cout << A << endl;
        A = PEach * A;

        cout << "After" << endl;
        cout << A << endl;

        PVec.push_back(PEach);

        int reFlag = ithGaussFact(i);
        if (reFlag == 0)
        {
            return 0;
        }

        cout << "After i Fact" << endl;
        cout << A << endl;

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
    RC.push_back(P,L,U);

    return RC;
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

    vector<Matrix<T>> LUMatrix = LUDeCompose();
    Matrix<T> U = LUMatrix[1];

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

    return EXIT_SUCCESS;
}

#endif
