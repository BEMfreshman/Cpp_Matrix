#ifndef __MATRIX_H_
#define __MATRIX_H_




#include <iostream>
#include <assert.h>
#include <cmath>

#define EPS 1e-10

using namespace std;

template<typename T>
class Matrix
{
public:
    Matrix();
    Matrix(int Row, int Col);
    Matrix(const Matrix<T>& mat);
    Matrix<T>& operator=(const Matrix<T>& mat);  //赋值操作
    Matrix(int Row, int Col, T** arr);
    Matrix(int Row, int Col, T* arr, const string& StorageType);
    //采用数组初始化
    //第一二个参数为本矩阵的行数和列数
    //第三个参数为矩阵初始化的数组
    //第四个为第三个参数中的元素个数

	//初始化操作符


    /***********
     *重载操作符*
     **********/

    const Matrix<T> operator +(const Matrix<T>& mat) const; //矩阵加矩阵
    const Matrix<T> operator +(const T num) const; //矩阵加标量
    const Matrix<T> operator -(const Matrix<T>& mat) const;
    const Matrix<T> operator -(const T num) const;
    const Matrix<T> operator *(const Matrix<T>& mat) const;
    const Matrix<T> operator *(const T num) const;

    Matrix<T>& operator +=(const Matrix<T>& mat);
    Matrix<T>& operator +=(const T num);
    Matrix<T>& operator -=(const Matrix<T>& mat);
    Matrix<T>& operator -=(const T num);
    Matrix<T>& operator *=(const T num);
	Matrix<T>& operator /= (const T num);

    T& operator ()(int index_row,int index_col);//将函数操作符重载，实现寻址操作符功能
    const T operator()(int index_row,int index_col) const;//供常对象使用

    friend ostream& operator << <T>(ostream& os,const Matrix<T>& mat);//重载输出操作符

	
	

    /**************
     *基本的矩阵操作*
     **************/
    void SetZeros();   //所有元素置零
    void IdentityMatrix();  // 单位矩阵
	void Resize(int Row, int Col); //重新分配


	Matrix<T>& TransPose();

	Matrix<T>& FirstTypeTransForm(int Row_One, int Row_Two);
	const Matrix<T> FirstTypeTransForm(int Row_One, int Row_Two) const;
	//第一类初等变换---交换两行

	Matrix<T>& SecondTypeTransForm(int Row, double Num);
	const Matrix<T> SecondTypeTransForm(int Row, double Num) const;
	//第二类初等变换---某一行乘以一个数
	// Row---待乘数的行
	// Num---需要乘以的数

	Matrix<T>& ThirdTypeTransForm(int Row_One, int Row_Two, double Num);
	const Matrix<T> ThirdTypeTransForm(int Row_One, int Row_Two, double Num) const;
	//第三类变换----某一行的倍数加到另一行
	//Row_One ----- 需要乘数字的一行
	//Row_Two ----- 待加的一行
	//Num    ------ 乘数


	const Matrix<T> ExtractBlock(int RowStart, int ColStart, int RowNumToBlock, int ColNumToBlock) const;
	void SetBlock(int RowStart, int ColStart, int RowNumToSet, int ColNumtoSet, const Matrix<T>& mat);
		


    inline int GetNumRow();
	inline int GetNumRow() const;
    inline int GetNumCol();
	inline int GetNumCol() const;
    inline int GetNumData();

    


	~Matrix();

private:
	int NumRow;
	int NumCol;
	int Size;
	T** p1;

	void Allocate(int Num_Row,int Num_Col);
	void DeAllocate();
	
};




template<typename T>
void Matrix<T>::Allocate(int NumRow,int NumCol)
{
    //分配内存
    if(NumRow <= 0 || NumCol <=0)
    {
        //表示没有设置矩阵的维数，打印信息，不进行分配
        cout << "没有设置矩阵的维数，请检查" << endl;
        exit(0);
    }
    else
    {
        Size = NumRow * NumCol;
        p1 = new T* [NumRow];
        assert(p1 != NULL);
        for (int i = 0; i < NumRow; ++i)
        {
            p1[i] = new T[NumCol]; // 指向二维数组每行的开头位置
            assert(p1[i] != NULL);
        }
		this->NumRow = NumRow;
		this->NumCol = NumCol;
    }
}

template<typename T>
void Matrix<T>::DeAllocate()
{
    //释放矩阵所占有的内存资源
    //并将NumRow，NumCol重置为0

    for(int i = 0;i < NumRow;++i)
    {
        delete [] p1[i];
        p1[i] = NULL;
    }
    delete [] p1;
    p1 = NULL;
    NumRow = 0;
    NumCol = 0;
    Size = 0;
}

template<typename T>
Matrix<T>::Matrix() :NumRow(0), NumCol(0), Size(0), p1(NULL)
{
    
}

template<typename T>
Matrix<T>::Matrix(int Row, int Col) :NumRow(Row), NumCol(Col)
{
    Allocate(NumRow,NumCol);
}
template<typename T>
Matrix<T>::Matrix(const Matrix& mat) : NumRow(mat.NumRow), NumCol(mat.NumCol)
{
    Allocate(NumRow,NumCol);
    for (int i = 0; i < NumRow; ++i)
    {
        for (int j = 0; j < NumCol; ++j)
        {
            p1[i][j] = mat.p1[i][j];
        }
    }
}

template<typename T>
Matrix<T>& Matrix<T>::operator =(const Matrix<T>& mat)
{
    //需深拷贝
    if(&mat != this)
    {
        //防止自赋值
        DeAllocate();
        Allocate(mat.NumRow,mat.NumCol);
        for (int i = 0; i < NumRow; ++i)
        {
            for (int j = 0; j < NumCol; ++j)
            {
                p1[i][j] = mat.p1[i][j];
            }
        }
    }
    return *this;
}

template<typename T>
Matrix<T>::Matrix(int Row, int Col, T** arr):NumRow(Row),NumCol(Col)
{

    Allocate(NumRow,NumCol);  //分配存储空间
    for (int i = 0;i < NumRow;++i)
    {
        for(int j = 0;j < NumCol;++j)
        {
            p1[i][j] = arr[i][j];
        }
    }
}

template<typename T>
Matrix<T>::Matrix(int Row,int Col,T * arr,const string& Storage)
    :NumRow(Row),NumCol(Col)
{
    Allocate(NumRow,NumCol);
    if(Storage == "Row")
    {
        for(int i= 0;i < Row;i++)
        {
            for(int j = 0;j < Col;j++)
            {
                p1[i][j] = arr[Col * i + j];
            }
        }
    }
}


template<typename T>
void Matrix<T>::SetZeros()
{
    for(int i = 0;i < NumRow;++i)
    {
        for (int j = 0;j < NumCol;++j)
        {
            p1[i][j] = T (0);
        }
    }
}

template<typename T>
void Matrix<T>::IdentityMatrix()
{
    this->SetZeros();
    //首先全部置为0；
    for(int i = 0;i < NumRow;i++)
    {
        p1[i][i] = T (1);
    }
}


template<typename T>
Matrix<T>::~Matrix()
{
    DeAllocate();
}

template<typename T>
const Matrix<T> Matrix<T>::operator +(const Matrix<T>& mat) const
{
    if(NumRow != mat.NumRow || NumCol != mat.NumCol)
    {
        cout << "两个矩阵之间的维数不相等，程序退出" << endl;
        exit(0);
    }
    Matrix<T> res_mat(NumRow,NumCol);

    for(int i = 0;i < NumRow;++i)
    {
        for(int j = 0;j < NumCol;++j)
        {
            res_mat.p1[i][j] = p1[i][j] + mat.p1[i][j];
        }
    }
    return res_mat;
}

template<typename T>
const Matrix<T> Matrix<T>::operator +(const T num) const
{
    //向量与标量做加法
    Matrix<T> res_mat(NumRow,NumCol);
    for(int i = 0;i < NumRow;++i)
    {
        for(int j = 0;j < NumCol;++j)
        {
            res_mat.p1[i][j] = p1[i][j] + num;
        }
    }
    return res_mat;
}

template<typename T>
const Matrix<T> Matrix<T>::operator -(const Matrix<T>& mat) const
{
    if(NumRow != mat.NumRow || NumCol != mat.NumCol)
    {
        cout << "两个矩阵之间的维数不相等，程序退出" << endl;
        exit(0);
    }
    Matrix<T> res_mat(NumRow,NumCol);

    for(int i = 0;i < NumRow;++i)
    {
        for(int j = 0;j < NumCol;++j)
        {
            res_mat.p1[i][j] = p1[i][j] - mat.p1[i][j];
        }
    }
    return res_mat;
}

template<typename T>
const Matrix<T> Matrix<T>::operator -(const T num) const
{
    //向量与标量做加法
    Matrix<T> res_mat(NumRow,NumCol);
    for(int i = 0;i < NumRow;++i)
    {
        for(int j = 0;j < NumCol;++j)
        {
            res_mat.p1[i][j] = p1[i][j] - num;
        }
    }
    return res_mat;
}

template<typename T>
const Matrix<T> Matrix<T>::operator *(const Matrix<T>& mat) const
{
    if(NumCol !=mat.NumRow)
    {
        cout << "第一个矩阵的列数与第二个矩阵的行数不相等，程序退出" << endl;
        exit(0);
    }

    Matrix<T> res_mat(NumRow,mat.NumCol);
	res_mat.SetZeros();

    for(int i = 0;i < NumRow;++i)
    {
        for(int j = 0;j < mat.NumCol;++j)
        {
            for(int k = 0;k < NumCol;++k)
            {
                res_mat.p1[i][j] += p1[i][k]*mat.p1[k][j];
            }
        }
    }
    return res_mat;
}

template<typename T>
const Matrix<T> Matrix<T>::operator *(const T num) const
{
    Matrix<T> res_mat(NumRow,NumCol);
    for(int i = 0;i < NumRow;++i)
    {
        for(int j = 0;j < NumCol;++j)
        {
            res_mat.p1[i][j] = p1[i][j]*num;
        }
    }
    return res_mat;
}

template<typename T>
Matrix<T>& Matrix<T>::operator +=(const Matrix<T>& mat)
{
    if(NumRow != mat.NumRow || NumCol != mat.NumCol)
    {
        cout << "两个矩阵之间的维数不相等，程序退出" << endl;
        exit(0);
    }
    for(int i = 0;i < NumRow;++i)
    {
        for(int j = 0;j < NumCol;++j)
        {
            p1[i][j] += mat.p1[i][j];
        }
    }
    return *this;
}

template<typename T>
Matrix<T>& Matrix<T>::operator +=(const T num)
{
    for(int i = 0;i < NumRow;++i)
    {
        for(int j = 0;j < NumCol;++j)
        {
            p1[i][j] += num;
        }
    }
    return *this;
}

template<typename T>
Matrix<T>& Matrix<T>::operator -=(const Matrix<T>& mat)
{
    if(NumRow != mat.NumRow || NumCol != mat.NumCol)
    {
        cout << "两个矩阵之间的维数不相等，程序退出" << endl;
        exit(0);
    }
    for(int i = 0;i < NumRow;++i)
    {
        for(int j = 0;j < NumCol;++j)
        {
            p1[i][j] -= mat.p1[i][j];
        }
    }
    return *this;
}

template<typename T>
Matrix<T>& Matrix<T>::operator -=(const T num)
{
    for(int i = 0;i < NumRow;++i)
    {
        for(int j = 0;j < NumCol;++j)
        {
            p1[i][j] -= num;
        }
    }
    return *this;
}


template<typename T>
Matrix<T>& Matrix<T>::operator *=(const T num)
{
    for(int i = 0;i < NumRow;++i)
    {
        for(int j = 0;j < NumCol;++j)
        {
            p1[i][j] *= num;
        }
    }
    return *this;
}

template<typename T>
Matrix<T>& Matrix<T>::operator /= (const T num)
{
	for (int i = 0; i < NumRow; ++i)
	{
		for (int j = 0; j < NumCol; ++j)
		{
			p1[i][j] /= num;
		}
	}
	return *this;
}

template<typename T>
T& Matrix<T>::operator ()(int index_row,int index_col)
{
    return p1[index_row][index_col];
}

template<typename T>
const T Matrix<T>::operator ()(int index_row,int index_col) const
{
    return p1[index_row][index_col];
}


template<typename T>
inline int Matrix<T>::GetNumRow()
{
    return NumRow;
}

template<typename T>
inline int Matrix<T>::GetNumRow() const
{
	return NumRow;
}

template<typename T>
inline int Matrix<T>::GetNumCol()
{
    return NumCol;
}

template<typename T>
inline int Matrix<T>::GetNumCol() const
{
	return NumCol;
}

template<typename T>
inline int Matrix<T>::GetNumData()
{
    return Size;
}

template<typename T>
const Matrix<T> Matrix<T>::ExtractBlock(int RowStart, int ColStart, int RowNumToBlock, int ColNumToBlock) const
{
	//抽取一个Block出来

	assert(RowStart + RowNumToBlock <= NumRow);
	assert(ColStart + ColNumToBlock <= NumCol);


	Matrix<T> tmp(RowNumToBlock, ColNumToBlock);
	for (int i = RowStart; i < RowStart+ RowNumToBlock; i++)
	{
		for (int j = ColStart; j < ColStart + ColNumToBlock; j++)
		{
			tmp.p1[i - RowStart][j - ColStart] = p1[i][j];
		}
	}
	return tmp;
}

template<typename T>
void Matrix<T>::SetBlock(int RowStart, int ColStart, int RowNumToSet, int ColNumToSet, const Matrix<T>& mat)
{
	assert(RowStart + RowNumToSet <= NumRow);
	assert(ColStart + ColNumToSet <= NumCol);

	assert(RowNumToSet <= mat.NumRow);
	assert(ColNumToSet <= mat.NumCol);

	for (int i = 0; i < RowNumToSet; i++)
	{
		for (int j = 0; j < ColNumToSet; j++)
		{
			p1[i + RowStart][j + ColStart] = mat.p1[i][j];
		}
	}
}


template<typename T>
Matrix<T>& Matrix<T>::FirstTypeTransForm(int Row_One, int Row_Two)
{
	//交换两行

	T * row_tmp = NULL;
	row_tmp = p1[Row_One];
	p1[Row_One] = p1[Row_Two];
	p1[Row_Two] = row_tmp;

	return *this;
}


template<typename T>
const Matrix<T> Matrix<T>::FirstTypeTransForm(int Row_One, int Row_Two) const
{
	//交换两行

	Matrix<T> tmpMat = *this;

	T * row_tmp = NULL;
	row_tmp = tmpMat.p1[Row_One];
	tmpMat.p1[Row_One] = tmpMat.p1[Row_Two];
	tmpMat.p1[Row_Two] = row_tmp;

	return tmpMat;
}



template<typename T>
Matrix<T>& Matrix<T>::SecondTypeTransForm(int Row, double Num)
{
//第二类初等变换---某一行乘以一个数
// Row---待乘数的行
// Num---需要乘以的数
	
	T * row_tmp = p1[Row];

	for (int i = 0; i < NumCol; i++)
	{
		row_tmp[i] *= Num;
	}
	return *this;
}


template<typename T>
const Matrix<T> Matrix<T>::SecondTypeTransForm(int Row, double Num) const
{
	//第二类初等变换---某一行乘以一个数
	// Row---待乘数的行
	// Num---需要乘以的数
	Matrix<T> tmpMat = *this;

	T * row_tmp = tmpMat.p1[Row];

	for (int i = 0; i < NumCol; i++)
	{
		row_tmp[i] *= Num;
	}
	return tmpMat;

}



template<typename T>
Matrix<T>& Matrix<T>::ThirdTypeTransForm(int Row_One, int Row_Two, double Num)
{
	//第三类变换----某一行的倍数加到另一行
	//Row_One ----- 需要乘数字的一行
	//Row_Two ----- 待加的一行
	//Num    ------ 乘数


	T *row_one = p1[Row_One];
	T *row_two = p1[Row_Two];

	for (int i = 0; i < NumCol; i++)
	{
		row_two[i] += row_one[i] * Num;
	}
	return *this;
}


template<typename T>
const Matrix<T> Matrix<T>::ThirdTypeTransForm(int Row_One, int Row_Two, double Num) const
{
	//第三类变换----某一行的倍数加到另一行
	//Row_One ----- 需要乘数字的一行
	//Row_Two ----- 待加的一行
	//Num    ------ 乘数

	Matrix<T> tmpMat = *this;

	T *row_one = tmpMat.p1[Row_One];
	T *row_two = tmpMat.p1[Row_Two];

	for (int i = 0; i < NumCol; i++)
	{
		row_two[i] += row_one[i] * Num;
	}
	return tmpMat;
}


template<typename T>
ostream& operator<<(ostream& os,const Matrix<T>& mat)
{
	int Row = mat.GetNumRow();
	int Col = mat.GetNumCol();

	for (int i = 0; i < Row; i++)
	{
		for (int j = 0; j < Col; j++)
		{
			os << mat(i, j) << " ";
		}
		os << endl;
	}
	return os;
}

template<typename T>
void Matrix<T>::Resize(int Row, int Col)
{
	if (p1 != NULL)
	{
		DeAllocate();
	}
	Allocate(Row, Col);
}

template<typename T>
Matrix<T>& Matrix<T>::TransPose()
{
	T tmp;
	for (int i = 0; i < NumRow; i++)
	{
		for (int j = i; j < NumCol; j++)
		{
			tmp = p1[i][j];
			p1[i][j] = p1[j][i];
			p1[j][i] = tmp;
		}
	}
	return *this;
}



#endif   // MATRIX_H_
