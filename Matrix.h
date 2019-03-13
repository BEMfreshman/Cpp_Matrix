#ifndef MATRIX_H_
#define MATRIX_H_

//#include <petscmat.h>

#include <iostream>
#include <assert.h>
#include <cmath>
#include <vector>
#include <algorithm>

#include "Vector.h"

#define EPS 1e-10

using namespace std;

template<typename T>
class Matrix
{
public:
    Matrix();
    Matrix(size_t Row, size_t Col);
    Matrix(const Matrix<T>& mat);
    Matrix<T>& operator=(const Matrix<T>& mat);  //赋值操作
    Matrix(size_t Row, size_t Col, T** arr);
    Matrix(size_t Row, size_t Col, T* arr, const string& StorageType);
    //采用数组初始化
    //第一二个参数为本矩阵的行数和列数
    //第三个参数为矩阵初始化的数组
    //第四个为第三个参数中的元素个数

    /***********
     *重载操作符*
     **********/

    const Matrix<T> operator +(const Matrix<T>& mat) const; //矩阵加矩阵
    const Matrix<T> operator +(const T& num) const; //矩阵加标量
    const Matrix<T> operator -(const Matrix<T>& mat) const;
    const Matrix<T> operator -(const T& num) const;
	const Matrix<T> operator -() const;
    const Matrix<T> operator *(const Matrix<T>& mat) const;
    const Vector<T> operator *(const Vector<T>& vec) const;
    const Matrix<T> operator *(const T& num) const;
    friend const Matrix<T> operator *(const T& num,const Matrix<T>& mat)
    {
        return mat*num;
    }



    Matrix<T>& operator +=(const Matrix<T>& mat);
    Matrix<T>& operator +=(const T& num);
    Matrix<T>& operator -=(const Matrix<T>& mat);
    Matrix<T>& operator -=(const T& num);
	Matrix<T>& operator *=(const Matrix<T>& mat);      //mat必须是方阵
    Matrix<T>& operator *=(const T& num);
	Matrix<T>& operator /= (const T& num);

    T& operator ()(size_t index_row,size_t index_col);//将函数操作符重载，实现寻址操作符功能
    const T operator()(size_t index_row,size_t index_col) const;//供常对象使用


    friend ostream& operator<< (ostream& os,const Matrix<T>& mat)//重载输出操作符
    {
        size_t Row = mat.GetNumRow();
        size_t Col = mat.GetNumCol();

        for (size_t i = 0; i < Row; i++)
        {
            for (size_t j = 0; j < Col; j++)
            {
                os << mat(i, j) << " ";
            }
            os << endl;
        }
        return os;
    }
	
	

    /**************
     *基本的矩阵操作*
     **************/
    void SetZeros();   //所有元素置零
    void IdentityMatrix();  // 单位矩阵
    void Resize(size_t Row, size_t Col); //重新分配


	//Matrix<T>& TransPose();
	const Matrix<T> TransPose() const;

    Matrix<T>& FirstTypeTransForm(size_t Row_One, size_t Row_Two);
    const Matrix<T> FirstTypeTransForm(size_t Row_One, size_t Row_Two) const;
	//第一类初等变换---交换两行

    Matrix<T>& SecondTypeTransForm(size_t Row, double Num);
    const Matrix<T> SecondTypeTransForm(size_t Row, double Num) const;
	//第二类初等变换---某一行乘以一个数
	// Row---待乘数的行
	// Num---需要乘以的数

    Matrix<T>& ThirdTypeTransForm(size_t Row_One, size_t Row_Two, double Num);
    const Matrix<T> ThirdTypeTransForm(size_t Row_One, size_t Row_Two, double Num) const;
	//第三类变换----某一行的倍数加到另一行
	//Row_One ----- 需要乘数字的一行
	//Row_Two ----- 待加的一行
	//Num    ------ 乘数


    const Matrix<T> ExtractBlock(size_t RowStart, size_t ColStart, size_t RowNumToBlock, size_t ColNumToBlock) const;
    void SetBlock(size_t RowStart, size_t ColStart, size_t RowNumToSet, size_t ColNumtoSet, const Matrix<T>& mat);

    void FindMax(size_t RowStart, size_t ColStart, size_t RowNumToFind, size_t ColNumToFind,T* Value,size_t* row,size_t* col);
    void FindMin(size_t RowStart, size_t ColStart, size_t RowNumToFind, size_t ColNumToFind,T* Value,size_t* row,size_t* col);
		


    double norm_1() const;
    double norm_2() const;
    double norm_Inf() const;

    inline size_t GetNumRow();
    inline size_t GetNumRow() const;
    inline size_t GetNumCol();
    inline size_t GetNumCol() const;
    inline size_t GetNumData();

	~Matrix();



private:
    size_t NumRow;
    size_t NumCol;
    size_t Size;
	T** p1;

//    Mat PetMat;            //Petsc矩阵表示方法


private:
    void Allocate(size_t NumRow,size_t NumCol);
	void DeAllocate();

	void Swap(Matrix<T>& mat);

//    void BuildPETSCMat();
};




template<typename T>
void Matrix<T>::Allocate(size_t NumRow,size_t NumCol)
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
        assert(p1 != nullptr);
        for (size_t i = 0; i < NumRow; ++i)
        {
            p1[i] = new T[NumCol]; // 指向二维数组每行的开头位置
            assert(p1[i] != nullptr);
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

    for(size_t i = 0;i < NumRow;++i)
    {
        delete [] p1[i];
        p1[i] = nullptr;
    }
    delete [] p1;
    p1 = nullptr;
    NumRow = 0;
    NumCol = 0;
    Size = 0;
}


template<typename T>
void Matrix<T>::Swap(Matrix<T>& mat)
{
	T tmp;
	
	assert(NumCol = mat.NumCol);
	assert(NumRow = mat.NumRow);

    for (size_t i = 0; i < NumRow; i++)
	{
        for (size_t j = 0; j < NumCol; j++)
		{
			tmp = mat.p1[i][j];
			mat.p1[i][j] = p1[i][j];
			p1[i][j] = tmp;
		}
	}
}

template<typename T>
Matrix<T>::Matrix() :NumRow(0), NumCol(0), Size(0), p1(nullptr)
{
    
}

template <typename T>
Matrix<T>::Matrix(size_t Row, size_t Col) :NumRow(Row), NumCol(Col)
{
    Allocate(NumRow,NumCol);
}
template <typename T>
Matrix<T>::Matrix(const Matrix& mat) : NumRow(mat.NumRow), NumCol(mat.NumCol)
{
    Allocate(NumRow,NumCol);
    for (size_t i = 0; i < NumRow; ++i)
    {
        for (size_t j = 0; j < NumCol; ++j)
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
        for (size_t i = 0; i < NumRow; ++i)
        {
            for (size_t j = 0; j < NumCol; ++j)
            {
                p1[i][j] = mat.p1[i][j];
            }
        }
    }


    return *this;
}

template <typename T>
Matrix<T>::Matrix(size_t Row, size_t Col, T** arr):NumRow(Row),NumCol(Col)
{

    Allocate(NumRow,NumCol);  //分配存储空间
    for (size_t i = 0;i < NumRow;++i)
    {
        for(size_t j = 0;j < NumCol;++j)
        {
            p1[i][j] = arr[i][j];
        }
    }


}

template <typename T>
Matrix<T>::Matrix(size_t Row,size_t Col,T * arr,const string& Storage)
    :NumRow(Row),NumCol(Col)
{
    Allocate(NumRow,NumCol);
    if(Storage == "Row")
    {
        for(size_t i= 0;i < Row;i++)
        {
            for(size_t j = 0;j < Col;j++)
            {
                p1[i][j] = arr[Col * i + j];
            }
        }
    }
}


template<typename T>
void Matrix<T>::SetZeros()
{
    for(size_t i = 0;i < NumRow;++i)
    {
        for (size_t j = 0;j < NumCol;++j)
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
    for(size_t i = 0;i < NumRow;i++)
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

    for(size_t i = 0;i < NumRow;++i)
    {
        for(size_t j = 0;j < NumCol;++j)
        {
            res_mat.p1[i][j] = p1[i][j] + mat.p1[i][j];
        }
    }
    return res_mat;
}

template<typename T>
const Matrix<T> Matrix<T>::operator +(const T& num) const
{
    //向量与标量做加法
    Matrix<T> res_mat(NumRow,NumCol);
    for(size_t i = 0;i < NumRow;++i)
    {
        for(size_t j = 0;j < NumCol;++j)
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

    for(size_t i = 0;i < NumRow;++i)
    {
        for(size_t j = 0;j < NumCol;++j)
        {
            res_mat.p1[i][j] = p1[i][j] - mat.p1[i][j];
        }
    }
    return res_mat;
}

template<typename T>
const Matrix<T> Matrix<T>::operator -(const T& num) const
{
    //向量与标量做加法
    Matrix<T> res_mat(NumRow,NumCol);
    for(size_t i = 0;i < NumRow;++i)
    {
        for(size_t j = 0;j < NumCol;++j)
        {
            res_mat.p1[i][j] = p1[i][j] - num;
        }
    }
    return res_mat;
}

template<typename T>
const Matrix<T> Matrix<T>::operator -() const
{
	Matrix<T> res_mat(NumRow, NumCol);
    for (size_t i = 0; i < NumRow; ++i)
	{
        for (size_t j = 0; j < NumCol; ++j)
		{
			res_mat.p1[i][j] = -p1[i][j];
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

    for(size_t i = 0;i < NumRow;++i)
    {
        for(size_t j = 0;j < mat.NumCol;++j)
        {
            for(size_t k = 0;k < NumCol;++k)
            {
                res_mat.p1[i][j] += p1[i][k]*mat.p1[k][j];
            }
        }
    }
    return res_mat;
}

template<typename T>
const Matrix<T> Matrix<T>::operator *(const T& num) const
{
    Matrix<T> res_mat(NumRow,NumCol);
    for(size_t i = 0;i < NumRow;++i)
    {
        for(size_t j = 0;j < NumCol;++j)
        {
            res_mat.p1[i][j] = p1[i][j]*num;
        }
    }
    return res_mat;
}

template <typename T>
const Vector<T> Matrix<T>::operator *(const Vector<T>& vec) const
{
    Vector<T> RC(NumRow);

    if(NumCol != vec.getNum())
    {
        throw runtime_error("dimension is not equal");
    }

    for(size_t i = 0 ; i < NumRow;i++)
    {
        Vector<T> RowVec(NumCol);
        for(size_t j = 0 ; j < NumCol;j++)
        {
            RowVec.newValue(j,p1[i][j]);
        }

        T res = RowVec * vec;
        RC.newValue(i,res);
    }

    return RC;
}

template<typename T>
Matrix<T>& Matrix<T>::operator +=(const Matrix<T>& mat)
{
    if(NumRow != mat.NumRow || NumCol != mat.NumCol)
    {
        cout << "两个矩阵之间的维数不相等，程序退出" << endl;
        exit(0);
    }
    for(size_t i = 0;i < NumRow;++i)
    {
        for(size_t j = 0;j < NumCol;++j)
        {
            p1[i][j] += mat.p1[i][j];
        }
    }
    return *this;
}

template<typename T>
Matrix<T>& Matrix<T>::operator +=(const T& num)
{
    for(size_t i = 0;i < NumRow;++i)
    {
        for(size_t j = 0;j < NumCol;++j)
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
    for(size_t i = 0;i < NumRow;++i)
    {
        for(size_t j = 0;j < NumCol;++j)
        {
            p1[i][j] -= mat.p1[i][j];
        }
    }
    return *this;
}

template<typename T>
Matrix<T>& Matrix<T>::operator -=(const T& num)
{
    for(size_t i = 0;i < NumRow;++i)
    {
        for(size_t j = 0;j < NumCol;++j)
        {
            p1[i][j] -= num;
        }
    }
    return *this;
}

template<typename T>
Matrix<T>& Matrix<T>::operator *=(const Matrix<T>& mat)
{
	Matrix<T> tmp = (*this) * mat;
	Swap(tmp);
	return *this;
}


template<typename T>
Matrix<T>& Matrix<T>::operator *=(const T& num)
{
    for(size_t i = 0;i < NumRow;++i)
    {
        for(size_t j = 0;j < NumCol;++j)
        {
            p1[i][j] *= num;
        }
    }
    return *this;
}

template<typename T>
Matrix<T>& Matrix<T>::operator /= (const T& num)
{
    for (size_t i = 0; i < NumRow; ++i)
	{
        for (size_t j = 0; j < NumCol; ++j)
		{
			p1[i][j] /= num;
		}
	}
	return *this;
}

template<typename T>
T& Matrix<T>::operator ()(size_t index_row,size_t index_col)
{
    if(index_row >= NumRow || index_col >= NumCol)
    {
        throw out_of_range("Out of dimension");
    }


    return p1[index_row][index_col];
}

template<typename T>
const T Matrix<T>::operator ()(size_t index_row,size_t index_col) const
{
    if(index_row >= NumRow || index_col >= NumCol)
    {
        throw out_of_range("Out of dimension");
    }
    return p1[index_row][index_col];
}


template<typename T>
inline size_t Matrix<T>::GetNumRow()
{
    return NumRow;
}

template<typename T>
inline size_t Matrix<T>::GetNumRow() const
{
	return NumRow;
}

template<typename T>
inline size_t Matrix<T>::GetNumCol()
{
    return NumCol;
}

template<typename T>
inline size_t Matrix<T>::GetNumCol() const
{
	return NumCol;
}

template<typename T>
inline size_t Matrix<T>::GetNumData()
{
    return Size;
}

template<typename T>
const Matrix<T> Matrix<T>::ExtractBlock(size_t RowStart, size_t ColStart, size_t RowNumToBlock, size_t ColNumToBlock) const
{
	//抽取一个Block出来

	assert(RowStart + RowNumToBlock <= NumRow);
	assert(ColStart + ColNumToBlock <= NumCol);


	Matrix<T> tmp(RowNumToBlock, ColNumToBlock);
    for (size_t i = RowStart; i < RowStart+ RowNumToBlock; i++)
	{
        for (size_t j = ColStart; j < ColStart + ColNumToBlock; j++)
		{
			tmp.p1[i - RowStart][j - ColStart] = p1[i][j];
		}
	}
	return tmp;
}

template<typename T>
void Matrix<T>::SetBlock(size_t RowStart, size_t ColStart, size_t RowNumToSet, size_t ColNumToSet, const Matrix<T>& mat)
{
	assert(RowStart + RowNumToSet <= NumRow);
	assert(ColStart + ColNumToSet <= NumCol);

	assert(RowNumToSet <= mat.NumRow);
	assert(ColNumToSet <= mat.NumCol);

    for (size_t i = 0; i < RowNumToSet; i++)
	{
        for (size_t j = 0; j < ColNumToSet; j++)
		{
			p1[i + RowStart][j + ColStart] = mat.p1[i][j];
		}
	}
}

template<typename T>
void Matrix<T>::FindMax(size_t RowStart, size_t ColStart, size_t RowNumToFind, size_t ColNumToFind, T* Value, size_t* row, size_t* col)
{
	assert(RowStart + RowNumToFind <= NumRow);
	assert(ColStart + ColNumToFind <= NumCol);

	T max = p1[RowStart][ColStart];
    for (size_t i = RowStart; i < RowStart + RowNumToFind; i++)
	{
        for (size_t j = ColStart; j < ColStart + ColNumToFind; j++)
		{
			if (max < p1[i][j])
			{
				max = p1[i][j];
				*row = i;
				*col = j;
			}
		}
	}

	*Value = max;
}

template<typename T>
void Matrix<T>::FindMin(size_t RowStart, size_t ColStart, size_t RowNumToFind, size_t ColNumToFind, T* Value, size_t* row, size_t* col)
{
	assert(RowStart + RowNumToFind <= NumRow);
	assert(ColStart + ColNumToFind <= NumCol);

	T min = p1[RowStart][ColStart];
    for (size_t i = RowStart; i < RowStart + RowNumToFind; i++)
	{
        for (size_t j = ColStart; j < ColStart + ColNumToFind; j++)
		{
			if (min > p1[i][j])
			{
				min = p1[i][j];
				*row = i;
				*row = j;
			}
		}
	}
	
	*Value = min;
}


template<typename T>
Matrix<T>& Matrix<T>::FirstTypeTransForm(size_t Row_One, size_t Row_Two)
{
	//交换两行

	T * row_tmp = nullptr;
	row_tmp = p1[Row_One];
	p1[Row_One] = p1[Row_Two];
	p1[Row_Two] = row_tmp;

	return *this;
}


template<typename T>
const Matrix<T> Matrix<T>::FirstTypeTransForm(size_t Row_One, size_t Row_Two) const
{
	//交换两行

	Matrix<T> tmpMat = *this;

	T * row_tmp = nullptr;
	row_tmp = tmpMat.p1[Row_One];
	tmpMat.p1[Row_One] = tmpMat.p1[Row_Two];
	tmpMat.p1[Row_Two] = row_tmp;

	return tmpMat;
}



template<typename T>
Matrix<T>& Matrix<T>::SecondTypeTransForm(size_t Row, double Num)
{
//第二类初等变换---某一行乘以一个数
// Row---待乘数的行
// Num---需要乘以的数
	
	T * row_tmp = p1[Row];

    for (size_t i = 0; i < NumCol; i++)
	{
		row_tmp[i] *= Num;
	}
	return *this;
}


template<typename T>
const Matrix<T> Matrix<T>::SecondTypeTransForm(size_t Row, double Num) const
{
	//第二类初等变换---某一行乘以一个数
	// Row---待乘数的行
	// Num---需要乘以的数
	Matrix<T> tmpMat = *this;

	T * row_tmp = tmpMat.p1[Row];

    for (size_t i = 0; i < NumCol; i++)
	{
		row_tmp[i] *= Num;
	}
	return tmpMat;

}



template<typename T>
Matrix<T>& Matrix<T>::ThirdTypeTransForm(size_t Row_One, size_t Row_Two, double Num)
{
	//第三类变换----某一行的倍数加到另一行
	//Row_One ----- 需要乘数字的一行
	//Row_Two ----- 待加的一行
	//Num    ------ 乘数


	T *row_one = p1[Row_One];
	T *row_two = p1[Row_Two];

    for (size_t i = 0; i < NumCol; i++)
	{
		row_two[i] += row_one[i] * Num;
	}
	return *this;
}


template<typename T>
const Matrix<T> Matrix<T>::ThirdTypeTransForm(size_t Row_One, size_t Row_Two, double Num) const
{
	//第三类变换----某一行的倍数加到另一行
	//Row_One ----- 需要乘数字的一行
	//Row_Two ----- 待加的一行
	//Num    ------ 乘数

	Matrix<T> tmpMat = *this;

	T *row_one = tmpMat.p1[Row_One];
	T *row_two = tmpMat.p1[Row_Two];

    for (size_t i = 0; i < NumCol; i++)
	{
		row_two[i] += row_one[i] * Num;
	}
	return tmpMat;
}

template<typename T>
void Matrix<T>::Resize(size_t Row, size_t Col)
{
	if (p1 != nullptr)
	{
		DeAllocate();
	}
	Allocate(Row, Col);
}

template<typename T>
const Matrix<T> Matrix<T>::TransPose() const
{
	Matrix<T> mat(NumCol,NumRow);
	T tmp;
    for (size_t i = 0; i < NumRow; i++)
	{
        for (size_t j = 0; j < NumCol; j++)
		{
			mat.p1[j][i] = p1[i][j];
		}
	}
	return mat;
}

template <typename T>
double Matrix<T>::norm_1() const
{
    vector<T> SUMCol;

    for(size_t i = 0 ; i < NumCol;i++)
    {
        double SUM = 0.0;

        for(size_t j = 0 ; j < NumRow;j++)
        {
            SUM += p1[j][i];
        }

        SUMCol.push_back(SUM);
    }

    return max_element(SUMCol.begin(),SUMCol.end());
}

template<typename T>
double Matrix<T>::norm_2() const
{

}

template<typename T>
double Matrix<T>::norm_Inf() const
{
    vector<T> SUMRow;

    for(size_t i = 0 ; i < NumRow;i++)
    {
        double SUM = 0.0;

        for(size_t j = 0 ; j < NumCol;j++)
        {
            SUM += p1[i][j];
        }

        SUMRow.push_back(SUM);
    }

    return max_element(SUMRow.begin(),SUMRow.end());
}



#endif   // MATRIX_H_
