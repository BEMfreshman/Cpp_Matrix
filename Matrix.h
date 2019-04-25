#ifndef MATRIX_H_
#define MATRIX_H_

//#include <petscmat.h>

#include <iostream>
#include <assert.h>
#include <cmath>
#include <vector>
#include <algorithm>

#include "NotImplementedExcetion.h"


#define EPS 1e-10

using namespace std;


template <typename T>
class LU;

template <typename T>
class GaussSolver;

template <typename T>
class Utility;

template<typename T>
class Matrix {
public:
    Matrix();

    Matrix(size_t Row, size_t Col);

    Matrix(const Matrix<T> &mat);

    Matrix<T> &operator=(const Matrix<T> &mat);  //赋值操作
    Matrix(size_t Row, size_t Col, T **arr);

    Matrix(size_t Row, size_t Col, T *arr, const string &StorageType);

    //采用数组初始化
    //第一二个参数为本矩阵的行数和列数
    //第三个参数为矩阵初始化的数组
    //第四个为第三个参数中的元素个数

    /***********
     *重载操作符*
     **********/

    const Matrix<T> operator+(const Matrix<T> &mat) const; //矩阵加矩阵
    const Matrix<T> operator+(const T &num) const; //矩阵加标量
    const Matrix<T> operator-(const Matrix<T> &mat) const;

    const Matrix<T> operator-(const T &num) const;

    const Matrix<T> operator-() const;

    const Matrix<T> operator*(const Matrix<T> &mat) const;

    const Matrix<T> operator*(const T &num) const;

    friend const Matrix<T> operator*(const T &num, const Matrix<T> &mat) {
        return mat * num;
    }


    Matrix<T> &operator+=(const Matrix<T> &mat);

    Matrix<T> &operator+=(const T &num);

    Matrix<T> &operator-=(const Matrix<T> &mat);

    Matrix<T> &operator-=(const T &num);

    Matrix<T> &operator*=(const Matrix<T> &mat);      //mat必须是方阵
    Matrix<T> &operator*=(const T &num);

    Matrix<T> &operator/=(const T &num);

    T &operator()(size_t index_row, size_t index_col);//将函数操作符重载，实现寻址操作符功能
    const T operator()(size_t index_row, size_t index_col) const;//供常对象使用
    const T get(size_t index_row, size_t index_col) const;

    friend ostream &operator<<(ostream &os, const Matrix<T> &mat)//重载输出操作符
    {
        size_t Row = mat.GetNumRow();
        size_t Col = mat.GetNumCol();

        for (size_t i = 0; i < Row; i++) {
            for (size_t j = 0; j < Col; j++) {
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
    void SetConstants(const T& value);     //置为常数
    void IdentityMatrix();  // 单位矩阵
    void Resize(size_t Row, size_t Col); //重新分配
    T det();

    const Matrix<T> TransPose() const;
    const Matrix<T> Inv() const;


    Matrix<T> &FirstTypeTransForm(size_t Row_One, size_t Row_Two);

    const Matrix<T> FirstTypeTransForm(size_t Row_One, size_t Row_Two) const;
    //第一类初等变换---交换两行

    Matrix<T> &SecondTypeTransForm(size_t Row, double Num);

    const Matrix<T> SecondTypeTransForm(size_t Row, double Num) const;
    //第二类初等变换---某一行乘以一个数
    // Row---待乘数的行
    // Num---需要乘以的数

    Matrix<T> &ThirdTypeTransForm(size_t Row_One, size_t Row_Two, double Num);

    const Matrix<T> ThirdTypeTransForm(size_t Row_One, size_t Row_Two, double Num) const;
    //第三类变换----某一行的倍数加到另一行
    //Row_One ----- 需要乘数字的一行
    //Row_Two ----- 待加的一行
    //Num    ------ 乘数


    const Matrix<T> ExtractBlock(size_t RowStart, size_t ColStart, size_t RowNumToBlock, size_t ColNumToBlock) const;

    void SetBlock(size_t RowStart, size_t ColStart, size_t RowNumToSet, size_t ColNumtoSet, const Matrix<T> &mat);
    void SetRow(size_t i, const Matrix<T>& RowMat);
    void SetCol(size_t i, const Matrix<T>& ColMat);



    void FindMax(size_t RowStart, size_t ColStart, size_t RowNumToFind, size_t ColNumToFind, T *Value, size_t *row,
                 size_t *col);

    void FindMin(size_t RowStart, size_t ColStart, size_t RowNumToFind, size_t ColNumToFind, T *Value, size_t *row,
                 size_t *col);


    double norm_1() const;

    double norm_2() const;

    double norm_Inf() const;

    double Getcond() const;       //计算条件数

    //矩阵判断
    inline bool isSquare() const {return NumRow == NumCol;};
    inline bool isUpTri() const;
    inline bool isLowTri() const;
    inline bool isDiag() const;
    inline bool isSymmetric() const;
    inline bool HasZerosDiag() const;
    bool operator==(const Matrix<T>& mat) const ;



    inline size_t GetNumRow();

    inline size_t GetNumRow() const;

    inline size_t GetNumCol();

    inline size_t GetNumCol() const;

    inline size_t GetNumData();

    //取子矩阵
    const Matrix<T> GetLowTriMatrix(bool needDiag) const;
    const Matrix<T> GetUpTriMatrix(bool needDiag) const;

    const Matrix<T> GetDiagMatrix() const;
    const Matrix<T> sign() const;


    //取最大值，最小值
    size_t GetMaxIdRow(size_t RowId);
    size_t GetMinIdRow(size_t RowId);

    size_t GetMaxIdCol(size_t ColId);
    size_t GetMinIdCol(size_t ColId);

    void GetMaxId(size_t* RowNo,size_t* ColNo);
    void GetMinId(size_t* RowNo,size_t* ColNo);




    ~Matrix();


    /*
     * 友元类
     */

    friend class LU<T>;
    friend class GaussSolver<T>;
    friend class Utility<T>;


private:
    size_t NumRow;
    size_t NumCol;
    size_t Size;
    T **p1;

//    Mat PetMat;            //Petsc矩阵表示方法


private:
    void Allocate(size_t NumRow, size_t NumCol);

    void DeAllocate();

    void Swap(Matrix<T> &mat);

    const Matrix<T> InvTri() const;

//    void BuildPETSCMat();
};


template<typename T>
void Matrix<T>::Allocate(size_t NumRow, size_t NumCol) {
    //分配内存
    if (NumRow <= 0 || NumCol <= 0) {
        //表示没有设置矩阵的维数，打印信息，不进行分配
        cout << "没有设置矩阵的维数，请检查" << endl;
        exit(0);
    } else {
        Size = NumRow * NumCol;
        p1 = new T *[NumRow];
        assert(p1 != nullptr);
        for (size_t i = 0; i < NumRow; ++i) {
            p1[i] = new T[NumCol]; // 指向二维数组每行的开头位置
            assert(p1[i] != nullptr);
        }
        this->NumRow = NumRow;
        this->NumCol = NumCol;
    }
}

template<typename T>
void Matrix<T>::DeAllocate() {
    //释放矩阵所占有的内存资源
    //并将NumRow，NumCol重置为0

    for (size_t i = 0; i < NumRow; ++i) {
        delete[] p1[i];
        p1[i] = nullptr;
    }
    delete[] p1;
    p1 = nullptr;
    NumRow = 0;
    NumCol = 0;
    Size = 0;
}


template<typename T>
void Matrix<T>::Swap(Matrix<T> &mat) {
    T tmp;

    assert(NumCol = mat.NumCol);
    assert(NumRow = mat.NumRow);

    for (size_t i = 0; i < NumRow; i++) {
        for (size_t j = 0; j < NumCol; j++) {
            tmp = mat.p1[i][j];
            mat.p1[i][j] = p1[i][j];
            p1[i][j] = tmp;
        }
    }
}

template <typename T>
bool Matrix<T>::HasZerosDiag() const
{
    size_t UB = min(GetNumRow(),GetNumCol());

    for(size_t i = 0 ; i < UB;i++)
    {
        if(abs(p1[i][i]) < EPS)
        {
            return true;
        }
    }

    return false;
}

template <typename T>
const Matrix<T> Matrix<T>::InvTri() const
{
    if(isUpTri())
    {
        Matrix<T> RC(GetNumRow(),GetNumCol());

        for(size_t i = 0 ; i < GetNumCol();i++)
        {
            Matrix<T> ei(GetNumRow(),1);
            ei.SetZeros();
            ei(i,0) = 1.0;

            Matrix<T> AnsCol(GetNumRow(),1);
            AnsCol.SetZeros();
            AnsCol = Utility<T>::UptriSolve(*this,ei);
            RC.SetCol(i,AnsCol);
        }

        return RC;
    }
    else if(isLowTri())
    {
        Matrix<T> RC(GetNumRow(),GetNumCol());

        for(size_t i = 0 ; i < GetNumCol();i++)
        {
            Matrix<T> ei(GetNumRow(),1);
            ei.SetZeros();
            ei(i,0) = 1.0;

            Matrix<T> AnsCol(GetNumRow(),1);
            AnsCol.SetZeros();
            AnsCol = Utility<T>::LowtriSolve(*this,ei);
            RC.SetCol(i,AnsCol);
        }

        return RC;
    }
    else
    {
        throw runtime_error("Not Tri Matrix");
    }
}

template<typename T>
Matrix<T>::Matrix() :NumRow(0), NumCol(0), Size(0), p1(nullptr) {

}

template<typename T>
Matrix<T>::Matrix(size_t Row, size_t Col) :NumRow(Row), NumCol(Col) {
    Allocate(NumRow, NumCol);
}

template<typename T>
Matrix<T>::Matrix(const Matrix &mat) : NumRow(mat.NumRow), NumCol(mat.NumCol) {
    Allocate(NumRow, NumCol);
    for (size_t i = 0; i < NumRow; ++i) {
        for (size_t j = 0; j < NumCol; ++j) {
            p1[i][j] = mat.p1[i][j];
        }
    }
}

template<typename T>
Matrix<T> &Matrix<T>::operator=(const Matrix<T> &mat) {
    //需深拷贝
    if (&mat != this) {
        //防止自赋值
        DeAllocate();
        Allocate(mat.NumRow, mat.NumCol);
        for (size_t i = 0; i < NumRow; ++i) {
            for (size_t j = 0; j < NumCol; ++j) {
                p1[i][j] = mat.p1[i][j];
            }
        }
    }


    return *this;
}

template<typename T>
Matrix<T>::Matrix(size_t Row, size_t Col, T **arr):NumRow(Row), NumCol(Col) {

    Allocate(NumRow, NumCol);  //分配存储空间
    for (size_t i = 0; i < NumRow; ++i) {
        for (size_t j = 0; j < NumCol; ++j) {
            p1[i][j] = arr[i][j];
        }
    }


}

template<typename T>
Matrix<T>::Matrix(size_t Row, size_t Col, T *arr, const string &Storage)
        :NumRow(Row), NumCol(Col) {
    Allocate(NumRow, NumCol);
    if (Storage == "Row") {
        for (size_t i = 0; i < Row; i++) {
            for (size_t j = 0; j < Col; j++) {
                p1[i][j] = arr[Col * i + j];
            }
        }
    }
}


template<typename T>
void Matrix<T>::SetZeros() {
    for (size_t i = 0; i < NumRow; ++i) {
        for (size_t j = 0; j < NumCol; ++j) {
            p1[i][j] = T(0);
        }
    }
}

template <typename T>
void Matrix<T>::SetConstants(const T& value)
{
    for (size_t i = 0; i < NumRow; ++i) {
        for (size_t j = 0; j < NumCol; ++j) {
            p1[i][j] = value;
        }
    }
}

template<typename T>
void Matrix<T>::IdentityMatrix() {
    this->SetZeros();
    //首先全部置为0；
    for (size_t i = 0; i < NumRow; i++) {
        p1[i][i] = T(1);
    }
}


template<typename T>
Matrix<T>::~Matrix() {
    DeAllocate();


}

template<typename T>
const Matrix<T> Matrix<T>::operator+(const Matrix<T> &mat) const {
    if (NumRow != mat.NumRow || NumCol != mat.NumCol) {
        cout << "两个矩阵之间的维数不相等，程序退出" << endl;
        exit(0);
    }
    Matrix<T> res_mat(NumRow, NumCol);

    for (size_t i = 0; i < NumRow; ++i) {
        for (size_t j = 0; j < NumCol; ++j) {
            res_mat.p1[i][j] = p1[i][j] + mat.p1[i][j];
        }
    }
    return res_mat;
}

template<typename T>
const Matrix<T> Matrix<T>::operator+(const T &num) const {
    //向量与标量做加法
    Matrix<T> res_mat(NumRow, NumCol);
    for (size_t i = 0; i < NumRow; ++i) {
        for (size_t j = 0; j < NumCol; ++j) {
            res_mat.p1[i][j] = p1[i][j] + num;
        }
    }
    return res_mat;
}

template<typename T>
const Matrix<T> Matrix<T>::operator-(const Matrix<T> &mat) const {
    if (NumRow != mat.NumRow || NumCol != mat.NumCol) {
        cout << "两个矩阵之间的维数不相等，程序退出" << endl;
        exit(0);
    }
    Matrix<T> res_mat(NumRow, NumCol);

    for (size_t i = 0; i < NumRow; ++i) {
        for (size_t j = 0; j < NumCol; ++j) {
            res_mat.p1[i][j] = p1[i][j] - mat.p1[i][j];
        }
    }
    return res_mat;
}

template<typename T>
const Matrix<T> Matrix<T>::operator-(const T &num) const {
    //向量与标量做加法
    Matrix<T> res_mat(NumRow, NumCol);
    for (size_t i = 0; i < NumRow; ++i) {
        for (size_t j = 0; j < NumCol; ++j) {
            res_mat.p1[i][j] = p1[i][j] - num;
        }
    }
    return res_mat;
}

template<typename T>
const Matrix<T> Matrix<T>::operator-() const {
    Matrix<T> res_mat(NumRow, NumCol);
    for (size_t i = 0; i < NumRow; ++i) {
        for (size_t j = 0; j < NumCol; ++j) {
            res_mat.p1[i][j] = -p1[i][j];
        }
    }
    return res_mat;
}

template<typename T>
const Matrix<T> Matrix<T>::operator*(const Matrix<T> &mat) const {
    if (NumCol != mat.NumRow) {
        cout << "第一个矩阵的列数与第二个矩阵的行数不相等，程序退出" << endl;

        throw runtime_error("The number of col of first matrix is not equal to the number of row of second matrix");
    }

    Matrix<T> res_mat(NumRow, mat.NumCol);
    res_mat.SetZeros();

    for (size_t i = 0; i < NumRow; ++i) {
        for (size_t j = 0; j < mat.NumCol; ++j) {
            for (size_t k = 0; k < NumCol; ++k) {
                res_mat.p1[i][j] += p1[i][k] * mat.p1[k][j];
            }
        }
    }
    return res_mat;
}

template<typename T>
const Matrix<T> Matrix<T>::operator*(const T &num) const {
    Matrix<T> res_mat(NumRow, NumCol);
    for (size_t i = 0; i < NumRow; ++i) {
        for (size_t j = 0; j < NumCol; ++j) {
            res_mat.p1[i][j] = p1[i][j] * num;
        }
    }
    return res_mat;
}

template<typename T>
Matrix<T> &Matrix<T>::operator+=(const Matrix<T> &mat) {
    if (NumRow != mat.NumRow || NumCol != mat.NumCol) {
        cout << "两个矩阵之间的维数不相等，程序退出" << endl;
        exit(0);
    }
    for (size_t i = 0; i < NumRow; ++i) {
        for (size_t j = 0; j < NumCol; ++j) {
            p1[i][j] += mat.p1[i][j];
        }
    }
    return *this;
}

template<typename T>
Matrix<T> &Matrix<T>::operator+=(const T &num) {
    for (size_t i = 0; i < NumRow; ++i) {
        for (size_t j = 0; j < NumCol; ++j) {
            p1[i][j] += num;
        }
    }
    return *this;
}

template<typename T>
Matrix<T> &Matrix<T>::operator-=(const Matrix<T> &mat) {
    if (NumRow != mat.NumRow || NumCol != mat.NumCol) {
        cout << "两个矩阵之间的维数不相等，程序退出" << endl;
        exit(0);
    }
    for (size_t i = 0; i < NumRow; ++i) {
        for (size_t j = 0; j < NumCol; ++j) {
            p1[i][j] -= mat.p1[i][j];
        }
    }
    return *this;
}

template<typename T>
Matrix<T> &Matrix<T>::operator-=(const T &num) {
    for (size_t i = 0; i < NumRow; ++i) {
        for (size_t j = 0; j < NumCol; ++j) {
            p1[i][j] -= num;
        }
    }
    return *this;
}

template<typename T>
Matrix<T> &Matrix<T>::operator*=(const Matrix<T> &mat) {
    Matrix<T> tmp = (*this) * mat;
    Swap(tmp);
    return *this;
}


template<typename T>
Matrix<T> &Matrix<T>::operator*=(const T &num) {
    for (size_t i = 0; i < NumRow; ++i) {
        for (size_t j = 0; j < NumCol; ++j) {
            p1[i][j] *= num;
        }
    }
    return *this;
}

template<typename T>
Matrix<T> &Matrix<T>::operator/=(const T &num) {
    for (size_t i = 0; i < NumRow; ++i) {
        for (size_t j = 0; j < NumCol; ++j) {
            p1[i][j] /= num;
        }
    }
    return *this;
}

template<typename T>
T& Matrix<T>::operator()(size_t index_row, size_t index_col) {
    if (index_row >= NumRow || index_col >= NumCol) {
        throw out_of_range("Out of dimension");
    }


    return p1[index_row][index_col];
}

template<typename T>
const T Matrix<T>::operator()(size_t index_row, size_t index_col) const {
    if (index_row >= NumRow || index_col >= NumCol) {
        throw out_of_range("Out of dimension");
    }
    return p1[index_row][index_col];
}

template<typename T>
const T Matrix<T>::get(size_t index_row, size_t index_col) const {
    return operator()(index_row, index_col);
}


template<typename T>
inline size_t Matrix<T>::GetNumRow() {
    return NumRow;
}

template<typename T>
inline size_t Matrix<T>::GetNumRow() const {
    return NumRow;
}

template<typename T>
inline size_t Matrix<T>::GetNumCol() {
    return NumCol;
}

template<typename T>
inline size_t Matrix<T>::GetNumCol() const {
    return NumCol;
}

template<typename T>
inline size_t Matrix<T>::GetNumData() {
    return Size;
}

template<typename T>
const Matrix<T>
Matrix<T>::ExtractBlock(size_t RowStart, size_t ColStart, size_t RowNumToBlock, size_t ColNumToBlock) const {
    //抽取一个Block出来

    assert(RowStart + RowNumToBlock <= NumRow);
    assert(ColStart + ColNumToBlock <= NumCol);


    Matrix<T> tmp(RowNumToBlock, ColNumToBlock);
    for (size_t i = RowStart; i < RowStart + RowNumToBlock; i++) {
        for (size_t j = ColStart; j < ColStart + ColNumToBlock; j++) {
            tmp.p1[i - RowStart][j - ColStart] = p1[i][j];
        }
    }
    return tmp;
}

template<typename T>
void
Matrix<T>::SetBlock(size_t RowStart, size_t ColStart, size_t RowNumToSet, size_t ColNumToSet, const Matrix<T> &mat) {
    assert(RowStart + RowNumToSet <= NumRow);
    assert(ColStart + ColNumToSet <= NumCol);

    assert(RowNumToSet <= mat.NumRow);
    assert(ColNumToSet <= mat.NumCol);

    for (size_t i = 0; i < RowNumToSet; i++) {
        for (size_t j = 0; j < ColNumToSet; j++) {
            p1[i + RowStart][j + ColStart] = mat.p1[i][j];
        }
    }
}

template <typename T>
void Matrix<T>::SetRow(size_t index, const Matrix<T>& RowMat)
{
    if(RowMat.GetNumRow() != 1)
    {
        throw runtime_error("The parameter is not Row Matrix");
    }

    if(index >= GetNumRow())
    {
        throw out_of_range("index is out of range");
    }

    if(RowMat.GetNumCol() != GetNumCol())
    {
        throw runtime_error("The number of col of parameter matrix is not equal to"
                            "this matrix");
    }

    for(size_t i = 0 ; i < RowMat.GetNumCol();i++)
    {
        p1[index][i] = RowMat(0,i);
    }
}

template <typename T>
void Matrix<T>::SetCol(size_t index, const Matrix<T>& ColMat)
{
    if(ColMat.GetNumCol() != 1)
    {
        throw runtime_error("The parameter is not Col Matrix");
    }

    if(index >= GetNumCol())
    {
        throw out_of_range("index is out of range");
    }

    if(ColMat.GetNumRow() != GetNumRow())
    {
        throw runtime_error("The number of row of parameter matrix is not equal to"
                            "this matrix");
    }

    for(size_t i = 0 ; i < ColMat.GetNumRow();i++)
    {
        p1[i][index] = ColMat(i,0);
    }
}

template<typename T>
void
Matrix<T>::FindMax(size_t RowStart, size_t ColStart, size_t RowNumToFind, size_t ColNumToFind, T *Value, size_t *row,
                   size_t *col) {
    assert(RowStart + RowNumToFind <= NumRow);
    assert(ColStart + ColNumToFind <= NumCol);

    T max = p1[RowStart][ColStart];
    for (size_t i = RowStart; i < RowStart + RowNumToFind; i++) {
        for (size_t j = ColStart; j < ColStart + ColNumToFind; j++) {
            if (max < p1[i][j]) {
                max = p1[i][j];
                *row = i;
                *col = j;
            }
        }
    }

    *Value = max;
}

template<typename T>
void
Matrix<T>::FindMin(size_t RowStart, size_t ColStart, size_t RowNumToFind, size_t ColNumToFind, T *Value, size_t *row,
                   size_t *col) {
    assert(RowStart + RowNumToFind <= NumRow);
    assert(ColStart + ColNumToFind <= NumCol);

    T min = p1[RowStart][ColStart];
    for (size_t i = RowStart; i < RowStart + RowNumToFind; i++) {
        for (size_t j = ColStart; j < ColStart + ColNumToFind; j++) {
            if (min > p1[i][j]) {
                min = p1[i][j];
                *row = i;
                *row = j;
            }
        }
    }

    *Value = min;
}


template<typename T>
Matrix<T> &Matrix<T>::FirstTypeTransForm(size_t Row_One, size_t Row_Two) {
    //交换两行

    T *row_tmp = nullptr;
    row_tmp = p1[Row_One];
    p1[Row_One] = p1[Row_Two];
    p1[Row_Two] = row_tmp;

    return *this;
}


template<typename T>
const Matrix<T> Matrix<T>::FirstTypeTransForm(size_t Row_One, size_t Row_Two) const {
    //交换两行

    Matrix<T> tmpMat = *this;

    T *row_tmp = nullptr;
    row_tmp = tmpMat.p1[Row_One];
    tmpMat.p1[Row_One] = tmpMat.p1[Row_Two];
    tmpMat.p1[Row_Two] = row_tmp;

    return tmpMat;
}


template<typename T>
Matrix<T> &Matrix<T>::SecondTypeTransForm(size_t Row, double Num) {
//第二类初等变换---某一行乘以一个数
// Row---待乘数的行
// Num---需要乘以的数

    T *row_tmp = p1[Row];

    for (size_t i = 0; i < NumCol; i++) {
        row_tmp[i] *= Num;
    }
    return *this;
}


template<typename T>
const Matrix<T> Matrix<T>::SecondTypeTransForm(size_t Row, double Num) const {
    //第二类初等变换---某一行乘以一个数
    // Row---待乘数的行
    // Num---需要乘以的数
    Matrix<T> tmpMat = *this;

    T *row_tmp = tmpMat.p1[Row];

    for (size_t i = 0; i < NumCol; i++) {
        row_tmp[i] *= Num;
    }
    return tmpMat;

}


template<typename T>
Matrix<T> &Matrix<T>::ThirdTypeTransForm(size_t Row_One, size_t Row_Two, double Num) {
    //第三类变换----某一行的倍数加到另一行
    //Row_One ----- 需要乘数字的一行
    //Row_Two ----- 待加的一行
    //Num    ------ 乘数


    T *row_one = p1[Row_One];
    T *row_two = p1[Row_Two];

    for (size_t i = 0; i < NumCol; i++) {
        row_two[i] += row_one[i] * Num;
    }
    return *this;
}


template<typename T>
const Matrix<T> Matrix<T>::ThirdTypeTransForm(size_t Row_One, size_t Row_Two, double Num) const {
    //第三类变换----某一行的倍数加到另一行
    //Row_One ----- 需要乘数字的一行
    //Row_Two ----- 待加的一行
    //Num    ------ 乘数

    Matrix<T> tmpMat = *this;

    T *row_one = tmpMat.p1[Row_One];
    T *row_two = tmpMat.p1[Row_Two];

    for (size_t i = 0; i < NumCol; i++) {
        row_two[i] += row_one[i] * Num;
    }
    return tmpMat;
}

template<typename T>
void Matrix<T>::Resize(size_t Row, size_t Col) {
    if (p1 != nullptr) {
        DeAllocate();
    }
    Allocate(Row, Col);
}

template <typename T>
T Matrix<T>::det()
{
    LU<T> ludes(*this);

    //求解行列式
    size_t RowNum = GetNumRow();
    size_t ColNum = GetNumCol();

    if (RowNum != ColNum)
    {
        cout << "错误，非方阵，不可计算行列式" << endl;
        exit(1);
    }
    else if (RowNum == 1)
    {
        //只有一个元素
        return static_cast<T>(*this(0, 0));
    }

    vector<Matrix<T>> LUMatrix = ludes.LUDeCompose();
    Matrix<T> U = LUMatrix[1];

    T Val = U(0,0);

    for (int i = 1; i < RowNum; i++)
    {
        Val *= U(i, i);
    }

    if (ludes.getFirstTransFormTimes() % 2 != 0)
    {
        //做了奇数次第一类变换
        Val = -Val;
    }

    return Val;
}

template<typename T>
const Matrix<T> Matrix<T>::TransPose() const {
    Matrix<T> mat(NumCol, NumRow);
    T tmp;
    for (size_t i = 0; i < NumRow; i++) {
        for (size_t j = 0; j < NumCol; j++) {
            mat.p1[j][i] = p1[i][j];
        }
    }
    return mat;
}

template <typename T>
const Matrix<T> Matrix<T>::Inv() const
{
    //矩阵求逆


    //一般矩阵采用LU分解法求逆
    if(!isSquare())
    {
        throw runtime_error("This Matrix is not Square");
    }


    if(isDiag())
    {
        if(HasZerosDiag())
        {
            throw runtime_error("There is zero Item on Diag in Diag Matrix");
        }

        Matrix<T> RC(GetNumRow(),GetNumCol());
        RC.SetZeros();
        for(size_t i = 0; i < GetNumRow(); i++)
        {
            RC(i,i) = 1.0 / p1[i][i];
        }

        return RC;
    }
    else if(isLowTri())
    {
        return InvTri();
    }
    else if(isUpTri())
    {
        return InvTri();
    }
    else
    {
        Matrix<T> Imat(GetNumRow(),GetNumCol());
        Imat.IdentityMatrix();

        LU<T> lu(*this);
        vector<Matrix<T>> PLUmat = lu.PLUDeCompose();

        Matrix<T> P = PLUmat[0];
        Matrix<T> L = PLUmat[1];
        Matrix<T> U = PLUmat[2];

        Matrix<T> invL = L.Inv();
        Matrix<T> invU = U.Inv();

        Matrix<T> RC(GetNumRow(),GetNumCol());


        RC = invU * invL * P * Imat;

        return RC;
    }
}

template<typename T>
double Matrix<T>::norm_1() const {
    vector<T> SUMCol;

    for (size_t i = 0; i < NumCol; i++) {
        double SUM = 0.0;

        for (size_t j = 0; j < NumRow; j++) {
            SUM += abs(p1[j][i]);
        }

        SUMCol.push_back(SUM);
    }

    return *(max_element(SUMCol.begin(), SUMCol.end()));
}

template<typename T>
double Matrix<T>::norm_2() const
{
    throw NotImplementedExcetion("2-norm of matrix has not been implemented yet");
}

template<typename T>
double Matrix<T>::norm_Inf() const
{
    vector<T> SUMRow;

    for (size_t i = 0; i < NumRow; i++)
    {
        double SUM = 0.0;

        for (size_t j = 0; j < NumCol; j++)
        {
            SUM += abs(p1[i][j]);
        }

        SUMRow.push_back(SUM);
    }

    return *(max_element(SUMRow.begin(), SUMRow.end()));
}

template <typename T>
inline bool Matrix<T>::isUpTri() const
{
    bool isUpTri = true;
    for(size_t i = 0;i < NumRow;i++)
    {
        for(size_t j = 0;j < NumCol;j++)
        {
            if(i > j)
            {
                isUpTri = (abs(p1[i][j]) < EPS);

                if(!isUpTri)
                {
                    return isUpTri;
                }
            }
        }
    }

    return isUpTri;
}



template <typename T>
inline bool Matrix<T>::isLowTri() const
{
    bool isLowTri = true;
    for(size_t i = 0;i < NumRow;i++)
    {
        for(size_t j = 0;j < NumCol;j++)
        {
            if(i < j)
            {
                isLowTri = (abs(p1[i][j]) < EPS);

                if(!isLowTri)
                {
                    return isLowTri;
                }
            }
        }
    }

    return isLowTri;
}

template <typename T>
inline bool Matrix<T>::isDiag() const
{
    return isLowTri() && isUpTri();
}

template <typename T>
inline bool Matrix<T>::isSymmetric() const
{
    return (*this) == (*this).TransPose();
}

template <typename T>
bool Matrix<T>::operator==(const Matrix<T>& mat) const
{
    if(NumRow != mat.NumRow || NumCol != mat.NumCol)
    {
        return false;
    }

    for(size_t i = 0 ; i < NumRow;i++)
    {
        for(size_t j = 0 ; j < NumCol;j++)
        {
            if(abs(p1[i][j] - mat(i,j)) > EPS)
            {
                return false;
            }
        }
    }

    return true;
}


template <typename T>
size_t Matrix<T>::GetMaxIdRow(size_t RowId)
{
    if(RowId > NumRow)
    {
        throw out_of_range("RowId is greater than RowNum");
    }

    T Value = p1[RowId][0];
    size_t index = 0;
    for(size_t i = 1 ; i < NumCol;i++)
    {
        if(Value  < p1[RowId][i])
        {
            index = i;
        }
    }

    return index;
}

template <typename T>
size_t Matrix<T>::GetMinIdRow(size_t RowId)
{
    if(RowId > NumRow)
    {
        throw out_of_range("RowId is greater than RowNum");
    }

    T Value = p1[RowId][0];
    size_t index = 0;
    for(size_t i = 1 ; i < NumCol;i++)
    {
        if(Value > p1[RowId][i])
        {
            index = i;
        }
    }

    return index;
}

template <typename T>
size_t Matrix<T>::GetMaxIdCol(size_t ColId)
{
    if(ColId > NumCol)
    {
        throw out_of_range("ColId is greater than ColNum");
    }

    T Value = p1[0][ColId];
    size_t index = 0;
    for(size_t i = 1 ; i < NumRow;i++)
    {
        if(Value  < p1[i][ColId])
        {
            index = i;
        }
    }

    return index;
}

template <typename T>
size_t Matrix<T>::GetMinIdCol(size_t ColId)
{
    if(ColId > NumCol)
    {
        throw out_of_range("ColId is greater than ColNum");
    }

    T Value = p1[0][ColId];
    size_t index = 0;
    for(size_t i = 1 ; i < NumRow;i++)
    {
        if(Value  > p1[0][i])
        {
            index = i;
        }
    }

    return index;
}

template <typename T>
void Matrix<T>::GetMaxId(size_t* RowNo,size_t* ColNo)
{
    T Value = p1[0][0];
    *RowNo = 0;
    *ColNo = 0;

    for(size_t i = 0 ; i < NumRow;i++)
    {
        for(size_t j = 0 ;j < NumCol;j++)
        {
            if(i == 0 && j == 0)
            {
                continue;
            }
            if(Value < p1[i][j])
            {
                Value = p1[i][j];
                *RowNo = i;
                *ColNo = j;
            }
        }
    }
}

template <typename T>
void Matrix<T>::GetMinId(size_t* RowNo,size_t* ColNo)
{
    T Value = p1[0][0];
    *RowNo = 0;
    *ColNo = 0;

    for(size_t i = 0 ; i < NumRow;i++)
    {
        for(size_t j = 0 ;j < NumCol;j++)
        {
            if(i == 0 && j == 0)
            {
                continue;
            }
            if(Value > p1[i][j])
            {
                Value = p1[i][j];
                *RowNo = i;
                *ColNo = j;
            }
        }
    }
}

template <typename T>
const Matrix<T> Matrix<T>::GetLowTriMatrix(bool needDiag) const
{
    if (!isSquare()) {
        throw runtime_error("This Matrix is not square matrix");
    }

    Matrix<T> RC(NumRow, NumCol);
    RC.SetZeros();

    if (needDiag)
    {
        for (size_t i = 0; i < NumRow; i++)
        {
            for (size_t j = 0; j < NumCol; j++)
            {

                if (i <= j)
                {
                    RC(i, j) = p1[i][j];
                }
            }

        }
    }
    else
    {
        for (size_t i = 0; i < NumRow; i++)
        {
            for (size_t j = 0; j < NumCol; j++)
            {

                if (i < j)
                {
                    RC(i, j) = p1[i][j];
                }
            }

        }
    }

    return RC;
}

template <typename T>
const Matrix<T> Matrix<T>::GetUpTriMatrix(bool needDiag) const
{
    if (!isSquare()) {
        throw runtime_error("This Matrix is not square matrix");
    }

    Matrix<T> RC(NumRow, NumCol);
    RC.SetZeros();

    if (needDiag)
    {
        for (size_t i = 0; i < NumRow; i++)
        {
            for (size_t j = 0; j < NumCol; j++)
            {

                if (i >= j)
                {
                    RC(i, j) = p1[i][j];
                }
            }

        }
    }
    else
    {
        for (size_t i = 0; i < NumRow; i++)
        {
            for (size_t j = 0; j < NumCol; j++)
            {

                if (i > j)
                {
                    RC(i, j) = p1[i][j];
                }
            }

        }
    }

    return RC;
}

template <typename T>
const Matrix<T> Matrix<T>::GetDiagMatrix() const
{
    if(!isSquare())
    {
        throw runtime_error("This Matrix is not square");
    }

    Matrix<T> RC(NumRow, NumCol);
    RC.SetZeros();

    for(size_t i = 0 ; i < NumRow;i++)
    {
        RC(i,i) = p1[i][i];
    }

    return RC;
}

template <typename T>
const Matrix<T> Matrix<T>::sign() const
{
    Matrix<T> RC(NumRow,NumCol);

    for(size_t i = 0 ; i < NumRow;i++)
    {
        for(size_t j = 0 ; j < NumCol;j++)
        {
            RC(i,j) = static_cast<T>(Utility<T>::sgn(p1[i][j]));
        }
    }

    return RC;
}

template <typename T>
double Matrix<T>::Getcond() const
{
    //条件数的计算公式：RES = norm(inv(A)) * norm(A);
    //这个算法重点在于如何计算norm(inv(A))
    //算法实践方案采用了原始文献（见下一个条目）中的practice algorithm中的算法描述，该方案也是LAPACK的实践方案
    //算法理论参见 《数值线性代数》 第二版 P69
    // 原始文献：FORTRAN codes for estimating the one-norm of a real or complex matrix, with applications to condition estimation
    //该算法在LAPACK中有运用，参见https://www.netlib.org/lapack/lug/node38.html中xyycon条目的描述



    int kflag = 1;

    Matrix<T> b(NumRow,1);
    Matrix<T> nu(NumRow,1);


    b.SetConstants(static_cast<T>(1.0/NumRow));
    GaussSolver<T> GS(*this,b);

    try
    {
        nu = GS.ColPivotSolve();
        cout << "nu" << endl;
        cout << nu << endl;
    }
    catch (runtime_error& e)
    {
        throw runtime_error("GS solver for w failed");
    }


    if(NumRow == 1)
    {
        //快速返回
        return abs(static_cast<double>(nu(0,0)));
    }

    double invnorm = nu.norm_1();
    Matrix<T> xi = nu.sign();


    GaussSolver<T> GS_Another(this->TransPose(),xi);

    Matrix<T> x(NumRow,1);
    try
    {
        x = GS_Another.ColPivotSolve();
        cout << "x" << endl;
        cout << x << endl;
    }
    catch (runtime_error& e)
    {
        throw runtime_error("GS solver for z failed");
    }



    int k  = 2;
    while (k <= 5)
    {
        size_t j = x.GetMaxIdCol(0);
        b.SetZeros();
        b(j,0) = 1.0;
        GaussSolver<T> InternalGS1(*this,b);
        try
        {
            nu = InternalGS1.Solve();
            cout << "nu" << endl;
            cout << nu << endl;

        }
        catch (runtime_error& e)
        {
            throw runtime_error("Condition Estimate Internal GS solver 1 failed");
        }

        double invnormOld = invnorm;
        invnorm = nu.norm_1();

        if(invnorm <= invnormOld || xi == nu.sign())
        {
            //表示该迭代步的计算结果未大于上一次的计算结果，或者已经开始死循环
            //准备退出

            b.SetZeros();
            for(size_t i = 0; i < b.GetNumRow();i++)
            {
                size_t indexBaseOne = i + 1;
                b(i, 0) = pow(-1, indexBaseOne + 1) * (1 + (indexBaseOne - 1.0) / (b.GetNumRow() - 1));
            }

            GaussSolver<T> FinalGS(*this,b);
            try
            {
                x = FinalGS.ColPivotSolve();
                cout << "x" << endl;
                cout << x << endl;
            }
            catch(runtime_error& e)
            {
                throw runtime_error("Condition Estimate Final GS solver failed");
            }

            if( 2.0 * x.norm_1() / ( 3 * x.GetNumRow()) > invnorm)
            {
                nu = x;
                invnorm = 2.0 * x.norm_1() / ( 3 * x.GetNumRow());
                break;
            }

        }

        xi = nu.sign();
        GaussSolver<T> InternalGS2(this->TransPose(),xi);
        try
        {
            x = InternalGS2.ColPivotSolve();
        }
        catch (runtime_error& e)
        {
            throw runtime_error("Condition Estimate Internal GS solver 2 failed");
        }
        k = k + 1;
    }

    double norm = this->norm_1();
    return norm * invnorm;

}



#endif   // MATRIX_H_
