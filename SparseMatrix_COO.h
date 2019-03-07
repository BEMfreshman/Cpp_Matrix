//
// Created by Yangwenlu on 2019/3/7.
//

#ifndef CPP_MATRIX_SPARSEMATRIX_COO_H
#define CPP_MATRIX_SPARSEMATRIX_COO_H

#include <iostream>
#include <cmath>
#include <vector>

//COO格式的SparseMatrix实现
//RowIndex记录行号
//ColIndex记录列号码
//

using namespace std;

template<typename T>
class SparseMatrix_COO
{
public:
    SparseMatrix_COO(size_t MaxRowNum,size_t MaxColNum);

    void insertValue(size_t Row,size_t Col,const T Value);

    void resize(size_t MaxRowNum,size_t MaxColNum);

private:
    size_t nonZeros;
    size_t MaxRowNum;
    size_t MaxColNum;

    vector<size_t> RowIndex;
    vector<size_t> ColIndex;
    vector<T> data;
};

template<typename T>
SparseMatrix_COO<T>::SparseMatrix_COO(size_t maxRowNum,size_t maxColNum)
        :nonZeros(0),MaxRowNum(maxRowNum),MaxColNum(maxColNum)
{

}

template<typename T>
void SparseMatrix_COO<T>::insertValue(size_t Row, size_t Col, const T Value)
{

}

template <typename T>
void SparseMatrix_COO<T>::resize(size_t MaxRowNum, size_t MaxColNum)
{

}


#endif //CPP_MATRIX_SPARSEMATRIX_COO_H
