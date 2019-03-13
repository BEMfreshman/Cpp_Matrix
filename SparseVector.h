//
// Created by Yangwenlu on 2019/3/7.
//

#ifndef CPP_MATRIX_SPARSEVECTOR_H
#define CPP_MATRIX_SPARSEVECTOR_H

#include <vector>
#include <algorithm>
#include <cmath>

using namespace std;

template<typename T>
class SparseVector
{
public:
    explicit SparseVector(size_t maxNum);

    void resize(size_t maxNum);
    void insertValue(size_t index,const T& value);

private:

    size_t maxNum;

    vector<size_t> indexes;
    vector<T> data;
};

template <typename T>
SparseVector<T>::SparseVector(size_t max_Num)
    :maxNum(max_Num)
{

}


template <typename T>
void SparseVector<T>::insertValue(size_t index, const T& value)
{
    vector<size_t>::iterator it = find(indexes.begin(),indexes.end(),index);
    if(it == indexes.end())
    {
        indexes.push_back(index);
        data.push_back(value);
    } else
    {
        data[it-indexes.begin()] = value;
    }
}

template <typename T>
void SparseVector<T>::resize(size_t maxNum)
{
    this->maxNum = maxNum;

    vector<int>().swap(indexes);
    vector<T>().swap(data);
}



#endif //CPP_MATRIX_SPARSEVECTOR_H
