#ifndef VECTOR_H
#define VECTOR_H

//向量类
#include "storagetype.h"
#include <stdlib.h>
#include <iostream>
#include <cmath>

using namespace std;

template <typename T>
class Vector
{
public:
    Vector();
    Vector(size_t num);
    Vector(const Vector<T >& vec);
    Vector<T >& operator=(const Vector<T >& vec);

    Vector(size_t num,T* arr);

    ~Vector();

    /*
     * 重载操作符号
     */

    const Vector<T> operator+(const Vector<T>& vec) const;
    const Vector<T> operator+(const T n) const;
    const Vector<T> operator-(const Vector<T>& vec) const;
    const Vector<T> operator-(const T n) const;
    const Vector<T> operator-() const;
    const Vector<T> operator*(const T n)const;

    T& operator()(size_t i);
    const T operator()(size_t i) const;

private:
    void Allocate(size_t num);
    void DeAllocate();

private:
    T* data;
    int num;

};

template <typename T   >
void Vector<T >::Allocate(size_t num)
{
    if(data != nullptr)
    {
        delete [] data;
    }
    data = new T[num];
}


template<typename T   >
void Vector<T >::DeAllocate()
{
    if(data != nullptr)
    {
        delete [] data;
        data = nullptr;
    }
}


template <typename T   >
Vector<T >::Vector():num(0),data(nullptr)
{

}

template <typename T   >
Vector<T >::Vector(size_t dataNum):num(dataNum)
{
    Allocate(dataNum);
}

template <typename T   >
Vector<T >::Vector(const Vector<T >& vec):num(vec.num)
{
    Allocate(num);
    for(size_t i = 0 ; i < num;i++)
    {
        data[i] = vec.data[i];
    }
}

template <typename T   >
Vector<T >& Vector<T >::operator=(const Vector<T >& vec)
{
    if(&vec != this)
    {
        DeAllocate();
        Allocate(vec.num);
        for(size_t i = 0 ; i < num;i++)
        {
            data[i] = vec.data[i];
        }
    }

    return *this;
}

template<typename T   >
Vector<T >::~Vector<T >()
{
    DeAllocate();
}

template <typename T   >
Vector<T >::Vector(size_t DataNum,T* arr)
    :num(DataNum)
{
    Allocate(DataNum);

    for(size_t i = 0 ; i < DataNum;i++)
    {
        data[i] = arr[i];
    }
}

template <typename T   >
const Vector<T > Vector<T >::operator+(const Vector<T >& vec) const
{
    if(vec.num != num)
    {
        throw runtime_error("The dimension is not equal");
    }

    Vector<T > RC(num);

    for(size_t i = 0 ; i < num;i++)
    {
        RC.data[i] = data[i] + vec.data[i];
    }

    return RC;
}


template <typename T   >
const Vector<T > Vector<T >::operator+(const T n) const
{
    Vector<T> RC(num);

    for(size_t i = 0 ; i < num;i++)
    {
        RC(i) = data[i] + n;
    }

    return RC;
}

template <typename T>
const Vector<T> Vector<T>::operator-(const Vector<T>& vec) const
{
    Vector<T> RC(num);
    RC  = &this + (-vec);
    return RC;

}

template <typename T>
const Vector<T> Vector<T>::operator-(const T n) const
{
    Vector<T> RC(num);
    RC = &this + (-n);
    return RC;
}

template <typename T>
const Vector<T> Vector<T>::operator-() const
{
    Vector<T> RC(num);
    for(size_t i = 0 ; i < num;i++)
    {
        RC.data[i] = -RC.data[i];
    }

    return RC;
}

template <typename T>
const Vector<T> Vector<T>::operator*(const T n)const
{
    Vector<T> RC(num);

    for(size_t i = 0 ; i < num;i++)
    {
        RC(i) = data[i] * n;
    }

    return RC;
}

template <typename T>
T& Vector<T>::operator()(size_t i)
{
    if(i >= num)
    {
        throw out_of_range("i is greater than dimension");
    }



    return data[i];
}

template <typename T>
const T Vector<T>::operator()(size_t i) const
{
    if(i >= num)
    {
        throw out_of_range("i is greater than dimension");
    }

   return data[i];
}

#endif // VECTOR_H
