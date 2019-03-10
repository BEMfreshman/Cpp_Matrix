#ifndef VECTOR_H
#define VECTOR_H

//向量类
#include <stdlib.h>
#include <iostream>
#include <cmath>

using namespace std;

template <typename T>
class Vector
{
public:
    Vector();
    explicit Vector(size_t num);
    Vector(const Vector<T >& vec);
    Vector<T >& operator=(const Vector<T >& vec);

    Vector(size_t num,T* arr);

    ~Vector();

    /*
     * 重载操作符号
     */

    const Vector<T> operator+(const Vector<T>& vec) const;
    const Vector<T> operator+(const T& n) const;
    const Vector<T> operator-(const Vector<T>& vec) const;
    const Vector<T> operator-(const T& n) const;
    const Vector<T> operator-() const;
    const Vector<T> operator*(const T& n)const;

    T& operator()(size_t i);
    const T operator()(size_t i) const;

    const T dot(const Vector<T>& vec) const;

    friend ostream& operator<<(ostream& os,const Vector<T>& vec)
    {
        size_t Num = vec.num;
        for(size_t i = 0 ; i < Num;i++)
        {
            os << vec(i);
        }

        os << endl;
        return os;
    }

    size_t getNum() const;


    void newValue(size_t i,const T& value);
    void insertValue(size_t i, const T& value);


    //norm
    double norm_1() const;     //1范数
    double norm_2() const;    //2范数
    double norm_Inf() const;  //Inf范数

private:
    void Allocate(size_t num);
    void DeAllocate();

private:
    T* data;
    size_t num;

};

template <typename T   >
void Vector<T >::Allocate(size_t num)
{
    if(data != nullptr)
    {
        delete [] data;
    }
    data = new T[num];
    for(size_t i = 0 ; i < num;i++)
    {
        data[i] = (T)0;
    }
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
const Vector<T > Vector<T >::operator+(const T& n) const
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
const Vector<T> Vector<T>::operator-(const T& n) const
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
const Vector<T> Vector<T>::operator*(const T& n)const
{
    Vector<T> RC(num);

    for(size_t i = 0 ; i < num;i++)
    {
        RC(i) = data[i] * n;
    }

    return RC;
}
template<typename T>
const T Vector<T>::dot(const Vector<T>& vec) const
{
    T RC = (T)0;

    if(num != vec.num)
    {
        throw runtime_error("dimension is not equal");
    }

    for(size_t i = 0 ; i < num;i++)
    {
        RC += data[i]*vec(i);
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

template <typename T>
double Vector<T>::norm_1() const//1范数
{
    double RC = 0.0;
    for(size_t i = 0 ; i < num;i++)
    {
        RC += abs(data[i]);
    }

    return RC;

}

template <typename T>
double Vector<T>::norm_2() const//2范数
{
    double SUM = 0.0;
    for(size_t i = 0 ; i < num;i++)
    {
        SUM += pow(data[i],2);
    }

    return sqrt(SUM);
}

template <typename T>
double Vector<T>::norm_Inf() const//Inf范数
{
    double maxRC = 0.0;
    for(size_t i = 0 ; i < num ; i++)
    {
        maxRC = max(maxRC,data[i]);
    }

    return maxRC;
}

template <typename T>
size_t Vector<T>::getNum() const
{
    return num;
}

template <typename T>
void Vector<T>::newValue(size_t i,const T& value)
{
    if(i >= num)
    {
        throw out_of_range("i is more than num");
    }

    data[i] = value;
}

template <typename T>
void Vector<T>::insertValue(size_t i, const T& value)
{
    if(i >= num)
    {
        throw out_of_range("i is more than num");
    }

    data[i] += value;
}

#endif // VECTOR_H
