#ifndef VECTOR_H
#define VECTOR_H

//向量类
#include <stdlib.h>
#include <iostream>
#include <cmath>

#include "Matrix.h"

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

    const Matrix<T> ToMatrix(const string& Type);

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
    const Vector<T> operator/(const T& n)const;


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


    //变换
    vector<Matrix<T>> House();    //HouseHolder变换


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

template <typename T>
const Vector<T> Vector<T>::operator / (const T& n) const
{
    return operator*(1/n);
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

template <typename T>
vector<Matrix<T>> Vector<T>::House()
{
    vector<Matrix<T>> RC;

    double InfNorm = this->norm_Inf();
    Matrix<T> mat = this->ToMatrix("Col");

    Matrix<T> muMat(num,1);
    muMat.SetBlock(1,0,muMat.GetNumRow() - 1,muMat.GetNumCol(),mat.ExtractBlock(1,0,num - 1,1));

    Matrix<T> sigmaMat = mat.ExtractBlock(1,0,num - 1,1).TransPose()
            * mat.ExtractBlock(1,0,num - 1,1);

    auto sigma = static_cast<double>(sigmaMat(0,0));


    Matrix<T> betaMat(1,1);
    if(abs(sigma) <= EPS)
    {
        //sigma =0

        betaMat(0,0) = 0;

        RC.push_back(muMat,betaMat);
        return RC;
    }
    else
    {
        auto alpha = static_cast<double>((pow(mat(0,0),2) + sigma));
        if(mat(0,0) < 0 || abs(data[0])<= EPS)
        {
            muMat(0,0) = mat(0,0) - alpha;
        }
        else
        {
            muMat(0,0) = static_cast<T>(-sigma/(mat(0,0) + alpha));
        }

        betaMat(0,0) = 2  / (sigma + pow(muMat(0,0),2));
//        muMat = muMat / muMat(0,0);
        RC.push_back(muMat,betaMat);

        return RC;
    }
}


template <typename T>
const Matrix<T> Vector<T>::ToMatrix(const string& StorageType)
{
    size_t num = this->getNum();
    if(StorageType == "Row")
    {
        Matrix<T> RC(1,num);
        for(int i = 0 ; i < num;i++)
        {
            RC(0,i) = data[i];
        }
        return RC;
    }
    else if(StorageType == "Col")
    {
        Matrix<T> RC(num,1);
        for( int i = 0 ; i < num ; i++)
        {
            RC(i,0) = data[i];
        }
        return RC;
    }
    else
    {
        throw runtime_error("Wrong StorageType");
    }
}

#endif // VECTOR_H
