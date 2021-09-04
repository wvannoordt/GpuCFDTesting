#ifndef V3_H
#define V3_H
#include "CudaHeaders.h"
template <typename T> struct v3
{
    T data[3];
    
    _f_hybrid v3(T t1, T t2, T t3)
    {
        data[0] = t1;
        data[1] = t2;
        data[2] = t3;
    }
    
    _f_hybrid v3(const v3<T> & rhs)
    {
        data[0] = rhs.data[0];
        data[1] = rhs.data[1];
        data[2] = rhs.data[2];
    }
    
    _f_hybrid T& operator [] (const int i) {return data[i];}
};

#endif