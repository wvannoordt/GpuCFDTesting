#ifndef M9_H
#define M9_H
#include "CudaHeaders.h"
template <typename T> struct m9
{
    T data[9];
    
    _f_hybrid m9(void)
    {
        data[0] = 0.0;
        data[1] = 0.0;
        data[2] = 0.0;
        data[3] = 0.0;
        data[4] = 0.0;
        data[5] = 0.0;
        data[6] = 0.0;
        data[7] = 0.0;
        data[8] = 0.0;
    }
    
    _f_hybrid m9(const m9<T> & rhs)
    {
        data[0] = rhs.data[0];
        data[1] = rhs.data[1];
        data[2] = rhs.data[2];
        data[3] = rhs.data[3];
        data[4] = rhs.data[4];
        data[5] = rhs.data[5];
        data[6] = rhs.data[6];
        data[7] = rhs.data[7];
        data[8] = rhs.data[8];
    }
    
    _f_hybrid T& operator () (const int i, const int j) {return data[i+3*j];}
};

#endif