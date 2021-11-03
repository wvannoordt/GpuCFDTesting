#ifndef REDUCE_H
#define REDUCE_H

#include "FlowField.h"

template <class ReduceOp, class Transform> __global__ void K_Reduce(MdArray<double, 4> arr, const Transform& callable_transform)
{
    
}

template <class ReduceOp, class Transform> void Reduce(FlowField& flow, const Transform& callable_transform)
{
    
}

#endif