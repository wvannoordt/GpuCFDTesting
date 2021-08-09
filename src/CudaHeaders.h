#ifndef CUDA_HEADERS_H
#define CUDA_HEADERS_H
#include <cuda.h>
#include "cuda_device_runtime_api.h"
#include "cuda_runtime_api.h"
#include <cuda_runtime.h>

#ifdef __CUDACC__
#define _f_hybrid __host__ __device__
#define _f_host __host__
#define _f_device __device__
#else
#define _f_hybrid
#define _f_host
#define _f_device
#endif

#define Cu_Check(mycode) mycode

#endif