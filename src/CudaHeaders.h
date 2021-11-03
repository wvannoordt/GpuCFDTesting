#ifndef CUDA_HEADERS_H
#define CUDA_HEADERS_H
#include <cuda.h>
#include "cuda_device_runtime_api.h"
#include "cuda_runtime_api.h"
#include <cuda_runtime.h>
#include "CallError.h"
#include <string>
#ifdef __CUDACC__
#define _f_hybrid __host__ __device__
#define _f_host __host__
#define _f_device __device__
#define _f_no_inline __noinline__
#else
#define _f_hybrid
#define _f_host
#define _f_device
#define _f_no_inline
#endif

#define Cu_Check(mycode) {cudaError_t ____code = mycode; if (____code!=cudaSuccess) {CallError(std::string("CUDA Runtime Error. Message:\n") + cudaGetErrorString(____code));}}
#define CU_BLOCK_SIZE 4

#endif