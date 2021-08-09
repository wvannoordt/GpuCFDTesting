#include "Fill.cuh"
#include "print.h"
#include "MdArray.h"
__global__ void K_FillTgv(MdArray<double, 4> flow, GasSpec gas, TgvSpec tgv, int nguard, int nvars)
{
    int i = blockIdx.x*blockDim.x+threadIdx.x+nguard;
    int j = blockIdx.y*blockDim.y+threadIdx.y+nguard;
    int k = blockIdx.z*blockDim.z+threadIdx.z+nguard;
    // if (i<flow.dims[0]-nguard && j<flow.dims[1]-nguard && k<flow.dims[2]-nguard && i >= nguard && j >= nguard && k >= nguard)
    if (i<flow.dims[0] && j<flow.dims[1] && k<flow.dims[2])
    {
        for (int y = 0; y < nvars; y++)
        {
            flow(i, j, k, y) = 10.4+y;
        }
    }
}

void FillTgv(FlowField& flow, const GasSpec& gas, const TgvSpec& tgv, const GpuConfig& config)
{
    dim3 grid = config.GridConfig();
    dim3 block = config.BlockConfig();
    for (int lb = 0; lb < flow.numBlocks; lb++)
    {
        auto array = flow.GetBlock(lb);
        K_FillTgv<<<grid, block>>>(array, gas, tgv, flow.nguard, flow.numVars);
    }
    cudaDeviceSynchronize();
}

// void FillConst(FlowField& arr,  double val)
// {
// 
// }