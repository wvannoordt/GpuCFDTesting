#include "Fill.cuh"
#include "print.h"
#include "MdArray.h"
__global__ void K_FillTgv(MdArray<double, 4> flow, GasSpec gas, TgvSpec tgv, int nguard, int nvars, Box box)
{
    int i = blockIdx.x*blockDim.x+threadIdx.x;
    int j = blockIdx.y*blockDim.y+threadIdx.y;
    int k = blockIdx.z*blockDim.z+threadIdx.z;
    if (i<flow.dims[0]-nguard && j<flow.dims[1]-nguard && k<flow.dims[2]-nguard && i >= nguard && j >= nguard && k >= nguard)
    {
        double x = (box.bounds[0]+i*box.dx[0])/tgv.L;
        double y = (box.bounds[2]+j*box.dx[1])/tgv.L;
        double z = (box.bounds[4]+k*box.dx[2])/tgv.L;
        double p = tgv.P0+((tgv.rho0*tgv.V0*tgv.V0)*(cos(2*x)+cos(2*y))*(cos(2*z)+2.0))/16.0;
        double T = p/tgv.rho0*gas.R;
        double u = tgv.V0*sin(x)*cos(y)*cos(z);
        double v = -tgv.V0*cos(x)*sin(y)*cos(z);
        double w = 0.0;
        flow(i, j, k, nvars-1) = w;
        flow(i, j, k, 0) = p;
        flow(i, j, k, 1) = T;
        flow(i, j, k, 2) = u;
        flow(i, j, k, 3) = v;
    }
}

void FillTgv(FlowField& flow, const GasSpec& gas, const TgvSpec& tgv, const GpuConfig& config)
{
    dim3 grid = config.GridConfig();
    dim3 block = config.BlockConfig();
    for (int lb = 0; lb < flow.numBlocks; lb++)
    {
        auto array = flow.GetBlock(lb);
        auto box = flow.GetBox(lb);
        K_FillTgv<<<grid, block>>>(array, gas, tgv, flow.nguard, flow.numVars, box);
    }
    cudaDeviceSynchronize();
}

// void FillConst(FlowField& arr,  double val)
// {
// 
// }