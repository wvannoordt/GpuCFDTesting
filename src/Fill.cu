#include "Fill.h"
#include "print.h"
#include "MdArray.h"
#include "v3.h"
#include "Metric.h"
__global__ void K_FillTgv(MdArray<double, 4> flow, GasSpec gas, TgvSpec tgv, int nguard, Box box)
{
    int i = blockIdx.x*blockDim.x+threadIdx.x;
    int j = blockIdx.y*blockDim.y+threadIdx.y;
    int k = blockIdx.z*blockDim.z+threadIdx.z;
    int nvars=flow.dims[3];
    if (i<flow.dims[0]-nguard && j<flow.dims[1]-nguard && k<flow.dims[2]-nguard && i >= nguard && j >= nguard && k >= nguard)
    {
        v3<double> eta;
        eta[0] = (box.bounds[0]+((double)(i-nguard)+0.5)*box.dx[0]);
        eta[1] = (box.bounds[2]+((double)(j-nguard)+0.5)*box.dx[1]);
        eta[2] = (box.bounds[4]+((double)(k-nguard)+0.5)*box.dx[2]);
        v3<double> xyz;
        GetCoords(eta, xyz);
        for (int d = 0; d < 3; d++) xyz[d]/=tgv.L;
        double p = tgv.P0+((tgv.rho0*tgv.V0*tgv.V0)*(cos(2*xyz[0])+cos(2*xyz[1]))*(cos(2*xyz[2])+2.0))/16.0;
        double T = p/tgv.rho0*gas.R;
        double u = tgv.V0*sin(xyz[0])*cos(xyz[1])*cos(xyz[2]);
        double v = -tgv.V0*cos(xyz[0])*sin(xyz[1])*cos(xyz[2]);
        double w = 0.0;
        flow(i, j, k, nvars-1) = w;
        flow(i, j, k, 0) = p;
        flow(i, j, k, 1) = T;
        flow(i, j, k, 2) = u;
        flow(i, j, k, 3) = v;
    }
}

void FillTgv(FlowField& flow, const GasSpec& gas, const TgvSpec& tgv)
{
    dim3 grid = flow.GridConfig();
    dim3 block = flow.BlockConfig();
    for (int lb = 0; lb < flow.numBlocks; lb++)
    {
        auto array = flow.GetBlock(lb);
        auto box = flow.GetBox(lb);
        K_FillTgv<<<grid, block>>>(array, gas, tgv, flow.nguard, box);
    }
    cudaDeviceSynchronize();
    Cu_Check(cudaGetLastError());
}

__global__ void K_FillConst(MdArray<double, 4> flow, int nguard, double val)
{
    int i = blockIdx.x*blockDim.x+threadIdx.x;
    int j = blockIdx.y*blockDim.y+threadIdx.y;
    int k = blockIdx.z*blockDim.z+threadIdx.z;
    int nvars=flow.dims[3];
    if (i<flow.dims[0] && j<flow.dims[1] && k<flow.dims[2])
    {
        for (int y = 0; y < nvars; y++)
        {
            flow(i, j, k, y) = val;
        }
    }
}

void FillConst(FlowField& arr,  double val)
{
    dim3 grid = arr.GridConfig();
    dim3 block = arr.BlockConfig();
    for (int lb = 0; lb < arr.numBlocks; lb++)
    {
        auto array = arr.GetBlock(lb);
        K_FillConst<<<grid, block>>>(array, arr.nguard, val);
    }
    cudaDeviceSynchronize();
    Cu_Check(cudaGetLastError());
}