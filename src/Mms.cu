#include "Mms.h"
#include "v3.h"
#include "Metric.h"
#include "Conv.h"
#include "Fill.h"
#include "Output.h"
__global__ void K_MmsConvRhs(NavierStokesMms mms, MdArray<double, 4> rhs, int nguard, Box box)
{
    int i = blockIdx.x*blockDim.x+threadIdx.x;
    int j = blockIdx.y*blockDim.y+threadIdx.y;
    int k = blockIdx.z*blockDim.z+threadIdx.z;
    int nvars=rhs.dims[3];
    if (i<rhs.dims[0]-nguard && j<rhs.dims[1]-nguard && k<rhs.dims[2]-nguard && i >= nguard && j >= nguard && k >= nguard)
    {
        v3<double> eta;
        eta[0] = (box.bounds[0]+((double)(i-nguard)+0.5)*box.dx[0]);
        eta[1] = (box.bounds[2]+((double)(j-nguard)+0.5)*box.dx[1]);
        eta[2] = (box.bounds[4]+((double)(k-nguard)+0.5)*box.dx[2]);
        v3<double> xyz;
        GetCoords(eta, xyz);
        double rhsAna[5];
        mms.conv_rhs(xyz[0], xyz[1], xyz[2], rhsAna);
        for (int nv = 0; nv < nvars; nv++) rhs(i, j, k, nv) -= rhsAna[nv];
    }
}

__global__ void K_MmsTestFcn(NavierStokesMms mms, MdArray<double, 4> prims, int nguard, Box box)
{
    int i = blockIdx.x*blockDim.x+threadIdx.x;
    int j = blockIdx.y*blockDim.y+threadIdx.y;
    int k = blockIdx.z*blockDim.z+threadIdx.z;
    int nvars=prims.dims[3];
    if (i<prims.dims[0] && j<prims.dims[1] && k<prims.dims[2])
    {
        v3<double> eta;
        eta[0] = (box.bounds[0]+((double)(i-nguard)+0.5)*box.dx[0]);
        eta[1] = (box.bounds[2]+((double)(j-nguard)+0.5)*box.dx[1]);
        eta[2] = (box.bounds[4]+((double)(k-nguard)+0.5)*box.dx[2]);
        v3<double> xyz;
        GetCoords(eta, xyz);
        double primsAna[5];
        mms.testFcn(xyz[0], xyz[1], xyz[2], primsAna);
        for (int nv = 0; nv < nvars; nv++) prims(i, j, k, nv) = primsAna[nv];
    }
}

void AnalyticalConvRhs(const NavierStokesMms& mms, FlowField& rhs, const GpuConfig& config)
{
    dim3 grid = config.GridConfig();
    dim3 block = config.BlockConfig();
    for (int lb = 0; lb < rhs.numBlocks; lb++)
    {
        auto array = rhs.GetBlock(lb);
        auto box = rhs.GetBox(lb);
        K_MmsConvRhs<<<grid, block>>>(mms, array, rhs.nguard, box);
    }
    cudaDeviceSynchronize();
    Cu_Check(cudaGetLastError());
}

void AnalyticalFcn(const NavierStokesMms& mms, FlowField& prims, const GpuConfig& config)
{
    dim3 grid = config.GridConfig();
    dim3 block = config.BlockConfig();
    for (int lb = 0; lb < prims.numBlocks; lb++)
    {
        auto array = prims.GetBlock(lb);
        auto box = prims.GetBox(lb);
        K_MmsTestFcn<<<grid, block>>>(mms, array, prims.nguard, box);
    }
    cudaDeviceSynchronize();
    Cu_Check(cudaGetLastError());
}

void RunMMS(FlowField& prims, FlowField& rhs, const GasSpec& gas, const GpuConfig& config)
{
    print("MMS");
    NavierStokesMms mms(gas);
    print("Fill analytical function");
    AnalyticalFcn(mms, prims, config);
    
    print("Zero RHS");
    FillConst(rhs, 0.0, config);
    
    print("Compute RHS");
    ComputeConvRhs(rhs, prims, gas, config);
    
    print("Analytical RHS");
    AnalyticalConvRhs(mms, rhs, config);
    Output(rhs, "output", "RHS_conv_error");
    CallError("Finished MMS");
}