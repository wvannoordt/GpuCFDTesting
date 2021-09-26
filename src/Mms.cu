#include "Mms.h"

__global__ void K_MmsRhs(NavierStokesMms mms, MdArray<double, 4> rhs, int nguard, Box box)
{
    int i = blockIdx.x*blockDim.x+threadIdx.x;
    int j = blockIdx.y*blockDim.y+threadIdx.y;
    int k = blockIdx.z*blockDim.z+threadIdx.z;
    int nvars=rhs.dims[3];
    if (i<rhs.dims[0]-nguard && j<rhs.dims[1]-nguard && k<rhs.dims[2]-nguard && i >= nguard && j >= nguard && k >= nguard)
    {
        double x = (box.bounds[0]+((double)(i-nguard)+0.5)*box.dx[0]);
        double y = (box.bounds[2]+((double)(j-nguard)+0.5)*box.dx[1]);
        double z = (box.bounds[4]+((double)(k-nguard)+0.5)*box.dx[2]);
        double rhsAna[5];
        mms.rhs(x, y, z, rhsAna);
        for (int nv = 0; nv < nvars; nv++) rhs(i, j, k, nv) = rhsAna[nv];
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
        double x = (box.bounds[0]+((double)(i-nguard)+0.5)*box.dx[0]);
        double y = (box.bounds[2]+((double)(j-nguard)+0.5)*box.dx[1]);
        double z = (box.bounds[4]+((double)(k-nguard)+0.5)*box.dx[2]);
        double primsAna[5];
        mms.testFcn(x, y, z, primsAna);
        for (int nv = 0; nv < nvars; nv++) prims(i, j, k, nv) = primsAna[nv];
    }
}

void AnalyticalRhs(const NavierStokesMms& mms, FlowField& rhs, const GpuConfig& config)
{
    dim3 grid = config.GridConfig();
    dim3 block = config.BlockConfig();
    for (int lb = 0; lb < rhs.numBlocks; lb++)
    {
        auto array = rhs.GetBlock(lb);
        auto box = rhs.GetBox(lb);
        K_MmsRhs<<<grid, block>>>(mms, array, rhs.nguard, box);
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
