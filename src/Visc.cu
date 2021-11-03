#include "Visc.h"
#include "Metric.h"
#include "StaticLoop.h"

__global__ void K_Visc(MdArray<double, 4> rhsAr, MdArray<double, 4> flow, GasSpec gas, int nguard, Box box)
{
    int i = blockIdx.x*blockDim.x+threadIdx.x;
    int j = blockIdx.y*blockDim.y+threadIdx.y;
    int k = blockIdx.z*blockDim.z+threadIdx.z;
    int nvars=flow.dims[3];
    int dim = nvars-2;
    double invdx[3] = {1.0, 1.0, 1.0};
    for (int pp = 0; pp < dim; pp++) invdx[pp] = 1.0/box.dx[pp];
    double Rgas = gas.R;
    v3<double> eta;
    m9<double> deta_dxyz;
    DataView<double, Dims<3, 3, 2>> stencil;
    
    v3<int> ijk(i, j, k);
    double rhsVals[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
    if (i<flow.dims[0]-nguard && j<flow.dims[1]-nguard && k<flow.dims[2]-nguard && i >= nguard && j >= nguard && k >= nguard)
    {
        eta[0] = (box.bounds[0]+((double)(i-nguard)+0.5)*box.dx[0]);
        eta[1] = (box.bounds[2]+((double)(j-nguard)+0.5)*box.dx[1]);
        eta[2] = (box.bounds[4]+((double)(k-nguard)+0.5)*box.dx[2]);
        GetCoordsGrad(eta, deta_dxyz);
        double jac = deta_dxyz.det();
        
        //plusMinus = 0 -> negative face
        //plusMinus = 1 -> positive face
        for (int plusMinus=0; plusMinus<=1; plusMinus++)
        {
            //Loop x, y, z
            static_for<0,3>([&](auto i)
            {
                //idir  = normal direction
                //idir1 = tangent direction
                //idir2 = other tangent direction
                int idir = i.value;
                // permutation<i.value, mod(i.value+1, 3), mod(i.value+2, 3)> perm;
            });
        }
    }
}

void ComputeViscRhs(FlowField& rhs, FlowField& flow, const GasSpec& gas)
{
    dim3 grid = rhs.GridConfig();
    dim3 block = rhs.BlockConfig();
    for (int lb = 0; lb < flow.numBlocks; lb++)
    {
        auto flowArray = flow.GetBlock(lb);
        auto rhsArray = rhs.GetBlock(lb);
        auto box = flow.GetBox(lb);
        K_Visc<<<grid, block>>>(rhsArray, flowArray, gas, flow.nguard, box);
    }
    cudaDeviceSynchronize(); 
    Cu_Check(cudaGetLastError());
}