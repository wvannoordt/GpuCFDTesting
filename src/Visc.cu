#include "Visc.h"
#include "Metric.h"

#define f_DivSplit(q,j,l,v1)           (0.500*(q((v1),(j)) +   q((v1),(j)+(l))))
#define f_DivSplit2(q,j,l,v1,v2)          (0.500*(q((v1),(j))*q((v2),(j)) +   q((v1),(j)+(l))*q((v2),(j)+(l))))
#define fg_QuadSplit(q,j,l,v1,v2)      (0.250*(q((v1),(j)) +   q((v1),(j)+(l)))*(q((v2),(j)) + q((v2),(j)+(l))))
#define fg_CubeSplit(q,j,l,v1,v2,v3)   (0.125*(q((v1),(j)) +   q((v1),(j)+(l)))*(q((v2),(j)) + q((v2),(j)+(l))) * (q((v3),(j)) + q((v3),(j)+(l))))
#define fg_DivSplit(q,j,l,v1,v2)       (0.500*((q((v1),(j)+(l))*q((v2),(j))) +   (q((v1),(j)) * q((v2),(j)+(l)))))

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
    DataView<double, Dims<3, 3, 3>> stencil;
    if (i<flow.dims[0]-nguard && j<flow.dims[1]-nguard && k<flow.dims[2]-nguard && i >= nguard && j >= nguard && k >= nguard)
    {
        eta[0] = (box.bounds[0]+((double)(i-nguard)+0.5)*box.dx[0]);
        eta[1] = (box.bounds[2]+((double)(j-nguard)+0.5)*box.dx[1]);
        eta[2] = (box.bounds[4]+((double)(k-nguard)+0.5)*box.dx[2]);
        GetCoordsGrad(eta, deta_dxyz);
        double jac = deta_dxyz.det();
        
        //Loop x, y, z
        for (int idir = 0; idir < dim; idir++)
        {
        
        }
    }
}

void ComputeViscRhs(FlowField& rhs, FlowField& flow, const GasSpec& gas, const GpuConfig& config)
{
    dim3 grid = config.GridConfig();
    dim3 block = config.BlockConfig();
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