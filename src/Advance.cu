#include "Advance.h"

__device__ void prim2cons(double (&prim)[5], double (&cons)[5], const GasSpec& gas)
{
    double p = prim[0];
    double T = prim[1];
    double u = prim[2];
    double v = prim[3];
    double w = prim[4];
    double Rgas = gas.cp*(gas.gamma - 1.0)/gas.gamma;
    double rho = p / (Rgas*T);
    double rhoU2 = rho*(u*u+v*v+w*w);
    double rhoE = 0.5*rhoU2 + (p/((gas.gamma - 1.0)));
    double rhoU = rho*u;
    double rhoV = rho*v;
    double rhoW = rho*w;
    cons[0] = rho;
    cons[1] = rhoE;
    cons[2] = rhoU;
    cons[3] = rhoV;
    cons[4] = rhoW;
}

__device__ void cons2prim(double (&cons)[5], double (&prim)[5], const GasSpec& gas)
{
    double rho = cons[0];
    double invrho = 1.0/rho;
    double u = invrho*cons[2];
    double v = invrho*cons[3];
    double w = invrho*cons[4];
    double rhoU2 = rho*(u*u+v*v+w*w);
    double p = (gas.gamma - 1.0)*(cons[1] - 0.5*rhoU2);
    double Rgas = gas.cp*(gas.gamma - 1.0)/gas.gamma;
    double T = p/(Rgas*rho);
    prim[0] = p;
    prim[1] = T;
    prim[2] = u;
    prim[3] = v;
    prim[4] = w;
}

__global__ void K_Advance(MdArray<double, 4> rhs, MdArray<double, 4> flow, GasSpec gas, double timestep, int nguard, Box box)
{
    int i = blockIdx.x*blockDim.x+threadIdx.x;
    int j = blockIdx.y*blockDim.y+threadIdx.y;
    int k = blockIdx.z*blockDim.z+threadIdx.z;
    int nvars=flow.dims[3];
    double prims[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
    double cons[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
    if (i<flow.dims[0]-nguard && j<flow.dims[1]-nguard && k<flow.dims[2]-nguard && i >= nguard && j >= nguard && k >= nguard)
    {
        for (int nv = 0; nv < nvars; nv++) prims[nv] = flow(i, j, k, nv);
        prim2cons(prims, cons, gas);
        for (int nv = 0; nv < nvars; nv++) cons[nv] += timestep*rhs(i, j, k, nv);
        cons2prim(cons, prims, gas);
        for (int nv = 0; nv < nvars; nv++) flow(i, j, k, nv) = prims[nv];
    }
}

void Advance(FlowField& rhs, FlowField& flow, double timestep, const GasSpec& gas)
{
    dim3 grid = rhs.GridConfig();
    dim3 block = rhs.BlockConfig();
    for (int lb = 0; lb < flow.numBlocks; lb++)
    {
        auto flowArray = flow.GetBlock(lb);
        auto rhsArray = rhs.GetBlock(lb);
        auto box = flow.GetBox(lb);
        K_Advance<<<grid, block>>>(rhsArray, flowArray, gas, timestep, flow.nguard, box);
    }
    cudaDeviceSynchronize();
    Cu_Check(cudaGetLastError());
}