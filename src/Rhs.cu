#include "Rhs.h"

#define f_DivSplit(q,j,l,v1)          (0.500*(q((v1),(j)) +   q((v1),(j)+(l))))
#define fg_QuadSplit(q,j,l,v1,v2)     (0.250*(q((v1),(j)) +   q((v1),(j)+(l)))*(q((v2),(j)) + q((v2),(j)+(l))))
#define fg_CubeSplit(q,j,l,v1,v2,v3)  (0.125*(q((v1),(j)) +   q((v1),(j)+(l)))*(q((v2),(j)) + q((v2),(j)+(l))) * (q((v3),(j)) + q((v3),(j)+(l))))
#define fg_DivSplit(q,j,l,v1,v2)     (0.500*((q((v1),(j)+(l))*q((v2),(j))) +   (q((v1),(j)) * q((v2),(j)+(l)))))

__global__ void K_Rhs(MdArray<double, 4> rhsAr, MdArray<double, 4> flow, GasSpec gas, int nguard, Box box)
{
    int i = blockIdx.x*blockDim.x+threadIdx.x;
    int j = blockIdx.y*blockDim.y+threadIdx.y;
    int k = blockIdx.z*blockDim.z+threadIdx.z;
    int nvars=flow.dims[3];
    int dim = nvars-2;
    int dijk[3];
    double invdx[3] = {1.0, 1.0, 1.0};
    for (int pp = 0; pp < dim; pp++) invdx[pp] = box.dx[pp];
    const int centOrder = 4;
    const int stencilWid = centOrder/2;
    double coeff[9] = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
    //ie,ke,T,P,rho,u,v,w
    StaticArray<double, 5, 1+centOrder> stencilData;
    double Rgas = gas.cp*(gas.gamma - 1.0)/gas.gamma;
    if (i<flow.dims[0]-nguard && j<flow.dims[1]-nguard && k<flow.dims[2]-nguard && i >= nguard && j >= nguard && k >= nguard)
    {
        for (int idir = 0; idir < dim; idir++)
        {
            dijk[idir] = 1;
            
            double rhs[2+3] = {0.0};
            double C[2]     = {0.0};
            double M[3*2] = {0.0};
            double PGRAD[2] = {0.0};
            double IE[2]    = {0.0};
            double KE[2]    = {0.0};
            double PDIFF[2] = {0.0};
        
            for (int n = 0; n < centOrder + 1; n++)
            {
                for (int v = 3; v < (5+dim); v++)
                {
                    int ii = i+dijk[0]*(n-stencilWid);
                    int jj = j+dijk[1]*(n-stencilWid);
                    int kk = k+dijk[2]*(n-stencilWid);
                    stencilData(v,n) = flow(ii, jj, kk, v-3);
                }
                // T
                stencilData(2,n) = stencilData(3,n)/(Rgas*stencilData(4,n));
        
                // IE = P/(rho*(gamma - 1))
                stencilData(0,n) = stencilData(3,n)/(stencilData(4,n)*(gas.gamma - 1.0));
        
                // ke (don't care)
                stencilData(1,n) = 0.0;
        
                // Not needed per se starts
                for (int vel_comp = 0; vel_comp < dim; vel_comp ++)
                {
                    stencilData(1,n) += 0.5*stencilData(5+vel_comp,n)*stencilData(5+vel_comp,n);
                }
                // Not needed per se ends
            }
            // Mass conservation                              
            for (int l = 1; l <= stencilWid; l++)
            {
                double al = coeff[l-1];
                int jf = stencilWid;
                for (int m = 0; m <= (l-1); m++)
                {
                    C[1] += 2.0*al*fg_QuadSplit(stencilData,jf-m, l,4,5+idir);
                    C[0] += 2.0*al*fg_QuadSplit(stencilData,jf+m,-l,4,5+idir);
                    for (int idir_mom = 0; idir_mom < dim; idir_mom++)
                    {
                        M[idir_mom      ] += 2.0*al*fg_CubeSplit(stencilData,jf-m, l,4,5+idir,5+idir_mom);
                        M[idir_mom + dim] += 2.0*al*fg_CubeSplit(stencilData,jf+m,-l,4,5+idir,5+idir_mom);
                    }
        
                    PGRAD[1] += 2.0*al*f_DivSplit(stencilData,jf-m, l,3);
                    PGRAD[0] += 2.0*al*f_DivSplit(stencilData,jf+m,-l,3);
        
                    for (int vel_comp = 0;  vel_comp < dim; vel_comp ++)
                    {
                        KE[1] += 2.0*al*fg_QuadSplit(stencilData,jf-m, l,4,5+idir)*0.5*(stencilData(5+vel_comp,jf-m)*stencilData(5+vel_comp,jf-m+l));
                        KE[0] += 2.0*al*fg_QuadSplit(stencilData,jf+m,-l,4,5+idir)*0.5*(stencilData(5+vel_comp,jf+m)*stencilData(5+vel_comp,jf+m-l));
                    }
        
                    IE[1] += 2.0*al*fg_CubeSplit(stencilData,jf-m, l,4,0,5+idir);
                    IE[0] += 2.0*al*fg_CubeSplit(stencilData,jf+m,-l,4,0,5+idir);
        
                    PDIFF[1] += 2.0*al*fg_DivSplit(stencilData,jf-m, l,5+idir,3);
                    PDIFF[0] += 2.0*al*fg_DivSplit(stencilData,jf+m,-l,5+idir,3);
                }
            }
        
            rhs[0] += invdx[idir]*(C[1] - C[0]);
            rhs[1] += -invdx[idir]*(IE[1] + KE[1] + PDIFF[1] - IE[0] - KE[0] - PDIFF[0]);
            rhs[2+idir] += -invdx[idir]*(PGRAD[1] - PGRAD[0]);
            for (int rhs_vel_comp = 0; rhs_vel_comp < dim; rhs_vel_comp++)
            {
                rhs[2+rhs_vel_comp] -= invdx[idir]*(M[rhs_vel_comp] - M[rhs_vel_comp+dim]);
            }
            for (int nv = 0; nv < nvars; nv++) rhsAr(i, j, k, nv) +=  rhs[nv];
            dijk[idir] = 0;
        }
    }
}

void ComputeRhs(FlowField& rhs, FlowField& flow, const GasSpec& gas, const GpuConfig& config)
{
    dim3 grid = config.GridConfig();
    dim3 block = config.BlockConfig();
    for (int lb = 0; lb < flow.numBlocks; lb++)
    {
        auto flowArray = flow.GetBlock(lb);
        auto rhsArray = rhs.GetBlock(lb);
        auto box = flow.GetBox(lb);
        K_Rhs<<<grid, block>>>(rhsArray, flowArray, gas, flow.nguard, box);
    }
    cudaDeviceSynchronize();
    Cu_Check(cudaGetLastError());
}