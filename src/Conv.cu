#include "Conv.h"
#include "Metric.h"
#include "StaticLoop.h"

#define f_DivSplit(q,j,l,v1)           (0.500*(q((v1),(j)) +   q((v1),(j)+(l))))
#define f_DivSplit2(q,j,l,v1,v2)          (0.500*(q((v1),(j))*q((v2),(j)) +   q((v1),(j)+(l))*q((v2),(j)+(l))))
#define fg_QuadSplit(q,j,l,v1,v2)      (0.250*(q((v1),(j)) +   q((v1),(j)+(l)))*(q((v2),(j)) + q((v2),(j)+(l))))
#define fg_CubeSplit(q,j,l,v1,v2,v3)   (0.125*(q((v1),(j)) +   q((v1),(j)+(l)))*(q((v2),(j)) + q((v2),(j)+(l))) * (q((v3),(j)) + q((v3),(j)+(l))))
#define fg_DivSplit(q,j,l,v1,v2)       (0.500*((q((v1),(j)+(l))*q((v2),(j))) +   (q((v1),(j)) * q((v2),(j)+(l)))))

__global__ void K_Conv_Central(MdArray<double, 4> rhsAr, MdArray<double, 4> flow, GasSpec gas, int nguard, Box box)
{
    int i = blockIdx.x*blockDim.x+threadIdx.x;
    int j = blockIdx.y*blockDim.y+threadIdx.y;
    int k = blockIdx.z*blockDim.z+threadIdx.z;
    int nvars=flow.dims[3];
    int dim = nvars-2;
    int dijk[3];
    double invdx[3] = {1.0, 1.0, 1.0};
    for (int pp = 0; pp < dim; pp++) invdx[pp] = 1.0/box.dx[pp];
    const int centOrder = 4;
    const int stencilWid = centOrder/2;
    const double coeff[9] = {2.0/3.0, -1.0/12.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    //ie,ke,T,P,rho,u,v,w,eta_x, eta_y, eta_z
    DataView<double, Dims<11,1+centOrder>> stencilData;
    double Rgas = gas.R;
    v3<double> eta;
    m9<double> deta_dxyz;
    if (i<flow.dims[0]-nguard && j<flow.dims[1]-nguard && k<flow.dims[2]-nguard && i >= nguard && j >= nguard && k >= nguard)
    {
        for (int idir = 0; idir < dim; idir++)
        {
            dijk[idir] = 1;
            double rhsArr[5] = {0.0};
        
            static_for<0,centOrder + 1>([&](auto nloop)
            {
                int n = nloop.value;
                int ii = i+dijk[0]*(n-stencilWid);
                int jj = j+dijk[1]*(n-stencilWid);
                int kk = k+dijk[2]*(n-stencilWid);
                eta[0] = (box.bounds[0]+((double)(ii-nguard)+0.5)*box.dx[0]);
                eta[1] = (box.bounds[2]+((double)(jj-nguard)+0.5)*box.dx[1]);
                eta[2] = (box.bounds[4]+((double)(kk-nguard)+0.5)*box.dx[2]);
                GetCoordsGrad(eta, deta_dxyz);
                double stencilJac = deta_dxyz.det();
                
                stencilData(8, n)  = deta_dxyz(idir, 0)/stencilJac;
                stencilData(9, n)  = deta_dxyz(idir, 1)/stencilJac;
                stencilData(10, n) = deta_dxyz(idir, 2)/stencilJac;
                //ie,U_Contra,T,P,rho,u,v,w
                //0  1        2 3 4   5 6 7
                for (int v = 3; v < (5+dim); v++)
                {
                    stencilData(v,n) = flow(ii, jj, kk, v-3);
                }
                // T
                stencilData(2,n) = stencilData(4, n);
                
                //rho
                stencilData(4, n) = stencilData(3,n)/(Rgas*stencilData(2,n));
        
                // IE = P/(rho*(gamma - 1))
                stencilData(0,n) = stencilData(3,n)/(stencilData(4,n)*(gas.gamma - 1.0));
        
                // U_Contra
                stencilData(1,n) = 0.0;
                for (int vel_comp = 0; vel_comp < dim; vel_comp ++)
                {
                    stencilData(1,n) += stencilData(5+vel_comp,n)*deta_dxyz(idir, vel_comp)/stencilJac;
                }
            });
            eta[0] = (box.bounds[0]+((double)(i-nguard)+0.5)*box.dx[0]);
            eta[1] = (box.bounds[2]+((double)(j-nguard)+0.5)*box.dx[1]);
            eta[2] = (box.bounds[4]+((double)(k-nguard)+0.5)*box.dx[2]);
            GetCoordsGrad(eta, deta_dxyz);
            double jac = deta_dxyz.det();
            // Mass conservation
            for (int l = 1; l <= stencilWid; l++)
            {
                double al = coeff[l-1];
                int jf = stencilWid;
                for (int m = 0; m <= (l-1); m++)
                {
                    rhsArr[0] += 2.0*al*fg_QuadSplit(stencilData,jf-m, l,4,1);
                    rhsArr[0] -= 2.0*al*fg_QuadSplit(stencilData,jf+m,-l,4,1);
                    for (int idir_mom = 0; idir_mom < dim; idir_mom++)
                    {
                        rhsArr[2+idir_mom] += 2.0*al*fg_CubeSplit(stencilData,jf-m, l,4,1,5+idir_mom);
                        rhsArr[2+idir_mom] -= 2.0*al*fg_CubeSplit(stencilData,jf+m,-l,4,1,5+idir_mom);
                        
                        rhsArr[2+idir_mom] += 2.0*al*f_DivSplit2(stencilData,jf-m, l,3,8+idir_mom);
                        rhsArr[2+idir_mom] -= 2.0*al*f_DivSplit2(stencilData,jf+m,-l,3,8+idir_mom);
                    }
        
                    for (int vel_comp = 0;  vel_comp < dim; vel_comp ++)
                    {
                        rhsArr[1] += 2.0*al*fg_QuadSplit(stencilData,jf-m, l,4,1)*0.5*(stencilData(5+vel_comp,jf-m)*stencilData(5+vel_comp,jf-m+l));
                        rhsArr[1] -= 2.0*al*fg_QuadSplit(stencilData,jf+m,-l,4,1)*0.5*(stencilData(5+vel_comp,jf+m)*stencilData(5+vel_comp,jf+m-l));
                    }
        
                    rhsArr[1] += 2.0*al*fg_CubeSplit(stencilData,jf-m, l,4,0,1);
                    rhsArr[1] -= 2.0*al*fg_CubeSplit(stencilData,jf+m,-l,4,0,1);
        
                    rhsArr[1] += 2.0*al*fg_DivSplit(stencilData,jf-m, l,1,3);
                    rhsArr[1] -= 2.0*al*fg_DivSplit(stencilData,jf+m,-l,1,3);
                }
            }
            for (int p = 0; p < nvars; p++) rhsAr(i, j, k, p) -= jac*invdx[idir]*rhsArr[p];
            
            dijk[idir] = 0;
        }
    }
}

void ComputeConvRhs(FlowField& rhs, FlowField& flow, const GasSpec& gas)
{
    dim3 grid = rhs.GridConfig();
    dim3 block = rhs.BlockConfig();
    for (int lb = 0; lb < flow.numBlocks; lb++)
    {
        auto flowArray = flow.GetBlock(lb);
        auto rhsArray = rhs.GetBlock(lb);
        auto box = flow.GetBox(lb);
        K_Conv_Central<<<grid, block>>>(rhsArray, flowArray, gas, flow.nguard, box);
    }
    cudaDeviceSynchronize(); 
    Cu_Check(cudaGetLastError());
}