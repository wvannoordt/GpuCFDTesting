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
    DataView<double, Dims<3, 2, 3, 3>> stencil;
    
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
            static_for<0,3>([&](auto idirLoop)
            {
                //idir0  = normal direction
                //idir1 = tangent direction
                //idir2 = other tangent direction
                int idir0 = idirLoop.value;
                int idir1 = (idir0+1)%3;
                int idir2 = (idir1+1)%3;
                for (int di0 = -plusMinus; di0 <= 1-plusMinus; di0++)
                {
                    ijk[idir0] += di0;
                    for (int di1 = -1; di1 <= 1; di1++)
                    {
                        ijk[idir1] += di1;
                        for (int di2 = -1; di2 <= 1; di2++)
                        {
                            ijk[idir2] += di2;
                            for (int vv = 0; vv < dim; vv++)
                            {
                                stencil(vv, di0+plusMinus, 1+di1, 1+di2) = flow(ijk[0], ijk[1], ijk[2], 2+vv);
                            }
                            ijk[idir2] -= di2;
                        }
                        ijk[idir1] -= di1;
                    }
                    ijk[idir0] -= di0;
                    
                    m9<double> faceVelGradComp;
                    
                    //du0/dxi0
                    for (int vv = 0; vv < dim; vv++)
                    {
                        faceVelGradComp(vv, idir0) = (stencil(vv, 1, 1, 1)-stencil(vv, 0, 1, 1))*invdx[idir0];
                        faceVelGradComp(vv, idir1) = 0.25*(stencil(vv, 1, 2, 1)-stencil(vv, 1, 0, 1) + stencil(vv, 0, 2, 1) - stencil(vv, 0, 0, 1))*invdx[idir1];
                        faceVelGradComp(vv, idir1) = 0.25*(stencil(vv, 1, 1, 2)-stencil(vv, 1, 1, 0) + stencil(vv, 0, 1, 2) - stencil(vv, 0, 1, 0))*invdx[idir2];
                    }
                    
                    m9<double> faceAverageMetrics;
                    m9<double> metricsNeighbor;
                    eta[idir0] += (1-2*plusMinus)*box.dx[idir0];
                    GetCoordsGrad(eta, metricsNeighbor);
                    eta[idir0] -= (1-2*plusMinus)*box.dx[idir0];
                    
                    for (int i2 = 0; i2<3; i2++)
                    {
                        for (int i1 = 0; i1<3; i1++)
                        {
                            faceAverageMetrics(i1,i2) = 0.5*(metricsNeighbor(i1,i2)+deta_dxyz(i1, i2));
                        }
                    }
                    
                    m9<double> faceVelGradPhys;
                    for (int i2 = 0; i2<3; i2++)
                    {
                        for (int i1 = 0; i1<3; i1++)
                        {
                            faceVelGradPhys(i1,i2) = faceVelGradComp(i1,i2)*faceAverageMetrics(i2,i1);
                        }
                    }
                    
                    m9<double> tau;
                    for (int i2 = 0; i2<3; i2++)
                    {
                        for (int i1 = 0; i1<3; i1++)
                        {
                            tau(i1, i2) = gas.visc*(faceVelGradPhys(i1,i2)+faceVelGradPhys(i2,i1));
                        }
                    }
                    
                    for (int i2 = 0; i2<3; i2++)
                    {
                        for (int i1 = 0; i1<3; i1++)
                        {
                            tau(i2, i2) += gas.beta*(faceVelGradPhys(i1,i1));
                        }
                    }
                    
                    for (int i1 = 0; i1<3; i1++)
                    {
                        rhsVals[2+idir0] -= (1-2*plusMinus)*faceAverageMetrics(idir0, i1)*tau(idir0, i1)*invdx[idir0];
                    }
                }
            });
        }
        
        for (int f = 0; f < 2+dim; f++)
        {
            rhsAr(i, j, k, f) -= jac*rhsVals[f];
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