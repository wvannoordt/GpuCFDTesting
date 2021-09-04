#include "Exchange.h"
#include "print.h"
#include "CallError.h"
#include "v3.h"
int GetNeighbor(FlowField& arr, int lb, int di, int dj, int dk)
{
    int ni = arr.blockDim[0];
    int nj = arr.blockDim[1];
    int nk = 1;
    if (arr.is3D) nk = arr.blockDim.back();
    int i = lb%ni;
    int j = ((lb-i)/ni)%nj;
    int k = ((lb-i-j*ni)/(ni*nj))%nk;
    i = (i + di + ni)%ni;
    j = (j + dj + nj)%nj;
    k = (k + dk + nk)%nk;
    return i+j*ni+k*ni*nj;
}

__global__ void K_ExchangeRegion(MdArray<double, 4> send, MdArray<double, 4> recv, v3<int> ijkminSend, v3<int> ijkminRecv, v3<int> exchangeRegionSize)
{
    int i = blockIdx.x*blockDim.x+threadIdx.x;
    int j = blockIdx.y*blockDim.y+threadIdx.y;
    int k = blockIdx.z*blockDim.z+threadIdx.z;
    int nvars=send.dims[3];
    if (i<exchangeRegionSize[0]&&j<exchangeRegionSize[1]&&k<exchangeRegionSize[2])
    {
        for (int v = 0; v < nvars; v++)
        {
            recv(i+ijkminRecv[0], j+ijkminRecv[1], k+ijkminRecv[2], v) = send(i+ijkminSend[0], j+ijkminSend[1], k+ijkminSend[2], v);
        }
    }
}

void ExchangeSingleBlock(FlowField& arr, int send, int recv, int dijk[3], GpuConfig& config)
{
    //send + dijk = recv
    auto sendBlock = arr.GetBlock(send);
    auto recvBlock = arr.GetBlock(recv);
    v3<int> ijkminRecv(0, 0, 0);
    v3<int> ijkminSend(0, 0, 0);
    v3<int> exchangeRegionSize(1, 1, 1);
    int ndim = 2;
    if (arr.is3D) ndim = 3;
    for (int i = 0; i < ndim; i++)
    {
        bool dijk0 = (dijk[i]==0);
        switch (dijk[i])
        {
            case -1:
            {
                exchangeRegionSize[i] = arr.nguard;
                ijkminSend[i] = arr.nguard;
                ijkminRecv[i] = arr.blockSize[i]+arr.nguard;
                break;
            }
            case 0:
            {
                exchangeRegionSize[i] = arr.blockSize[i];
                ijkminSend[i] = arr.nguard;
                ijkminRecv[i] = arr.nguard;
                break;
            }
            case 1:
            {
                exchangeRegionSize[i] = arr.nguard;
                ijkminSend[i] = arr.blockSize[i];
                ijkminRecv[i] = 0;
                break;
            }
        }            
    }
    
    dim3 block = config.BlockConfig();
    dim3 grid(1,1,1);
    grid.x = (exchangeRegionSize[0] + block.x - 1)/block.x;
    grid.y = (exchangeRegionSize[1] + block.y - 1)/block.y;
    grid.z = 1;
    if (arr.is3D)
    {
        grid.z = (exchangeRegionSize[2] + block.z - 1)/block.z;
    }
    
    K_ExchangeRegion<<<grid, block>>>(sendBlock, recvBlock, ijkminSend, ijkminRecv, exchangeRegionSize);
}

void Exchange(FlowField& arr, GpuConfig& config)
{
    int dijk[3] = {0, 0, 0};
    int numNeighbors = 9;
    if (arr.is3D) numNeighbors = 27;
    for (int lb = 0; lb < arr.numBlocks; lb++)
    {
        for (int neighId = 0; neighId < numNeighbors; neighId++)
        {
            dijk[0] = neighId%3;
            dijk[1] = ((neighId-dijk[0])/3)%3;
            dijk[2] = ((neighId-dijk[0]-3*dijk[1])/9)%3;            
            dijk[0] -= 1;
            dijk[1] -= 1;
            if (arr.is3D) dijk[2] -= 1;
            if (dijk[0]!=0||dijk[1]!=0||dijk[2]!=0)
            {
                int lbNeigh = GetNeighbor(arr, lb, dijk[0], dijk[1], dijk[2]);
                ExchangeSingleBlock(arr, lb, lbNeigh, dijk, config);
            }
        }
    }
    cudaDeviceSynchronize();
    Cu_Check(cudaGetLastError());
}