#ifndef FLOW_FIELD_H
#define FLOW_FIELD_H
#include "CudaHeaders.h"
#include "CallError.h"
#include <vector>
#include "MdArray.h"
#include "Box.h"
#include "print.h"

struct FlowField
{
    double* data_d = NULL;
    double* dataBlock = NULL;
    
    int imin, imax, jmin, jmax, kmin, kmax;
    std::vector<std::string> varNames;
    int nguard;
    int numVars;
    std::vector<int> blockDim;
    std::vector<int> blockSize;
    size_t blockSizeBytes;
    int numBlocks;
    bool is3D;
        
    FlowField(const std::vector<int>& blockDim_in, const std::vector<int>& blockSize_in, const int numVars_in, const int nguard_in)
    {
        nguard = nguard_in;
        numVars = numVars_in;
        blockDim = blockDim_in;
        blockSize = blockSize_in;
        
        if ((blockDim.size() != 2) && (blockDim.size() != 3)) CallError("Invalid vector size for blockDim and blockSize");
        imin = 0;
        jmin = 0;
        kmin = 0;
        imax = blockSize[0]+2*nguard;
        jmax = blockSize[1]+2*nguard;
        kmax = 1;
        if (blockSize.size()==3)
        {
            kmax = blockSize[2]+2*nguard;
        }
        is3D = kmax!=1;
        
        // (v, i, j, k, lb) vs (i, j, k, v, lb)
        this->blockSizeBytes = sizeof(double);
        for (auto p: blockSize)
        {
            blockSizeBytes *= p+2*nguard;
        }
        blockSizeBytes *= numVars;
        numBlocks = 1;
        for (auto p: blockDim) numBlocks *= p;
        size_t totalAllocSize = blockSizeBytes * numBlocks;
        cudaMallocHost((void**)(&dataBlock), blockSizeBytes);
        cudaMalloc((void**)(&data_d), totalAllocSize);
        for (int i = 0; i < numVars; i++)
        {
            std::string name = "var"+std::to_string(i);
            varNames.push_back(name);
        }
    }
    
    
    MdArray<double, 4> GetBlock(int lb) const
    {
        MdArray<double, 4> output(imax, jmax, kmax, numVars);
        output.data = data_d + lb*(blockSizeBytes/sizeof(double));
        return output;
    }
    
    MdArray<double, 4> UnloadBlock(int lb) const
    {
        MdArray<double, 4> output(imax, jmax, kmax, numVars);
        Cu_Check(cudaMemcpy((void*)dataBlock, (void*)(data_d + lb*(blockSizeBytes/sizeof(double))), this->blockSizeBytes, cudaMemcpyDeviceToHost));
        output.data = dataBlock;
        return output;
    }
    
    Box GetBox(int lb) const
    {
        int ijkbox[3] = {0};
        ijkbox[0] = lb%blockDim[0];
        ijkbox[1] = ((lb-ijkbox[0])/blockDim[0])%blockDim[1];
        ijkbox[2] = 0;
        if (is3D)
        {
            ijkbox[2] = (lb-ijkbox[0]-blockDim[0]*ijkbox[1])/(blockDim[0]*blockDim[1]);
        }
        double dxdomain[3] = {1.0, 1.0, 1.0};
        dxdomain[2] = 1.0;
        for (int i = 0; i < (is3D?3:2); i++)
        {
            dxdomain[i] = 1.0/blockDim[i];
        }
        Box output;
        output.bounds[4] = 0.0;
        output.bounds[5] = 1.0;
        for (int i = 0; i < (is3D?3:2); i++)
        {
            output.bounds[2*i] = ijkbox[i]*dxdomain[i];
            output.bounds[2*i+1] = (ijkbox[i]+1)*dxdomain[i];
        }
        output.dx[2] = 1.0;
        for (int i = 0; i < (is3D?3:2); i++)
        {
            output.dx[i] = (output.bounds[2*i+1]-output.bounds[2*i])/blockSize[i];
        }
        return output;
    }
    
    dim3 BlockConfig(void) const
    {
        dim3 output;
        output.x = CU_BLOCK_SIZE;
        output.y = CU_BLOCK_SIZE;
        output.z = 1;
        if (this->is3D)
        {
            output.z = CU_BLOCK_SIZE;
        }
        return output;
    }
    
    dim3 GridConfig(void) const
    {
        dim3 output;
        dim3 gridConf = this->BlockConfig();
        output.x = (this->blockSize[0]+2*this->nguard + gridConf.x - 1)/gridConf.x;
        output.y = (this->blockSize[1]+2*this->nguard + gridConf.y - 1)/gridConf.y;
        output.z = 1;
        if (this->is3D)
        {
            output.z = (this->blockSize[2]+2*this->nguard + gridConf.z - 1)/gridConf.z;
        }
        return output;
    }
    
    ~FlowField(void)
    {
        cudaFreeHost(dataBlock);
        cudaFree(data_d);
    }
};

#endif