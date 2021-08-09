#ifndef FLOW_FIELD_H
#define FLOW_FIELD_H
#include "CudaHeaders.h"
#include "CallError.h"
#include <vector>
#include "MdArray.h"
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
        dataBlock = (double*)malloc(blockSizeBytes);
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
    
    ~FlowField(void)
    {
        free(dataBlock);
        cudaFree(data_d);
    }
};

#endif