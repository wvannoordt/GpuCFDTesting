#ifndef GPU_CONFIG_H
#define GPU_CONFIG_H
#include "PTL.h"
#include <vector>
#include "CudaHeaders.h"
#include "InputClass.h"
struct GpuConfig
{
    std::vector<int> threadDimGpu;
    InputClass* input = NULL;
    GpuConfig(InputClass& input_in)
    {
        input = &input_in;
    }
    void Read(PTL::PropertySection& section)
    {
        section["threadDimGpu"].MapTo(&threadDimGpu) = new PTL::PTLIntegerVector("Dimensions of each individual GPU thread block");
        section.StrictParse();

        if ((threadDimGpu.size() != 2) && (threadDimGpu.size() != 3)) CallError("Invalid vector size for blockDimGpu");
        if (threadDimGpu.size() != input->dim) CallError("threadDimGpu and gridDimGpu do not agree with the provided dimensions for mesh");
    }
    
    dim3 BlockConfig(void) const
    {
        dim3 output;
        output.x = threadDimGpu[0];
        output.y = threadDimGpu[1];
        output.z = 1;
        if (input->dim==3)
        {
            output.z = threadDimGpu[2];
        }
        return output;
    }
    
    dim3 GridConfig(void) const
    {
        dim3 output;
        dim3 gridConf = this->BlockConfig();
        output.x = (input->blockSize[0]+2*input->nguard + gridConf.x - 1)/gridConf.x;
        output.y = (input->blockSize[1]+2*input->nguard + gridConf.y - 1)/gridConf.y;
        output.z = 1;
        if (input->dim==3)
        {
            output.z = (input->blockSize[2]+2*input->nguard + gridConf.z - 1)/gridConf.z;
        }
        return output;
    }
};

#endif