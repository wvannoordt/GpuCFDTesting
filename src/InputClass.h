#ifndef INPUT_CLASS_H
#define INPUT_CLASS_H
#include "PTL.h"
#include <vector>
#include "CallError.h"
struct InputClass
{
    std::vector<int> blockDim;
    std::vector<int> blockSize;
    
    void Read(PTL::PropertySection& section)
    {
        section["blockDim"].MapTo(&blockDim)   = new PTL::PTLIntegerVector("Dimensions of the block configuration");
        section["blockSize"].MapTo(&blockSize) = new PTL::PTLIntegerVector("Dimensions of each individual block");
        section.StrictParse();
        
        if (blockDim.size() != blockSize.size()) CallError("blockDim and blockSize have different dimensions!");
    }
};
#endif