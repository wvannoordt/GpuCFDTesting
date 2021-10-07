#ifndef INPUT_CLASS_H
#define INPUT_CLASS_H
#include "PTL.h"
#include <vector>
#include "CallError.h"
struct InputClass
{
    std::vector<int> blockDim;
    std::vector<int> blockSize;
    int nguard;
    int numVars;
    bool is3D;
    int dim;
    void Read(PTL::PropertySection& section)
    {
        section["blockDim"].MapTo(&blockDim)         = new PTL::PTLIntegerVector("Dimensions of the block configuration");
        section["blockSize"].MapTo(&blockSize)       = new PTL::PTLIntegerVector("Dimensions of each individual block");
        section["nguard"].MapTo(&nguard)             = new PTL::PTLInteger(2, "Number of exchange cells");
        section.StrictParse();
        
        if (blockDim.size() != blockSize.size()) CallError("blockDim and blockSize have different dimensions!");
        if ((blockDim.size() != 2) && (blockDim.size() != 3)) CallError("Invalid vector size for blockDim and blockSize");
        
        dim = 2;
        is3D = false;
        numVars = 4;
        
        if (blockDim.size()==3)
        {
            dim = 3;
            is3D = true;
            numVars = 5;
        }
        
    }
};
#endif