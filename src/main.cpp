#include <iostream>
#include <string>
#include "InputClass.h"
#include "print.h"
#include "FlowField.h"

std::string GetInputFilename(int argc, char** argv)
{
    if (argc<2) return "input.ptl";
    return std::string(argv[1]);
}

int main(int argc, char** argv)
{
    std::string filename = GetInputFilename(argc, argv);
    PTL::PropertyTree tree;
    tree.Read(filename);
    
    InputClass input;
    input.Read(tree["Grid"]);
    
    FlowField flow(input.blockDim, input.blockSize, input.numVars, input.nguard);
    flow.varNames[0] = "P";
    flow.varNames[1] = "T";
    flow.varNames[2] = "U";
    flow.varNames[3] = "V";
    if (input.is3D) flow.varNames[4] = "W";
    
    FlowField rhs(input.blockDim, input.blockSize, input.numVars, input.nguard);
    flow.varNames[0] = "Continuity";
    flow.varNames[1] = "Energy";
    flow.varNames[2] = "X-Momentum";
    flow.varNames[3] = "Y-Momentum";
    if (input.is3D) flow.varNames[4] = "Z-Momentum";
    
    // FillTgv(flow, input);
    // FillConst(rhs, 0.0);
    
    
    return 0;
}