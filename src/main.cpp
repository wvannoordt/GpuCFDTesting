#include <iostream>
#include <string>
#include "InputClass.h"
#include "print.h"
#include "FlowField.h"
#include "Fill.cuh"
#include "TgvSpec.h"
#include "GasSpec.h"
#include "GpuConfig.h"
#include "Output.h"

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
    
    GasSpec gas;
    gas.Read(tree["Fluid"]);
    print(gas);
    
    TgvSpec tgv;
    tgv.Read(tree["TGV"], gas);
    print(tgv);
    
    GpuConfig config(input);
    config.Read(tree["GPU"]);
    
    FlowField flow(input.blockDim, input.blockSize, input.numVars, input.nguard, input.domainBounds);
    flow.varNames[0] = "P";
    flow.varNames[1] = "T";
    flow.varNames[2] = "U";
    flow.varNames[3] = "V";
    if (input.is3D) flow.varNames[4] = "W";
    
    FlowField rhs(input.blockDim, input.blockSize, input.numVars, input.nguard, input.domainBounds);
    rhs.varNames[0] = "Continuity";
    rhs.varNames[1] = "Energy";
    rhs.varNames[2] = "X-Momentum";
    rhs.varNames[3] = "Y-Momentum";
    if (input.is3D) rhs.varNames[4] = "Z-Momentum";
    
    FillTgv(flow, gas, tgv, config);
    Output(flow, "output", "ini");
    // FillConst(rhs, 0.0);
    
    
    return 0;
}