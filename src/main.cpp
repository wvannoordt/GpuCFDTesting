#include <iostream>
#include <string>
#include "InputClass.h"
#include "print.h"
#include "FlowField.h"
#include "Fill.h"
#include "TgvSpec.h"
#include "GasSpec.h"
#include "Output.h"
#include "Exchange.h"
#include "TimeControl.h"
#include "OutputProps.h"
#include "Conv.h"
#include "Visc.h"
#include "Advance.h"
#include "ScopeTimer.h"
#include "Mms.h"
#include "Stats.h"
#include "StatsProps.h"
#include "TimeTraceSeries.h"

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
    
    TimeControl time;
    time.Read(tree["Time"]);
    
    OutputProps outputProps;
    outputProps.Read(tree["Output"]);
    
    
    FlowField flow(input.blockDim, input.blockSize, input.numVars, input.nguard);
    flow.varNames[0] = "P";
    flow.varNames[1] = "T";
    flow.varNames[2] = "U";
    flow.varNames[3] = "V";
    if (input.is3D) flow.varNames.back() = "W";
    
    FlowField rhs(input.blockDim, input.blockSize, input.numVars, input.nguard);
    rhs.varNames[0] = "Continuity";
    rhs.varNames[1] = "Energy";
    rhs.varNames[2] = "X-Momentum";
    rhs.varNames[3] = "Y-Momentum";
    if (input.is3D) rhs.varNames.back() = "Z-Momentum";
    
    TimeTraceSeries series;
    StatsProps stats;
    stats.Read(tree["Stats"]);
    DefineStats(stats, flow, gas, series);
    
    if (outputProps.doMMS) RunMMS(flow, rhs, gas);
    
    print("Set initial condition...");
    FillTgv(flow, gas, tgv);
    
    Exchange(flow);
    Output(flow, "output", "initialCondition");
    
    for (; time.nt <= time.numSteps; time++)
    {
        ScopeTimer tmr("timestep");
        print(time);
        series.ComputeAll(time);
        FillConst(rhs, 0.0);
        ComputeConvRhs(rhs, flow, gas);
        ComputeViscRhs(rhs, flow, gas);
        if (time.nt % outputProps.outputInterval == 0)
        {
            ScopeTimer outputTimer("output");
            std::string ts = zfill(time.nt, 8);
            std::string filename = "data"+ts;
            print("Outputting file", filename);
            Output(flow, "output", filename);
        }
        Advance(rhs, flow, time.timestep, gas);
        Exchange(flow);
    }
    
    return 0;
}