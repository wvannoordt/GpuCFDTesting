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
    input.Read(tree["CFD"]);
    
    return 0;
}