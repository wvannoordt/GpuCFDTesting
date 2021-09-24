#ifndef OUTPUT_PROPS_H
#define OUTPUT_PROPS_H
#include "PTL.h"
struct OutputProps
{
    int outputInterval;
    
    void Read(PTL::PropertySection& section)
    {
        section["outputInterval"].MapTo(&outputInterval) = new PTL::PTLInteger(100, "Interval at which to output data");
        section.StrictParse();
    }
};

#endif