#ifndef OUTPUT_PROPS_H
#define OUTPUT_PROPS_H
#include "PTL.h"
struct OutputProps
{
    int outputInterval;
    bool doMMS;
    
    void Read(PTL::PropertySection& section)
    {
        section["outputInterval"].MapTo(&outputInterval) = new PTL::PTLInteger(100, "Interval at which to output data");
        section["doMMS"].MapTo(&doMMS) = new PTL::PTLBoolean(false, "Controls method-of-manufactured-solutions mode");
        section.StrictParse();
    }
};

#endif