#ifndef STATS_PROPS_H
#define STATS_PROPS_H

#include "PTL.h"

struct StatsProps
{
    bool outputKineticEnergy;
    void Read(PTL::PropertySection& section)
    {
        section["outputKineticEnergy"].MapTo(&outputKineticEnergy) = new PTL::PTLBoolean(false, "Output integrated kinetic energy");
        section.StrictParse();
    }
};

#endif