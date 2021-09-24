#ifndef TIME_CONTROL_H
#define TIME_CONTROL_H

#include "PTL.h"

struct TimeControl
{
    int numSteps;
    double timestep;
    double time = 0.0;
    int nt = 0;
    
    void Read(PTL::PropertySection& section)
    {
        section["numSteps"].MapTo(&numSteps) = new PTL::PTLInteger(1000, "Number of timesteps to take");
        section["timestep"].MapTo(&timestep) = new PTL::PTLDouble(1e-4, "Physical timestep");
        section.StrictParse();
    }
    
    inline TimeControl& operator ++ (int dummy)
    {
        nt++;
        time+=timestep;
        return *this;
    }
    
};

static std::ostream & operator<<(std::ostream & os, const TimeControl & time)
{
    os << "Timestep: " << time.nt << ", time: " << time.time;
    return os;
}

#endif