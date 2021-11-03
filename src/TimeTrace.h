#ifndef TRACE_H
#define TRACE_H

#include "TimeControl.h"
#include <string>

class TimeTrace
{
    public:
        TimeTrace(std::string name_in)
        {
            name = name_in;
        }
        virtual void Compute(const TimeControl& time)=0;
        std::string GetName(void) {return name;}
    private:
        std::string name;
};

#endif