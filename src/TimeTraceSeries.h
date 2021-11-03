#ifndef TIME_TRACE_SERIES_H
#define TIME_TRACE_SERIES_H

#include <vector>
#include "TimeControl.h"
#include "TimeTrace.h"
#include "BaseClassContainer.h"

class TimeTraceSeries : public BaseClassContainer<TimeTrace>
{
    public:
        
        TimeTraceSeries()
        {
            
        }
        
        void ComputeAll(const TimeControl& time)
        {
            for (auto& trace: *this)
            {
                trace->Compute(time);
            }
        }
};

#endif