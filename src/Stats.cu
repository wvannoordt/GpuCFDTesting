#include "Stats.h"
#include "FlowField.h"
#include "TimeTrace.h"
#include "TimeControl.h"
#include "FlowField.h"
#include "GasSpec.h"

class IntegratedKineticEnergy : public TimeTrace
{
    public:
        IntegratedKineticEnergy(FlowField* prims_in) : TimeTrace("kineticEnergy")
        {
            prims = prims_in;
        }
        
        virtual void Compute(const TimeControl& time) final
        {
            print("compute KE [TOOO]");
        }
    private:
        FlowField* prims;
};

void DefineStats(const StatsProps& stats, FlowField& prims, GasSpec& gas, TimeTraceSeries& series)
{
    if (stats.outputKineticEnergy) series.Add<IntegratedKineticEnergy>(&prims);
}