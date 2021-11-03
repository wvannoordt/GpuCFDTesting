#ifndef STATS_H
#define STATS_H

#include "StatsProps.h"
#include "TimeTraceSeries.h"
#include "FlowField.h"
#include "GasSpec.h"

void DefineStats(const StatsProps& stats, FlowField& prims, GasSpec& gas, TimeTraceSeries& series);

#endif