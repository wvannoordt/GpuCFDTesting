#ifndef ADVANCE_H
#define ADVANCE_H

#include "FlowField.h"
#include "GasSpec.h"
#include "GpuConfig.h"

void Advance(FlowField& rhs, FlowField& flow, double timestep, const GasSpec& gas, const GpuConfig& config);

#endif