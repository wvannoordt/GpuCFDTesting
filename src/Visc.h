#ifndef VISC_H
#define VISC_H

#include "FlowField.h"
#include "GasSpec.h"
#include "GpuConfig.h"

void ComputeViscRhs(FlowField& rhs, FlowField& flow, const GasSpec& gas, const GpuConfig& config);

#endif