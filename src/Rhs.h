#ifndef RHS_H
#define RHS_H

#include "FlowField.h"
#include "GasSpec.h"
#include "GpuConfig.h"

void ComputeRhs(FlowField& rhs, FlowField& flow, const GasSpec& gas, const GpuConfig& config);

#endif