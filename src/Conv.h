#ifndef CONV_H
#define CONV_H

#include "FlowField.h"
#include "GasSpec.h"
#include "GpuConfig.h"

void ComputeConvRhs(FlowField& rhs, FlowField& flow, const GasSpec& gas, const GpuConfig& config);

#endif