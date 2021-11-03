#ifndef CONV_H
#define CONV_H

#include "FlowField.h"
#include "GasSpec.h"

void ComputeConvRhs(FlowField& rhs, FlowField& flow, const GasSpec& gas);

#endif