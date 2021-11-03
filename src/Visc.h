#ifndef VISC_H
#define VISC_H

#include "FlowField.h"
#include "GasSpec.h"

void ComputeViscRhs(FlowField& rhs, FlowField& flow, const GasSpec& gas);

#endif