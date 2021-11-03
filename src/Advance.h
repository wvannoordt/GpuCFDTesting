#ifndef ADVANCE_H
#define ADVANCE_H

#include "FlowField.h"
#include "GasSpec.h"

void Advance(FlowField& rhs, FlowField& flow, double timestep, const GasSpec& gas);

#endif