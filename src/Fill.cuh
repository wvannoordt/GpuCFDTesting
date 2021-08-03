#ifndef FILL_CUH
#define FILL_CUH

#include "FlowField.h"
#include "InputClass.h"
#include "TgvSpec.h"
#include "GasSpec.h"
#include "GpuConfig.h"
void FillTgv(FlowField& flow, const GasSpec& gas, const TgvSpec& tgv, const GpuConfig& config);
// void FillConst(FlowField& arr,  double val);

#endif