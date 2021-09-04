#ifndef EXCHANGE_H
#define EXCHANGE_H
#include "CudaHeaders.h"
#include "FlowField.h"
#include "GpuConfig.h"
int GetNeighbor(FlowField& arr, int lb, int di, int dj, int dk);
void Exchange(FlowField& arr, GpuConfig& config);

#endif