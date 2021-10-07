#ifndef METRIC_H
#define METRIC_H

#include "CudaHeaders.h"
#include "v3.h"
#include "m9.h"
static _f_hybrid inline void GetCoords(v3<double>& eta, v3<double>& xyz) // \eta(x)
{
    xyz[0] = eta[0];
    xyz[0] = (0.5*eta[0]*eta[0]-0.3333333333333*eta[0]*eta[0]*eta[0]) + 0.05*eta[0];
    xyz[1] = eta[1];
    xyz[2] = eta[2];
};

static _f_hybrid inline void GetCoordsGrad(v3<double>& eta, m9<double>& deta_dxyz) // \eta(x)
{
    //eta, x
    deta_dxyz(0,0) = 6.0*(eta[0]-eta[0]*eta[0]);
    deta_dxyz(0,0) = 1.0;
    deta_dxyz(0,1) = 0.0;
    deta_dxyz(0,2) = 0.0;
    
    deta_dxyz(1,0) = 0.0;
    deta_dxyz(1,1) = 1.0;
    deta_dxyz(1,2) = 0.0;
    
    deta_dxyz(2,0) = 0.0;
    deta_dxyz(2,1) = 0.0;
    deta_dxyz(2,2) = 1.0;
};

#endif