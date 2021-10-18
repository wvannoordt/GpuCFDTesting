#ifndef METRIC_H
#define METRIC_H

#include "CudaHeaders.h"
#include "v3.h"
#include "m9.h"

const double pi = 3.1415926535;
const double dxWall = 0.01;

static _f_hybrid inline void GetCoords(v3<double>& eta, v3<double>& xyz) // \eta(x)
{
    auto f = [&](double d) -> double {return 0.5-0.5*cos(pi*d);};
    double delta = (1.0-eta[0])*(-f(dxWall)) + eta[0]*f(dxWall);
    xyz[0] = f((1.0-eta[0])*(dxWall)+(eta[0])*(1.0-dxWall))+delta;
    xyz[0] = eta[0]*eta[0] + 0.1*eta[0];
    xyz[1] = eta[1];
    xyz[2] = eta[2];
};

static _f_hybrid inline void GetCoordsGrad(v3<double>& eta, m9<double>& deta_dxyz) // \eta(x)
{
    //xi
    deta_dxyz(0,0) = 1.0/(2.0*eta[0]+0.1);
    deta_dxyz(0,1) = 0.0;
    deta_dxyz(0,2) = 0.0;
    
    //eta
    deta_dxyz(1,0) = 0.0;
    deta_dxyz(1,1) = 1.0;
    deta_dxyz(1,2) = 0.0;
    
    //zeta
    deta_dxyz(2,0) = 0.0;
    deta_dxyz(2,1) = 0.0;
    deta_dxyz(2,2) = 1.0;
};

#endif