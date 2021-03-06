#ifndef GAS_SPEC_H
#define GAS_SPEC_H
#include "PTL.h"
struct GasSpec
{
    double gamma;
    double cp;
    double visc;
    double Tinf;
    double beta; // bulk viscosity
    double R;
    double prandtl;
    
    void Read(PTL::PropertySection& section)
    {
        section["gamma"].MapTo(&gamma) = new PTL::PTLDouble(1.4, "Specific heat ratio");
        section["cp"].MapTo(&cp) = new PTL::PTLDouble(1005.0, "Specific heat at Constant Pressure");
        section["visc"].MapTo(&visc) = new PTL::PTLDouble(1.8e-5, "Constant viscosity");
        section["Tinf"].MapTo(&Tinf) = new PTL::PTLDouble(300.0, "Freestream temperature");
        section["prandtl"].MapTo(&prandtl) = new PTL::PTLDouble(0.72, "Prandtl number");
        section.StrictParse();
        
        R = cp*(gamma-1)/gamma;
        beta = (-2.0/3.0)*visc; // stokes hypothesis
    }
};

static std::ostream & operator<<(std::ostream & os, const GasSpec & gas)
{
   os << "Gas properties:\n  gamma = " << gas.gamma << "\n  cp    = " << gas.cp << "\n  visc  = " << gas.visc << "\n  Tinf  = " << gas.Tinf;
   return os;
}

#endif