#ifndef TGV_SPEC_H
#define TGV_SPEC_H
#include "PTL.h"
#include "GasSpec.h"
#include <cmath>

struct TgvSpec
{
    double reynolds;
    double mach;
    double L;
    double rho0;
    double V0;
    double sos;
    double P0;
    
    void Read(PTL::PropertySection& section, const GasSpec& gas)
    {
        section["reynolds"].MapTo(&reynolds) = new PTL::PTLDouble(1600.0, "Specific heat ratio");
        section["mach"].MapTo(&mach)         = new PTL::PTLDouble(0.1,    "Specific heat at Constant Pressure");
        section["L"].MapTo(&L)               = new PTL::PTLDouble(1.0,    "Constant viscosity");
        section.StrictParse();
        
        sos = sqrt(gas.gamma*gas.R*gas.Tinf);
        V0 = mach*sos;
        rho0 = gas.visc*reynolds/(L*V0);
        P0 = rho0*gas.R*gas.Tinf;
    }
};

static std::ostream & operator<<(std::ostream & os, const TgvSpec & tgv)
{
   os << "TGV properties:\n";
   os << "  Re   = " << tgv.reynolds << "\n";
   os << "  M    = " << tgv.mach << "\n";
   os << "  L    = " << tgv.L << "\n";
   os << "  V0   = " << tgv.V0 << "\n";
   os << "  sos  = " << tgv.sos << "\n";
   os << "  rho0 = " << tgv.rho0 << "\n";
   os << "  P0   = " << tgv.P0 << "\n";
   return os;
}

#endif