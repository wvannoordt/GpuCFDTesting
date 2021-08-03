#ifndef TGV_SPEC_H
#define TGV_SPEC_H
#include "PTL.h"
struct TgvSpec
{
    double reynolds;
    double mach;
    double L;
    
    void Read(PTL::PropertySection& section)
    {
        section["reynolds"].MapTo(&reynolds) = new PTL::PTLDouble(1600.0, "Specific heat ratio");
        section["mach"].MapTo(&mach)         = new PTL::PTLDouble(0.1,    "Specific heat at Constant Pressure");
        section["L"].MapTo(&L)               = new PTL::PTLDouble(1.0,    "Constant viscosity");
        section.StrictParse();
    }
};

static std::ostream & operator<<(std::ostream & os, const TgvSpec & tgv)
{
   os << "TGV properties:\n  Re = " << tgv.reynolds << "\n  M  = " << tgv.mach << "\n  L  = " << tgv.L;
   return os;
}

#endif