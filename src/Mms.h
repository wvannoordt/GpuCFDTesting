#ifndef MMS_H
#define MMS_H

#include "FlowField.h"
#include "GasSpec.h"
#include "CudaHeaders.h"
#include "GpuConfig.h"
#include <cmath>


struct NavierStokesMms
{
    GasSpec gas;
    const double realpi = 3.1415926535;
    const double pi = 2.0*realpi;
    NavierStokesMms(const GasSpec& gas_in)
    {
        gas = gas_in;
    }

    _f_no_inline _f_hybrid double  P(const double& x, const double& y, const double& z)    const {return 5.0+2.0*cos(3.0*pi*x)*sin(2.0*pi*y)+sin(4.0*pi*z);}
    _f_no_inline _f_hybrid double dP_dx(const double& x, const double& y, const double& z) const {return -2.0*3.0*pi*sin(3.0*pi*x)*sin(2.0*pi*y);}
    _f_no_inline _f_hybrid double dP_dy(const double& x, const double& y, const double& z) const {return 2.0*2.0*pi*cos(3.0*pi*x)*cos(2.0*pi*y);}
    _f_no_inline _f_hybrid double dP_dz(const double& x, const double& y, const double& z) const {return 4.0*pi*cos(4.0*pi*z);}
    
    _f_no_inline _f_hybrid double  T(const double& x, const double& y, const double& z)    const {return 10.0+2.0*cos(2.0*pi*x)*sin(3.0*pi*y)+sin(4.0*pi*z);}
    _f_no_inline _f_hybrid double dT_dx(const double& x, const double& y, const double& z) const {return -2.0*2.0*pi*sin(2.0*pi*x)*sin(3.0*pi*y);}
    _f_no_inline _f_hybrid double dT_dy(const double& x, const double& y, const double& z) const {return  2.0*3.0*pi*cos(2.0*pi*x)*cos(3.0*pi*y);}
    _f_no_inline _f_hybrid double dT_dz(const double& x, const double& y, const double& z) const {return  4.0*pi*cos(4.0*pi*z);}
    
    _f_no_inline _f_hybrid double  U(const double& x, const double& y, const double& z)    const {return sin(3.0*pi*x)*cos(2.0*pi*y)*cos(2.0*pi*z);}
    _f_no_inline _f_hybrid double dU_dx(const double& x, const double& y, const double& z) const {return  3.0*pi*cos(3.0*pi*x)*cos(2.0*pi*y)*cos(2.0*pi*z);}
    _f_no_inline _f_hybrid double dU_dy(const double& x, const double& y, const double& z) const {return -2.0*pi*sin(3.0*pi*x)*sin(2.0*pi*y)*cos(2.0*pi*z);}
    _f_no_inline _f_hybrid double dU_dz(const double& x, const double& y, const double& z) const {return -2.0*pi*sin(3.0*pi*x)*cos(2.0*pi*y)*sin(2.0*pi*z);}
    
    _f_no_inline _f_hybrid double  V(const double& x, const double& y, const double& z)    const {return cos(3.0*pi*x)*cos(2.0*pi*y)*cos(3.0*pi*z);}
    _f_no_inline _f_hybrid double dV_dx(const double& x, const double& y, const double& z) const {return -3.0*pi*sin(3.0*pi*x)*cos(2.0*pi*y)*cos(3.0*pi*z);}
    _f_no_inline _f_hybrid double dV_dy(const double& x, const double& y, const double& z) const {return -2.0*pi*cos(3.0*pi*x)*sin(2.0*pi*y)*cos(3.0*pi*z);}
    _f_no_inline _f_hybrid double dV_dz(const double& x, const double& y, const double& z) const {return -3.0*pi*cos(3.0*pi*x)*cos(2.0*pi*y)*sin(3.0*pi*z);}
    
    _f_no_inline _f_hybrid double  W(const double& x, const double& y, const double& z)    const {return sin(3.0*pi*x)*sin(2.0*pi*y)*cos(4.0*pi*z);}
    _f_no_inline _f_hybrid double dW_dx(const double& x, const double& y, const double& z) const {return  3.0*pi*cos(3.0*pi*x)*sin(2.0*pi*y)*cos(4.0*pi*z);}
    _f_no_inline _f_hybrid double dW_dy(const double& x, const double& y, const double& z) const {return  2.0*pi*sin(3.0*pi*x)*cos(2.0*pi*y)*cos(4.0*pi*z);}
    _f_no_inline _f_hybrid double dW_dz(const double& x, const double& y, const double& z) const {return -4.0*pi*sin(3.0*pi*x)*sin(2.0*pi*y)*sin(4.0*pi*z);}
    
    _f_no_inline _f_hybrid double  Rho(const double& x, const double& y, const double& z) const {return P(x,y,z)/(gas.R*T(x,y,z));}
    _f_no_inline _f_hybrid double dRho_dx(const double& x, const double& y, const double& z) const
    {
        return (T(x,y,z)*dP_dx(x,y,z)-P(x,y,z)*dT_dx(x,y,z))/(T(x,y,z)*T(x,y,z)*gas.R);
    }
    
    _f_no_inline _f_hybrid double dRho_dy(const double& x, const double& y, const double& z) const
    {
        return (T(x,y,z)*dP_dy(x,y,z)-P(x,y,z)*dT_dy(x,y,z))/(T(x,y,z)*T(x,y,z)*gas.R);
    }
    
    _f_no_inline _f_hybrid double dRho_dz(const double& x, const double& y, const double& z) const
    {
        return (T(x,y,z)*dP_dz(x,y,z)-P(x,y,z)*dT_dz(x,y,z))/(T(x,y,z)*T(x,y,z)*gas.R);
    }
    
    _f_no_inline _f_hybrid double  H(const double& x, const double& y, const double& z) const
    {
        return 0.5*(U(x,y,z)*U(x,y,z) + V(x,y,z)*V(x,y,z) + W(x,y,z)*W(x,y,z)) + gas.gamma*P(x,y,z)/((gas.gamma-1.0)*Rho(x,y,z));
    }
    _f_no_inline _f_hybrid double dH_dx(const double& x, const double& y, const double& z) const
    {
        return U(x,y,z)*dU_dx(x,y,z)+V(x,y,z)*dV_dx(x,y,z)+W(x,y,z)*dW_dx(x,y,z)
            + gas.gamma*(Rho(x,y,z)*dP_dx(x,y,z)-dRho_dx(x,y,z)*P(x,y,z))/(Rho(x,y,z)*Rho(x,y,z)*(gas.gamma-1.0));
    }
    _f_no_inline _f_hybrid double dH_dy(const double& x, const double& y, const double& z) const
    {
        return U(x,y,z)*dU_dy(x,y,z)+V(x,y,z)*dV_dy(x,y,z)+W(x,y,z)*dW_dy(x,y,z)
            + gas.gamma*(Rho(x,y,z)*dP_dy(x,y,z)-dRho_dy(x,y,z)*P(x,y,z))/(Rho(x,y,z)*Rho(x,y,z)*(gas.gamma-1.0));
    }
    _f_no_inline _f_hybrid double dH_dz(const double& x, const double& y, const double& z) const
    {
        return U(x,y,z)*dU_dz(x,y,z)+V(x,y,z)*dV_dz(x,y,z)+W(x,y,z)*dW_dz(x,y,z)
            + gas.gamma*(Rho(x,y,z)*dP_dz(x,y,z)-dRho_dz(x,y,z)*P(x,y,z))/(Rho(x,y,z)*Rho(x,y,z)*(gas.gamma-1.0));
    }
    
    _f_no_inline _f_hybrid void conv_rhs(const double& x, const double& y, const double& z, double (&rhsAr)[5]) const
    {
        rhsAr[0] = 
            -U(x,y,z)*dRho_dx(x,y,z)-dU_dx(x,y,z)*Rho(x,y,z)
            -V(x,y,z)*dRho_dy(x,y,z)-dV_dy(x,y,z)*Rho(x,y,z)
            -W(x,y,z)*dRho_dz(x,y,z)-dW_dz(x,y,z)*Rho(x,y,z);
        rhsAr[1] = 
            -dRho_dx(x,y,z)*U    (x,y,z)*H    (x,y,z)
            -Rho    (x,y,z)*dU_dx(x,y,z)*H    (x,y,z)
            -Rho    (x,y,z)*U    (x,y,z)*dH_dx(x,y,z)
            -dRho_dy(x,y,z)*V    (x,y,z)*H    (x,y,z)
            -Rho    (x,y,z)*dV_dy(x,y,z)*H    (x,y,z)
            -Rho    (x,y,z)*V    (x,y,z)*dH_dy(x,y,z)
            -dRho_dz(x,y,z)*W    (x,y,z)*H    (x,y,z)
            -Rho    (x,y,z)*dW_dz(x,y,z)*H    (x,y,z)
            -Rho    (x,y,z)*W    (x,y,z)*dH_dz(x,y,z);
        rhsAr[2] = 
            -dRho_dx(x,y,z)*U    (x,y,z)*U    (x,y,z)
            -Rho    (x,y,z)*dU_dx(x,y,z)*U    (x,y,z)
            -Rho    (x,y,z)*U    (x,y,z)*dU_dx(x,y,z)
            -dRho_dy(x,y,z)*V    (x,y,z)*U    (x,y,z)
            -Rho    (x,y,z)*dV_dy(x,y,z)*U    (x,y,z)
            -Rho    (x,y,z)*V    (x,y,z)*dU_dy(x,y,z)
            -dRho_dz(x,y,z)*W    (x,y,z)*U    (x,y,z)
            -Rho    (x,y,z)*dW_dz(x,y,z)*U    (x,y,z)
            -Rho    (x,y,z)*W    (x,y,z)*dU_dz(x,y,z)
            -dP_dx  (x,y,z);
        rhsAr[3] = 
            -dRho_dx(x,y,z)*U    (x,y,z)*V    (x,y,z)
            -Rho    (x,y,z)*dU_dx(x,y,z)*V    (x,y,z)
            -Rho    (x,y,z)*U    (x,y,z)*dV_dx(x,y,z)
            -dRho_dy(x,y,z)*V    (x,y,z)*V    (x,y,z)
            -Rho    (x,y,z)*dV_dy(x,y,z)*V    (x,y,z)
            -Rho    (x,y,z)*V    (x,y,z)*dV_dy(x,y,z)
            -dRho_dz(x,y,z)*W    (x,y,z)*V    (x,y,z)
            -Rho    (x,y,z)*dW_dz(x,y,z)*V    (x,y,z)
            -Rho    (x,y,z)*W    (x,y,z)*dV_dz(x,y,z)
            -dP_dy  (x,y,z);
        rhsAr[4] = 
            -dRho_dx(x,y,z)*U    (x,y,z)*W    (x,y,z)
            -Rho    (x,y,z)*dU_dx(x,y,z)*W    (x,y,z)
            -Rho    (x,y,z)*U    (x,y,z)*dW_dx(x,y,z)
            -dRho_dy(x,y,z)*V    (x,y,z)*W    (x,y,z)
            -Rho    (x,y,z)*dV_dy(x,y,z)*W    (x,y,z)
            -Rho    (x,y,z)*V    (x,y,z)*dW_dy(x,y,z)
            -dRho_dz(x,y,z)*W    (x,y,z)*W    (x,y,z)
            -Rho    (x,y,z)*dW_dz(x,y,z)*W    (x,y,z)
            -Rho    (x,y,z)*W    (x,y,z)*dW_dz(x,y,z)
            -dP_dz  (x,y,z);
    }
    
    _f_hybrid void testFcn(const double& x, const double& y, const double& z, double (&primsAr)[5]) const
    {
        primsAr[0] = P(x, y, z);
        primsAr[1] = T(x, y, z);
        primsAr[2] = U(x, y, z);
        primsAr[3] = V(x, y, z);
        primsAr[4] = W(x, y, z);
    }
};

void AnalyticalConvRhs(const NavierStokesMms& mms, FlowField& rhs, const GpuConfig& config);
void AnalyticalFcn(const NavierStokesMms& mms, FlowField& prims, const GpuConfig& config);
void RunMMS(FlowField& prims, FlowField& rhs, const GasSpec& gas, const GpuConfig& config);
#endif