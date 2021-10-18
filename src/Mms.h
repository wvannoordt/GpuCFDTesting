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

    //P and derivatives
    _f_no_inline _f_hybrid double  P(const double& x, const double& y, const double& z)    const {return 5.0+2.0*cos(3.0*pi*x)*sin(2.0*pi*y)+sin(4.0*pi*z);}
    _f_no_inline _f_hybrid double dP_dx(const double& x, const double& y, const double& z) const {return -2.0*3.0*pi*sin(3.0*pi*x)*sin(2.0*pi*y);}
    _f_no_inline _f_hybrid double dP_dy(const double& x, const double& y, const double& z) const {return 2.0*2.0*pi*cos(3.0*pi*x)*cos(2.0*pi*y);}
    _f_no_inline _f_hybrid double dP_dz(const double& x, const double& y, const double& z) const {return 4.0*pi*cos(4.0*pi*z);}

    //T and derivatives
    _f_no_inline _f_hybrid double  T(const double& x, const double& y, const double& z)    const {return 10.0+2.0*cos(2.0*pi*x)*sin(3.0*pi*y)+sin(4.0*pi*z);}
    _f_no_inline _f_hybrid double dT_dx(const double& x, const double& y, const double& z) const {return -2.0*2.0*pi*sin(2.0*pi*x)*sin(3.0*pi*y);}
    _f_no_inline _f_hybrid double dT_dy(const double& x, const double& y, const double& z) const {return  2.0*3.0*pi*cos(2.0*pi*x)*cos(3.0*pi*y);}
    _f_no_inline _f_hybrid double dT_dz(const double& x, const double& y, const double& z) const {return  4.0*pi*cos(4.0*pi*z);}
    _f_no_inline _f_hybrid double d2T_dx2(const double& x, const double& y, const double& z) const {return -2.0*2.0*2.0*pi*pi*cos(2.0*pi*x)*sin(3.0*pi*y);}
    _f_no_inline _f_hybrid double d2T_dy2(const double& x, const double& y, const double& z) const {return  -2.0*3.0*3.0*pi*pi*cos(2.0*pi*x)*sin(3.0*pi*y);}
    _f_no_inline _f_hybrid double d2T_dz2(const double& x, const double& y, const double& z) const {return  -4.0*4.0*pi*pi*sin(4.0*pi*z);}
    
    //U and derivatives
    _f_no_inline _f_hybrid double  U(const double& x, const double& y, const double& z)    const {return sin(3.0*pi*x)*cos(2.0*pi*y)*cos(2.0*pi*z);}
    _f_no_inline _f_hybrid double dU_dx(const double& x, const double& y, const double& z) const {return  3.0*pi*cos(3.0*pi*x)*cos(2.0*pi*y)*cos(2.0*pi*z);}
    _f_no_inline _f_hybrid double dU_dy(const double& x, const double& y, const double& z) const {return -2.0*pi*sin(3.0*pi*x)*sin(2.0*pi*y)*cos(2.0*pi*z);}
    _f_no_inline _f_hybrid double dU_dz(const double& x, const double& y, const double& z) const {return -2.0*pi*sin(3.0*pi*x)*cos(2.0*pi*y)*sin(2.0*pi*z);}
    
    _f_no_inline _f_hybrid double d2U_dxx(const double& x, const double& y, const double& z) const {return  -3.0*3.0*pi*pi*sin(3.0*pi*x)*cos(2.0*pi*y)*cos(2.0*pi*z);}
    _f_no_inline _f_hybrid double d2U_dxy(const double& x, const double& y, const double& z) const {return  -3.0*2.0*pi*pi*cos(3.0*pi*x)*sin(2.0*pi*y)*cos(2.0*pi*z);}
    _f_no_inline _f_hybrid double d2U_dxz(const double& x, const double& y, const double& z) const {return  -3.0*2.0*pi*pi*cos(3.0*pi*x)*cos(2.0*pi*y)*sin(2.0*pi*z);}
    
    _f_no_inline _f_hybrid double d2U_dyx(const double& x, const double& y, const double& z) const {return -2.0*3.0*pi*pi*cos(3.0*pi*x)*sin(2.0*pi*y)*cos(2.0*pi*z);}
    _f_no_inline _f_hybrid double d2U_dyy(const double& x, const double& y, const double& z) const {return -2.0*2.0*pi*pi*sin(3.0*pi*x)*cos(2.0*pi*y)*cos(2.0*pi*z);}
    _f_no_inline _f_hybrid double d2U_dyz(const double& x, const double& y, const double& z) const {return  2.0*2.0*pi*pi*sin(3.0*pi*x)*sin(2.0*pi*y)*sin(2.0*pi*z);}
    
    _f_no_inline _f_hybrid double d2U_dzx(const double& x, const double& y, const double& z) const {return -2.0*3.0*pi*pi*cos(3.0*pi*x)*cos(2.0*pi*y)*sin(2.0*pi*z);}
    _f_no_inline _f_hybrid double d2U_dzy(const double& x, const double& y, const double& z) const {return  2.0*2.0*pi*pi*sin(3.0*pi*x)*sin(2.0*pi*y)*sin(2.0*pi*z);}
    _f_no_inline _f_hybrid double d2U_dzz(const double& x, const double& y, const double& z) const {return -2.0*2.0*pi*pi*sin(3.0*pi*x)*cos(2.0*pi*y)*cos(2.0*pi*z);}
    
    //V and derivatives
    _f_no_inline _f_hybrid double  V(const double& x, const double& y, const double& z)    const {return cos(3.0*pi*x)*cos(2.0*pi*y)*cos(3.0*pi*z);}
    _f_no_inline _f_hybrid double dV_dx(const double& x, const double& y, const double& z) const {return -3.0*pi*sin(3.0*pi*x)*cos(2.0*pi*y)*cos(3.0*pi*z);}
    _f_no_inline _f_hybrid double dV_dy(const double& x, const double& y, const double& z) const {return -2.0*pi*cos(3.0*pi*x)*sin(2.0*pi*y)*cos(3.0*pi*z);}
    _f_no_inline _f_hybrid double dV_dz(const double& x, const double& y, const double& z) const {return -3.0*pi*cos(3.0*pi*x)*cos(2.0*pi*y)*sin(3.0*pi*z);}
    
    _f_no_inline _f_hybrid double d2V_dxx(const double& x, const double& y, const double& z) const {return -3.0*3.0*pi*pi*cos(3.0*pi*x)*cos(2.0*pi*y)*cos(3.0*pi*z);}
    _f_no_inline _f_hybrid double d2V_dxy(const double& x, const double& y, const double& z) const {return  3.0*2.0*pi*pi*sin(3.0*pi*x)*sin(2.0*pi*y)*cos(3.0*pi*z);}
    _f_no_inline _f_hybrid double d2V_dxz(const double& x, const double& y, const double& z) const {return  3.0*3.0*pi*pi*sin(3.0*pi*x)*cos(2.0*pi*y)*sin(3.0*pi*z);}
    
    _f_no_inline _f_hybrid double d2V_dyx(const double& x, const double& y, const double& z) const {return  2.0*3.0*pi*pi*sin(3.0*pi*x)*sin(2.0*pi*y)*cos(3.0*pi*z);}
    _f_no_inline _f_hybrid double d2V_dyy(const double& x, const double& y, const double& z) const {return -2.0*2.0*pi*pi*cos(3.0*pi*x)*cos(2.0*pi*y)*cos(3.0*pi*z);}
    _f_no_inline _f_hybrid double d2V_dyz(const double& x, const double& y, const double& z) const {return  2.0*3.0*pi*pi*cos(3.0*pi*x)*sin(2.0*pi*y)*sin(3.0*pi*z);}
    
    _f_no_inline _f_hybrid double d2V_dzx(const double& x, const double& y, const double& z) const {return  3.0*3.0*pi*pi*sin(3.0*pi*x)*cos(2.0*pi*y)*sin(3.0*pi*z);}
    _f_no_inline _f_hybrid double d2V_dzy(const double& x, const double& y, const double& z) const {return  3.0*2.0*pi*pi*cos(3.0*pi*x)*sin(2.0*pi*y)*sin(3.0*pi*z);}
    _f_no_inline _f_hybrid double d2V_dzz(const double& x, const double& y, const double& z) const {return -3.0*3.0*pi*pi*cos(3.0*pi*x)*cos(2.0*pi*y)*cos(3.0*pi*z);}
    
    //W and derivatives
    _f_no_inline _f_hybrid double  W(const double& x, const double& y, const double& z)    const {return sin(3.0*pi*x)*sin(2.0*pi*y)*cos(4.0*pi*z);}
    _f_no_inline _f_hybrid double dW_dx(const double& x, const double& y, const double& z) const {return  3.0*pi*cos(3.0*pi*x)*sin(2.0*pi*y)*cos(4.0*pi*z);}
    _f_no_inline _f_hybrid double dW_dy(const double& x, const double& y, const double& z) const {return  2.0*pi*sin(3.0*pi*x)*cos(2.0*pi*y)*cos(4.0*pi*z);}
    _f_no_inline _f_hybrid double dW_dz(const double& x, const double& y, const double& z) const {return -4.0*pi*sin(3.0*pi*x)*sin(2.0*pi*y)*sin(4.0*pi*z);}
    
    _f_no_inline _f_hybrid double d2W_dxx(const double& x, const double& y, const double& z) const {return  -3.0*3.0*pi*pi*sin(3.0*pi*x)*sin(2.0*pi*y)*cos(4.0*pi*z);}
    _f_no_inline _f_hybrid double d2W_dxy(const double& x, const double& y, const double& z) const {return   3.0*2.0*pi*pi*cos(3.0*pi*x)*cos(2.0*pi*y)*cos(4.0*pi*z);}
    _f_no_inline _f_hybrid double d2W_dxz(const double& x, const double& y, const double& z) const {return  -3.0*4.0*pi*pi*cos(3.0*pi*x)*sin(2.0*pi*y)*sin(4.0*pi*z);}
    
    _f_no_inline _f_hybrid double d2W_dyx(const double& x, const double& y, const double& z) const {return   2.0*3.0*pi*pi*cos(3.0*pi*x)*cos(2.0*pi*y)*cos(4.0*pi*z);}
    _f_no_inline _f_hybrid double d2W_dyy(const double& x, const double& y, const double& z) const {return  -2.0*2.0*pi*pi*sin(3.0*pi*x)*sin(2.0*pi*y)*cos(4.0*pi*z);}
    _f_no_inline _f_hybrid double d2W_dyz(const double& x, const double& y, const double& z) const {return  -2.0*4.0*pi*pi*sin(3.0*pi*x)*cos(2.0*pi*y)*sin(4.0*pi*z);}
    
    _f_no_inline _f_hybrid double d2W_dzx(const double& x, const double& y, const double& z) const {return -4.0*3.0*pi*pi*cos(3.0*pi*x)*sin(2.0*pi*y)*sin(4.0*pi*z);}
    _f_no_inline _f_hybrid double d2W_dzy(const double& x, const double& y, const double& z) const {return -4.0*2.0*pi*pi*sin(3.0*pi*x)*cos(2.0*pi*y)*sin(4.0*pi*z);}
    _f_no_inline _f_hybrid double d2W_dzz(const double& x, const double& y, const double& z) const {return -4.0*4.0*pi*pi*sin(3.0*pi*x)*sin(2.0*pi*y)*cos(4.0*pi*z);}
    
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
    
    _f_no_inline _f_hybrid double Tau_xx(const double& x, const double& y, const double& z) const { return gas.visc*(dU_dx(x,y,z)+dU_dx(x,y,z))+gas.beta*(dU_dx(x,y,z)+dV_dy(x,y,z)+dW_dz(x,y,z)); }
    _f_no_inline _f_hybrid double Tau_xy(const double& x, const double& y, const double& z) const { return gas.visc*(dU_dy(x,y,z)+dV_dx(x,y,z)); }
    _f_no_inline _f_hybrid double Tau_xz(const double& x, const double& y, const double& z) const { return gas.visc*(dU_dz(x,y,z)+dW_dx(x,y,z)); }
    _f_no_inline _f_hybrid double Tau_yx(const double& x, const double& y, const double& z) const { return gas.visc*(dV_dx(x,y,z)+dU_dy(x,y,z)); }
    _f_no_inline _f_hybrid double Tau_yy(const double& x, const double& y, const double& z) const { return gas.visc*(dV_dy(x,y,z)+dV_dy(x,y,z))+gas.beta*(dU_dx(x,y,z)+dV_dy(x,y,z)+dW_dz(x,y,z)); }
    _f_no_inline _f_hybrid double Tau_yz(const double& x, const double& y, const double& z) const { return gas.visc*(dV_dz(x,y,z)+dW_dy(x,y,z)); }
    _f_no_inline _f_hybrid double Tau_zx(const double& x, const double& y, const double& z) const { return gas.visc*(dW_dx(x,y,z)+dU_dz(x,y,z)); }
    _f_no_inline _f_hybrid double Tau_zy(const double& x, const double& y, const double& z) const { return gas.visc*(dW_dy(x,y,z)+dV_dz(x,y,z)); }
    _f_no_inline _f_hybrid double Tau_zz(const double& x, const double& y, const double& z) const { return gas.visc*(dW_dz(x,y,z)+dW_dz(x,y,z))+gas.beta*(dU_dx(x,y,z)+dV_dy(x,y,z)+dW_dz(x,y,z)); }
    
    _f_no_inline _f_hybrid double dTau_xx_dx(const double& x, const double& y, const double& z) const { return gas.visc*(d2U_dxx(x,y,z)+d2U_dxx(x,y,z))+gas.beta*(d2U_dxx(x,y,z)+d2V_dyx(x,y,z)+d2W_dzx(x,y,z)); }
    _f_no_inline _f_hybrid double dTau_xx_dy(const double& x, const double& y, const double& z) const { return gas.visc*(d2U_dxy(x,y,z)+d2U_dxy(x,y,z))+gas.beta*(d2U_dxy(x,y,z)+d2V_dyy(x,y,z)+d2W_dzy(x,y,z)); }
    _f_no_inline _f_hybrid double dTau_xx_dz(const double& x, const double& y, const double& z) const { return gas.visc*(d2U_dxz(x,y,z)+d2U_dxz(x,y,z))+gas.beta*(d2U_dxz(x,y,z)+d2V_dyz(x,y,z)+d2W_dzz(x,y,z)); }
    
    _f_no_inline _f_hybrid double dTau_xy_dx(const double& x, const double& y, const double& z) const { return gas.visc*(d2U_dyx(x,y,z)+d2V_dxx(x,y,z)); }
    _f_no_inline _f_hybrid double dTau_xy_dy(const double& x, const double& y, const double& z) const { return gas.visc*(d2U_dyy(x,y,z)+d2V_dxy(x,y,z)); }
    _f_no_inline _f_hybrid double dTau_xy_dz(const double& x, const double& y, const double& z) const { return gas.visc*(d2U_dyz(x,y,z)+d2V_dxz(x,y,z)); }
    
    _f_no_inline _f_hybrid double dTau_xz_dx(const double& x, const double& y, const double& z) const { return gas.visc*(d2U_dzx(x,y,z)+d2W_dxx(x,y,z)); }
    _f_no_inline _f_hybrid double dTau_xz_dy(const double& x, const double& y, const double& z) const { return gas.visc*(d2U_dzy(x,y,z)+d2W_dxy(x,y,z)); }
    _f_no_inline _f_hybrid double dTau_xz_dz(const double& x, const double& y, const double& z) const { return gas.visc*(d2U_dzz(x,y,z)+d2W_dxz(x,y,z)); }
    
    _f_no_inline _f_hybrid double dTau_yx_dx(const double& x, const double& y, const double& z) const { return gas.visc*(d2V_dxx(x,y,z)+d2U_dyx(x,y,z)); }
    _f_no_inline _f_hybrid double dTau_yx_dy(const double& x, const double& y, const double& z) const { return gas.visc*(d2V_dxy(x,y,z)+d2U_dyy(x,y,z)); }
    _f_no_inline _f_hybrid double dTau_yx_dz(const double& x, const double& y, const double& z) const { return gas.visc*(d2V_dxz(x,y,z)+d2U_dyz(x,y,z)); }
    
    _f_no_inline _f_hybrid double dTau_yy_dx(const double& x, const double& y, const double& z) const { return gas.visc*(d2V_dyx(x,y,z)+d2V_dyx(x,y,z))+gas.beta*(d2U_dxx(x,y,z)+d2V_dyx(x,y,z)+d2W_dzx(x,y,z)); }
    _f_no_inline _f_hybrid double dTau_yy_dy(const double& x, const double& y, const double& z) const { return gas.visc*(d2V_dyy(x,y,z)+d2V_dyy(x,y,z))+gas.beta*(d2U_dxy(x,y,z)+d2V_dyy(x,y,z)+d2W_dzy(x,y,z)); }
    _f_no_inline _f_hybrid double dTau_yy_dz(const double& x, const double& y, const double& z) const { return gas.visc*(d2V_dyz(x,y,z)+d2V_dyz(x,y,z))+gas.beta*(d2U_dxz(x,y,z)+d2V_dyz(x,y,z)+d2W_dzz(x,y,z)); }
    
    _f_no_inline _f_hybrid double dTau_yz_dx(const double& x, const double& y, const double& z) const { return gas.visc*(d2V_dzx(x,y,z)+d2W_dyx(x,y,z)); }
    _f_no_inline _f_hybrid double dTau_yz_dy(const double& x, const double& y, const double& z) const { return gas.visc*(d2V_dzy(x,y,z)+d2W_dyy(x,y,z)); }
    _f_no_inline _f_hybrid double dTau_yz_dz(const double& x, const double& y, const double& z) const { return gas.visc*(d2V_dzz(x,y,z)+d2W_dyz(x,y,z)); }
    
    _f_no_inline _f_hybrid double dTau_zx_dx(const double& x, const double& y, const double& z) const { return gas.visc*(d2W_dxx(x,y,z)+d2U_dzx(x,y,z)); }
    _f_no_inline _f_hybrid double dTau_zx_dy(const double& x, const double& y, const double& z) const { return gas.visc*(d2W_dxy(x,y,z)+d2U_dzy(x,y,z)); }
    _f_no_inline _f_hybrid double dTau_zx_dz(const double& x, const double& y, const double& z) const { return gas.visc*(d2W_dxz(x,y,z)+d2U_dzz(x,y,z)); }
    
    _f_no_inline _f_hybrid double dTau_zy_dx(const double& x, const double& y, const double& z) const { return gas.visc*(d2W_dyx(x,y,z)+d2V_dzx(x,y,z)); }
    _f_no_inline _f_hybrid double dTau_zy_dy(const double& x, const double& y, const double& z) const { return gas.visc*(d2W_dyy(x,y,z)+d2V_dzy(x,y,z)); }
    _f_no_inline _f_hybrid double dTau_zy_dz(const double& x, const double& y, const double& z) const { return gas.visc*(d2W_dyz(x,y,z)+d2V_dzz(x,y,z)); }
    
    _f_no_inline _f_hybrid double dTau_zz_dx(const double& x, const double& y, const double& z) const { return gas.visc*(d2W_dzx(x,y,z)+d2W_dzx(x,y,z))+gas.beta*(d2U_dxx(x,y,z)+d2V_dyx(x,y,z)+d2W_dzx(x,y,z)); }
    _f_no_inline _f_hybrid double dTau_zz_dy(const double& x, const double& y, const double& z) const { return gas.visc*(d2W_dzy(x,y,z)+d2W_dzy(x,y,z))+gas.beta*(d2U_dxy(x,y,z)+d2V_dyy(x,y,z)+d2W_dzy(x,y,z)); }
    _f_no_inline _f_hybrid double dTau_zz_dz(const double& x, const double& y, const double& z) const { return gas.visc*(d2W_dzz(x,y,z)+d2W_dzz(x,y,z))+gas.beta*(d2U_dxz(x,y,z)+d2V_dyz(x,y,z)+d2W_dzz(x,y,z)); }
    
    _f_no_inline _f_hybrid double dQ_x_dx(const double& x, const double& y, const double& z) const { return (gas.visc/gas.prandtl)*d2T_dx2(x,y,z); }
    _f_no_inline _f_hybrid double dQ_y_dy(const double& x, const double& y, const double& z) const { return (gas.visc/gas.prandtl)*d2T_dy2(x,y,z); }
    _f_no_inline _f_hybrid double dQ_z_dz(const double& x, const double& y, const double& z) const { return (gas.visc/gas.prandtl)*d2T_dz2(x,y,z); }
    
    _f_no_inline _f_hybrid void visc_rhs(const double& x, const double& y, const double& z, double (&rhsAr)[5]) const
    {
        rhsAr[0] = 0.0;
        rhsAr[1] = 
            U(x,y,z)*dTau_xx_dx(x,y,z)+dU_dx(x,y,z)*Tau_xx(x,y,z)+
            V(x,y,z)*dTau_xy_dx(x,y,z)+dV_dx(x,y,z)*Tau_xy(x,y,z)+
            W(x,y,z)*dTau_xz_dx(x,y,z)+dW_dx(x,y,z)*Tau_xz(x,y,z)+
            U(x,y,z)*dTau_yx_dy(x,y,z)+dU_dy(x,y,z)*Tau_yx(x,y,z)+
            V(x,y,z)*dTau_yy_dy(x,y,z)+dV_dy(x,y,z)*Tau_yy(x,y,z)+
            W(x,y,z)*dTau_yz_dy(x,y,z)+dW_dy(x,y,z)*Tau_yz(x,y,z)+
            U(x,y,z)*dTau_zx_dz(x,y,z)+dU_dz(x,y,z)*Tau_zx(x,y,z)+
            V(x,y,z)*dTau_zy_dz(x,y,z)+dV_dz(x,y,z)*Tau_zy(x,y,z)+
            W(x,y,z)*dTau_zz_dz(x,y,z)+dW_dz(x,y,z)*Tau_zz(x,y,z)+
            dQ_x_dx(x,y,z)+dQ_y_dy(x,y,z)+dQ_z_dz(x,y,z);
        rhsAr[2] = dTau_xx_dx(x,y,z)+dTau_xy_dy(x,y,z)+dTau_xz_dz(x,y,z);
        rhsAr[3] = dTau_yx_dx(x,y,z)+dTau_yy_dy(x,y,z)+dTau_yz_dz(x,y,z);
        rhsAr[4] = dTau_zx_dx(x,y,z)+dTau_zy_dy(x,y,z)+dTau_zz_dz(x,y,z);
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