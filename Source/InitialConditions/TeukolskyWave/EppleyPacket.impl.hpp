/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(EPPLEYPACKET_HPP_)
#error "This file should only be included through EppleyPacket.hpp"
#endif

#ifndef EPPLEYPACKET_IMPL_HPP_
#define EPPLEYPACKET_IMPL_HPP_

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

// F function and its derivatives where x = r \pm t
EppleyPacketDerivs EppleyPacket::get_F_derivs(double x) const
{
    double A = this->m_amplitude, sigma = this->m_sigma, r0 = this->m_radial_offset;

    // --- paste generated temporaries here ---
    double exp_plus = exp(-pow(x + r0, 2) / pow(sigma, 2));
    double exp_minus = exp(-pow(x - r0, 2) / pow(sigma, 2));
    double sigma2 = sigma * sigma;
    double sigma4 = sigma2 * sigma2;
    double sigma6 = sigma4 * sigma2;
    double sigma8 = sigma4 * sigma4;

    // --- paste generated F0..F4 assignments here ---
    double F0, F1, F2, F3, F4;

    F0 = (exp_plus + exp_minus)*x;
    
    F1 = exp_minus + exp_plus + (2*exp_minus*(r0 - x)*x)/sigma2 - (2*exp_plus*x*(r0 \
    + x))/sigma2;

    F2 = (4*(exp_minus*(r0 - x) - exp_plus*(r0 + x)))/sigma2 - \
    (2*x*(exp_minus*sigma2 + exp_plus*sigma2 - 2*exp_minus*pow(r0 - \
    x,2) - 2*exp_plus*pow(r0 + x,2)))/sigma4;

    F3 = (-2*(3*sigma2*(exp_minus*sigma2 + exp_plus*sigma2 \
    - 2*exp_minus*pow(r0 - x,2) - 2*exp_plus*pow(r0 + x,2)) + \
    2*x*(3*exp_minus*sigma2*(r0 - x) - 2*exp_minus*pow(r0 - x,3) - \
    3*exp_plus*sigma2*(r0 + x) + 2*exp_plus*pow(r0 + \
    x,3))))/sigma6;

    F4 = (4*(-4*sigma2*(3*exp_minus*sigma2*(r0 - x) - \
    2*exp_minus*pow(r0 - x,3) - 3*exp_plus*sigma2*(r0 + x) + \
    2*exp_plus*pow(r0 + x,3)) + x*(3*exp_minus*sigma4 + \
    3*exp_plus*sigma4 - 12*exp_minus*sigma2*pow(r0 - x,2) + \
    4*exp_minus*pow(r0 - x,4) - 12*exp_plus*sigma2*pow(r0 + x,2) + \
    4*exp_plus*pow(r0 + x,4))))/sigma8;

    // F0 = pow(x, 9)*(exp_plus + exp_minus);
    
    // F1 = pow(x,8)*(9*(exp_plus + exp_minus) - (2*x*(exp_minus*(-r0 + x) + \
F4 = (4*(-4*sigma2*(3*expM*sigma2*(r0 - x) - \
2*expM*pow(r0 - x,3) - 3*expP*sigma2*(r0 + x) + \
2*expP*pow(r0 + x,3)) + x*(3*expM*sigma4 + \
3*expP*sigma4 - 12*expM*sigma2*pow(r0 - x,2) + \
4*expM*pow(r0 - x,4) - 12*expP*sigma2*pow(r0 + x,2) + \
4*expP*pow(r0 + x,4))))/sigma8;

    // F0 = pow(x, 9)*(exp_plus + exp_minus);
    
    // F1 = pow(x,8)*(9*(exp_plus + exp_minus) - (2*x*(exp_minus*(-r0 + x) + \
    //      exp_plus*(r0 + x)))/sigma2);

    // F2 = (2*pow(x,7)*(exp_minus*(36*sigma4 + sigma2*(18*r0 - \
    //       19*x)*x + 2*pow(r0 - x,2)*x*x) + exp_plus*(36*sigma4 + \
    //       2*x*x*pow(r0 + x,2) - sigma2*x*(18*r0 + \
    //       19*x))))/sigma4;

    // F3 = (2*pow(x,6)*(exp_minus*(252*sigma6 + \
    //       27*sigma4*(8*r0 - 9*x)*x + 4*pow(r0 - x,3)*pow(x,3) + \
    //       6*sigma2*x*x*(9*r0*r0 - 19*r0*x + 10*x*x)) \
    //       + exp_plus*(252*sigma6 - 4*pow(x,3)*pow(r0 + x,3) - \
    //       27*sigma4*x*(8*r0 + 9*x) + \
    //       6*sigma2*x*x*(9*r0*r0 + 19*r0*x + \
    //       10*x*x))))/sigma6;

    // F4 = (4*pow(x,5)*(exp_minus*(756*sigma8 + \
    //     72*sigma6*(14*r0 - 17*x)*x + 12*sigma2*(6*r0 - \
    //     7*x)*pow(r0 - x,2)*pow(x,3) + 4*pow(r0 - x,4)*pow(x,4) + \
    //     3*sigma4*x*x*(144*r0*r0 - 324*r0*x + \
    //     181*x*x)) + exp_plus*(756*sigma8 + 4*pow(x,4)*pow(r0 + \
    //     x,4) - 12*sigma2*pow(x,3)*pow(r0 + x,2)*(6*r0 + 7*x) - \
    //     72*sigma6*x*(14*r0 + 17*x) + \
    //     3*sigma4*x*x*(144*r0*r0 + 324*r0*x + \
    //     181*x*x))))/sigma8;

    return EppleyPacketDerivs{A*F0/2., A*F1/2., A*F2/2., A*F3/2., A*F4/2.};
}
// double EppleyPacket::get_F(double x) const
// {
//     return amplitude * x * exp(-x * x / (sigma * sigma));
// }

// double EppleyPacket::get_Fd1(double x) const
// {
//     return amplitude * (1. - 2. * x * x / (sigma * sigma)) * exp(-x * x / (sigma * sigma));
// }

// double EppleyPacket::get_Fd2(double x) const
// {
//     return amplitude * x * (4. * x * x / (sigma * sigma) - 6.) * exp(-x * x / (sigma * sigma)) / (sigma * sigma);
// }

// double EppleyPacket::get_Fd3(double x) const
// {
//     return amplitude * (-8. * pow(x, 4) / pow(sigma, 4)  +24. * x * x / (sigma * sigma) - 6.) * exp(-x * x / (sigma * sigma)) / (sigma * sigma);
// }

// double EppleyPacket::get_Fd4(double x) const
// {
//     return amplitude * x *(16. * pow(x, 4) / pow(sigma, 4) -80. * x*x/(sigma*sigma) + 60.) * exp(-x * x / (sigma * sigma)) / pow(sigma, 4);
// }

// Auxiliary functions. In the end we want the superposition, so these are also implemented as get_X_tot
EvenEppleyPacketCoefficients EvenEppleyPacket::get_ABC(double r) const
{
    double x_out = - r;
    double x_in = r;
    EppleyPacketDerivs F_derivs_out = get_F_derivs(x_out);
    EppleyPacketDerivs F_derivs_in = get_F_derivs(x_in);

    // Compute out coefficients
    double A_out = 3 * F_derivs_out.F2 / pow(r, 3) 
                   + 9. * F_derivs_out.F1 / pow(r, 4) 
                   + 9. * F_derivs_out.F0 / pow(r, 5);
    double B_out = -1. * F_derivs_out.F3 / (r*r)
                   -3. * F_derivs_out.F2 / pow(r, 3)
                   -6. * F_derivs_out.F1 / pow(r, 4)
                   -6. * F_derivs_out.F0 / pow(r, 5);
    double C_out = 0.25 * F_derivs_out.F4 / r
                   + 0.5 * F_derivs_out.F3 / (r*r)
                   + 2.25 * F_derivs_out.F2 / pow(r, 3)
                   + 5.25 * F_derivs_out.F1 / pow(r, 4)
                   + 5.25 * F_derivs_out.F0 / pow(r, 5);

    // Compute in coefficients
    double A_in = 3 * F_derivs_in.F2 / pow(r, 3) 
                  - 9. * F_derivs_in.F1 / pow(r, 4) 
                  + 9. * F_derivs_in.F0 / pow(r, 5); 
    double B_in = 1. * F_derivs_in.F3 / (r*r)
                  -3. * F_derivs_in.F2 / pow(r, 3)
                  +6. * F_derivs_in.F1 / pow(r, 4)
                  -6. * F_derivs_in.F0 / pow(r, 5);
    double C_in = 0.25 * F_derivs_in.F4 / r
                  - 0.5 * F_derivs_in.F3 / (r*r)
                  + 2.25 * F_derivs_in.F2 / pow(r, 3)
                  - 5.25 * F_derivs_in.F1 / pow(r, 4)
                  + 5.25 * F_derivs_in.F0 / pow(r, 5);

    return EvenEppleyPacketCoefficients{A_out - A_in, B_out - B_in, C_out - C_in};
}

// double EvenEppleyPacket::get_A(double r, double r0, int sign) const
// {
//     // t - sign * r = - sign * r at t=0
//     double x = - sign * r;
//     EppleyPacketDerivs F_derivs = get_F_derivs(x, r0);
//     return 3 * F_derivs.F2 / pow(r, 3) + sign * 9. * F_derivs.F1 / pow(r, 4) + 3. * F_derivs.F0 / pow(r, 5);
// }

// double EvenEppleyPacket::get_A_tot(double r, double t) const
// {
//     return get_A(r, t, 1) - get_A(r, t, -1);
// }

// double EvenEppleyPacket::get_B(double r, double t, int sign) const
// {
//     double x = t - sign * r;
//     return -1. *sign * get_Fd3(x) / (r*r)
//            -3. * get_Fd2(x) / pow(r, 3)
//            -6. * sign * get_Fd1(x) / pow(r, 4)
//            -6. * get_F(x) / pow(r, 5);
// }

// double EvenEppleyPacket::get_B_tot(double r, double t) const
// {
//     return get_B(r, t, 1) - get_B(r, t, -1);
// }

// double EvenEppleyPacket::get_C(double r, double t, int sign) const
// {
//     double x = t - sign * r;
//     return 0.25 * get_Fd4(x) / r
//            + 0.5 * sign * get_Fd3(x) / (r*r)
//            + 2.25 * get_Fd2(x) / pow(r, 3)
//            + 5.25 * sign * get_Fd1(x) / pow(r, 4)
//            + 5.25 * get_F(x) / pow(r, 5);
// }

// double EvenEppleyPacket::get_C_tot(double r, double t) const
// {
//     return get_C(r, t, 1) - get_C(r, t, -1);
// }

OddEppleyPacketCoefficients OddEppleyPacket::get_KL(double r) const
{
    double x_out = - r;
    double x_in = r;
    EppleyPacketDerivs F_derivs_out = get_F_derivs(x_out);
    EppleyPacketDerivs F_derivs_in = get_F_derivs(x_in);

    // Compute out coefficients
    double K_out = F_derivs_out.F2 / pow(r, 2) + 3. * F_derivs_out.F1 / pow(r, 3) + 3. * F_derivs_out.F0 / pow(r, 4);
    double L_out = F_derivs_out.F3 / r + 2. * F_derivs_out.F2 / pow(r, 2) + 3. * F_derivs_out.F1 / pow(r, 3) + 3. * F_derivs_out.F0 / pow(r, 4);

    // Compute in coefficients
    double K_in = F_derivs_in.F2 / pow(r, 2) - 3. * F_derivs_in.F1 / pow(r, 3) + 3. * F_derivs_in.F0 / pow(r, 4);
    double L_in = -1.*F_derivs_in.F3 / r + 2. * F_derivs_in.F2 / pow(r, 2) - 3. * F_derivs_in.F1 / pow(r, 3) + 3. * F_derivs_in.F0 / pow(r, 4);

    return OddEppleyPacketCoefficients{K_out - K_in, L_out - L_in};
}

// double OddEppleyPacket::get_K(double r, double t, int sign) const
// {
//     double x = t - sign * r;
//     return get_Fd2(x) / pow(r, 2) + sign * 3. * get_Fd1(x) / pow(r, 3) + 3. * get_F(x) / pow(r, 4);
// }

// double OddEppleyPacket::get_K_tot(double r, double t) const
// {
//     return get_K(r, t, 1) - get_K(r, t, -1);
// }

// double OddEppleyPacket::get_L(double r, double t, int sign) const
// {
//     double x = t - sign * r;
//     return sign * get_Fd3(x) / r + 2. * get_Fd2(x) / pow(r, 2) + 3. * sign * get_Fd1(x) / pow(r, 3) + 3. * get_F(x) / pow(r, 4);
// }

// double OddEppleyPacket::get_L_tot(double r, double t) const
// {
//     return get_L(r, t, 1) - get_L(r, t, -1);
// }

// ------------- m = 0 EppleyPacket -----------------

EppleyPacketMetricComponents EppleyPacketM0::get_metric_components(double x, double y, double z) const
{
    double r = sqrt(x * x + y * y + z * z);
    // regularize at the origin
    r = r + 1. * exp(-r*r);
    EvenEppleyPacketCoefficients coeffs = get_ABC(r);
    double A_tot = coeffs.A;
    double B_tot = coeffs.B;
    double C_tot = coeffs.C;
    EppleyPacketMetricComponents components;
    components.gxx = 1. + (-1. + 3.*y*y/(r*r) + 3.*x*x*z*z/pow(r,4)) * A_tot
                    - 6. * z*z * x*x * B_tot / pow(r, 4)
                    + 3. * (- y*y / (r*r) + x*x*z*z/pow(r,4)) * C_tot;
    components.gxy = 3. * x * y * (-1.*A_tot* (x*x + y*y) - 2*z*z*B_tot + (r*r+z*z)*C_tot) / pow(r, 4);
    components.gxz = 3. * x * z * (z*z*A_tot + (x*x + y*y - z*z)*B_tot - (x*x+y*y)*C_tot) / pow(r,4);
    components.gyy = 1. + (-1. + 3.*x*x/(r*r) + 3.*y*y*z*z/pow(r,4)) * A_tot
                    - 6. * z*z * y*y *B_tot / pow(r, 4)
                    + 3. * (- x*x / (r*r) + y*y*z*z/pow(r,4)) * C_tot;
    components.gyz = 3. * y * z * (z*z*A_tot + (x*x + y*y - z*z)*B_tot - (x*x+y*y)*C_tot) / pow(r, 4);
    components.gzz = 1. + (-1. + 3.*pow(z, 4) / pow(r, 4)) * A_tot
                    + 6. * z*z * (x*x + y*y) * B_tot / pow(r, 4)
                    + 3. * (x*x + y*y)*(x*x + y*y) * C_tot / pow(r, 4);
    return components;
}

// double EppleyPacketM0::get_gxx(double x, double y, double z, double r, double t) const
// {
//     return 1. + (-1. + 3.*y*y/(r*r) + 3.*x*x*z*z/pow(r,4)) * get_A_tot(r, t)
//               - 6. * z*z * x*x * get_B_tot(r, t) / pow(r, 4)
//               + 3. * (- y*y / (r*r) + x*x*z*z/pow(r,4)) * get_C_tot(r, t);
// }

// double EppleyPacketM0::get_gxy(double x, double y, double z, double r, double t) const
// {
//     return 3. * x * y * (-1.*get_A_tot(r, t)* (x*x + y*y) - 2*z*z*get_B_tot(r, t) + (r*r+z*z)*get_C_tot(r, t)) / pow(r, 4);
// }

// double EppleyPacketM0::get_gxz(double x, double y, double z, double r, double t) const
// {
//     return 3. * x * z * (z*z*get_A_tot(r, t) + (x*x + y*y - z*z)*get_B_tot(r, t) - (x*x+y*y)*get_C_tot(r, t)) / pow(r, 4);
// }

// double EppleyPacketM0::get_gyy(double x, double y, double z, double r, double t) const
// {
//     return 1. + (-1. + 3.*x*x/(r*r) + 3.*y*y*z*z/pow(r,4)) * get_A_tot(r, t)
//               - 6. * z*z * y*y * get_B_tot(r, t) / pow(r, 4)
//               + 3. * (- x*x / (r*r) + y*y*z*z/pow(r,4)) * get_C_tot(r, t);
// }

// double EppleyPacketM0::get_gyz(double x, double y, double z, double r, double t) const
// {
//     return 3. * y * z * (z*z*get_A_tot(r, t) + (x*x + y*y - z*z)*get_B_tot(r, t) - (x*x+y*y)*get_C_tot(r, t)) / pow(r, 4);
// }

// double EppleyPacketM0::get_gzz(double x, double y, double z, double r, double t) const
// {
//     return 1. + (-1. + 3.*pow(z, 4) / pow(r, 4)) * get_A_tot(r, t)
//               + 6. * z*z * (x*x + y*y) * get_B_tot(r, t) / pow(r, 4)
//                 + 3. * (x*x + y*y)*(x*x + y*y) * get_C_tot(r, t) / pow(r, 4);
// }

// ------------- m = 2 EppleyPacket -----------------

EppleyPacketMetricComponents EppleyPacketM2::get_metric_components(double x, double y, double z) const
{
    double r = sqrt(x * x + y * y + z * z);
    EvenEppleyPacketCoefficients coeffs = get_ABC(r);
    double A_tot = coeffs.A;
    double B_tot = coeffs.B;
    double C_tot = coeffs.C;
    EppleyPacketMetricComponents components;
    components.gxx = 1. + ((x*x - z*z) / (r*r) - x*x * (z*z + 2*y*y)/pow(r,4)) * A_tot
                    + 2*x*x*(z*z + 2*y*y) * B_tot / pow(r, 4)
                    + ((y*y+2*z*z)/(r*r) - x*x * (z*z + 2*y*y)/pow(r,4)) * C_tot;
    components.gxy = x*y*(x*x - y*y) * (A_tot - 2*B_tot + C_tot) / pow(r, 4);
    components.gxz = x*z * ( (2*x*x + z*z)*A_tot + (z*z + 3*y*y - x*x)*B_tot - (x*x + 2*z*z + 3*y*y) * C_tot) / pow(r,4);
    components.gyy = 1. + ((z*z - y*y) / (r*r) + y*y * (z*z + 2*x*x)/pow(r,4)) * A_tot
                    - 2*y*y*(z*z + 2*x*x) * B_tot / pow(r, 4)
                    + (-(x*x+2*z*z)/(r*r) + y*y * (z*z + 2*x*x)/pow(r,4)) * C_tot;
    components.gyz = y*z * ( -(2*y*y + z*z)*A_tot - (z*z + 3*x*x - y*y)*B_tot + (y*y + 2*z*z + 3*x*x) * C_tot) / pow(r,4);
    components.gzz = 1. + ((pow(y,4) - pow(x,4)) * A_tot
                    - 2. * z*z * (x*x - y*y) * B_tot
                    + (x*x - y*y)*(r*r + z*z) * C_tot) / pow(r, 4);
    return components;
}

// double EppleyPacketM2::get_gxx(double x, double y, double z, double r, double t) const
// {
//     return 1. + ((x*x - z*z) / (r*r) - x*x * (z*z + 2*y*y)/pow(r,4)) * get_A_tot(r, t)
//               + 2*x*x*(z*z + 2*y*y) * get_B_tot(r, t) / pow(r, 4)
//               + ((y*y+2*z*z)/(r*r) - x*x * (z*z + 2*y*y)/pow(r,4)) * get_C_tot(r, t);
// }

// double EppleyPacketM2::get_gxy(double x, double y, double z, double r, double t) const
// {
//     return x*y*(x*x - y*y) * (get_A_tot(r, t) - 2*get_B_tot(r, t) + get_C_tot(r, t)) / pow(r, 4);
// }           

// double EppleyPacketM2::get_gxz(double x, double y, double z, double r, double t) const
// {
//     return x*z * ( (2*x*x + z*z)*get_A_tot(r, t) + (z*z + 3*y*y - x*x)*get_B_tot(r, t) - (x*x + 2*z*z + 3*y*y) * get_C_tot(r, t)) / pow(r,4);
// }

// double EppleyPacketM2::get_gyy(double x, double y, double z, double r, double t) const
// {
//     return 1. + ((z*z - y*y) / (r*r) + y*y * (z*z + 2*x*x)/pow(r,4)) * get_A_tot(r, t)
//               - 2*y*y*(z*z + 2*x*x) * get_B_tot(r, t) / pow(r, 4)
//               + (-(x*x+2*z*z)/(r*r) + y*y * (z*z + 2*x*x)/pow(r,4)) * get_C_tot(r, t);
// }

// double EppleyPacketM2::get_gyz(double x, double y, double z, double r, double t) const
// {
//     return y*z * ( -(2*y*y + z*z)*get_A_tot(r, t) - (z*z + 3*x*x - y*y)*get_B_tot(r, t) + (y*y + 2*z*z + 3*x*x) * get_C_tot(r, t)) / pow(r,4);
// }

// double EppleyPacketM2::get_gzz(double x, double y, double z, double r, double t) const
// {
//     return 1. + ((pow(y,4) - pow(x,4)) * get_A_tot(r, t)
//               - 2. * z*z * (x*x - y*y) * get_B_tot(r, t)
//               + (x*x - y*y)*(r*r + z*z) * get_C_tot(r, t)) / pow(r, 4);
// }

// -------------- m = 2 Odd parity EppleyPacket -----------------

EppleyPacketMetricComponents OddEppleyPacketM2::get_metric_components(double x, double y, double z) const
{
    double r = sqrt(x * x + y * y + z * z);
    OddEppleyPacketCoefficients coeffs = get_KL(r);
    double K_tot = coeffs.K;
    double L_tot = coeffs.L;
    EppleyPacketMetricComponents components;
    components.gxx = 1. + 2.*z * pow(r,-3) * (-4 * K_tot * x*x *(x*x - 3* y*y)*(x*x + y*y) - L_tot*r*r*(3* x*x * y*y + pow(y,4)) 
                    + L_tot* x*x *(2* y*y *(x*x + y*y) + (y*y-x*x)* z*z))*pow(x*x + y*y,-2);
    components.gxy = 16*K_tot * x*y*z * pow(r,-3) * (-x*x + y*y)*pow(x*x + y*y,-1);
    components.gxz = 2*x*pow(r,-3)*(2*L_tot*r*r*y*y - L_tot*y*y*(x*x + y*y) + 2*K_tot*(pow(x,4) - pow(y,4)) + ((-2*K_tot + L_tot)*x*x + (6*K_tot - L_tot)*y*y)*pow(z,2))*
                    pow(x*x + y*y,-1);
    components.gyy = 1. + 2*z*pow(r,-3)*(L_tot*x*x*(-2*y*y*(x*x + y*y) + r*r*(x*x + 3*y*y)) + 4*K_tot*y*y*(-3*pow(x,4) - 2*x*x*y*y + pow(y,4)) 
                    + L_tot*(-x*x + y*y) * y*y * z*z)*pow(x*x + y*y,-2);
    components.gyz = 2*y*pow(r,-3)*(-2*L_tot* r*r * x*x + 2*K_tot*pow(x,4) + L_tot*pow(x,4) + L_tot * x*x * y*y - 2*K_tot*pow(y,4) 
                    + ((-6*K_tot + L_tot)*x*x + (2*K_tot - L_tot)*y*y)*pow(z,2))*pow(x*x + y*y,-1);
    components.gzz = 1. + 2*(4*K_tot - L_tot)*(x*x - y*y)*z*pow(r,-3);
    return components;
}

// double OddEppleyPacketM2::get_gxx(double x, double y, double z, double r, double t) const
// {
//     return 1. + 2.*z * pow(r,-3) * (-4 * get_K_tot(r, t) * x*x *(x*x - 3* y*y)*(x*x + y*y) - get_L_tot(r, t)*r*r*(3* x*x * y*y + pow(y,4)) 
//     + get_L_tot(r, t)* x*x *(2* y*y *(x*x + y*y) + (y*y-x*x)* z*z))*pow(x*x + y*y,-2);
// }

// double OddEppleyPacketM2::get_gxy(double x, double y, double z, double r, double t) const
// {
//     return 16*get_K_tot(r, t) * x*y*z * pow(r,-3) * (-x*x + y*y)*pow(x*x + y*y,-1);
// }

// double OddEppleyPacketM2::get_gxz(double x, double y, double z, double r, double t) const
// {
//     return 2*x*pow(r,-3)*(2*get_L_tot(r, t)*r*r*y*y - get_L_tot(r, t)*y*y*(x*x + y*y) + 2*get_K_tot(r, t)*(pow(x,4) - pow(y,4)) + ((-2*get_K_tot(r, t) + get_L_tot(r, t))*x*x + (6*get_K_tot(r, t) - get_L_tot(r, t))*y*y)*pow(z,2))*
//    pow(x*x + y*y,-1);
// }

// double OddEppleyPacketM2::get_gyy(double x, double y, double z, double r, double t) const
// {
//     return 1 + 2*z*pow(r,-3)*(get_L_tot(r, t)*x*x*(-2*y*y*(x*x + y*y) + r*r*(x*x + 3*y*y)) + 4*get_K_tot(r, t)*y*y*(-3*pow(x,4) - 2*x*x*y*y + pow(y,4)) + 
//       get_L_tot(r, t)*(-x*x + y*y) * y*y * z*z)*pow(x*x + y*y,-2);
// }

// double OddEppleyPacketM2::get_gyz(double x, double y, double z, double r, double t) const
// {
//     return 2*y*pow(r,-3)*(-2*get_L_tot(r, t)* r*r * x*x + 2*get_K_tot(r, t)*pow(x,4) + get_L_tot(r, t)*pow(x,4) + get_L_tot(r, t) * x*x * y*y - 2*get_K_tot(r, t)*pow(y,4) + ((-6*get_K_tot(r, t) + get_L_tot(r, t))*x*x + (2*get_K_tot(r, t) - get_L_tot(r, t))*y*y)*pow(z,2))*
//    pow(x*x + y*y,-1);
// }

// double OddEppleyPacketM2::get_gzz(double x, double y, double z, double r, double t) const
// {
//     return 1 + 2*(4*get_K_tot(r, t) - get_L_tot(r, t))*(x*x - y*y)*z*pow(r,-3);
// }

#endif /* EPPLEYPACKET_IMPL_HPP_ */
