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

EppleyPacket::EppleyPacket(EppleyPacket_params_t m_params)
    : amplitude(m_params.amplitude),
      sigma(m_params.sigma)
{
    // Constructor implementation
}

// F function and its derivatives where x = r \pm t
double EppleyPacket::get_F(double x) const
{
    return amplitude * x * exp(-x * x / (sigma * sigma));
}

double EppleyPacket::get_Fd1(double x) const
{
    return amplitude * (1. - 2. * x * x / (sigma * sigma)) * exp(-x * x / (sigma * sigma));
}

double EppleyPacket::get_Fd2(double x) const
{
    return amplitude * x * (4. * x * x / (sigma * sigma) - 6.) * exp(-x * x / (sigma * sigma)) / (sigma * sigma);
}

double EppleyPacket::get_Fd3(double x) const
{
    return amplitude * (-8. * pow(x, 4) / pow(sigma, 4)  +24. * x * x / (sigma * sigma) - 6.) * exp(-x * x / (sigma * sigma)) / (sigma * sigma);
}

double EppleyPacket::get_Fd4(double x) const
{
    return amplitude * x *(16. * pow(x, 4) / pow(sigma, 4) -80. * x*x/(sigma*sigma) + 60.) * exp(-x * x / (sigma * sigma)) / pow(sigma, 4);
}

// Auxiliary functions. In the end we want the superposition, so these are also implemented as get_X_tot

double EppleyPacket::get_A(double r, double t, int sign) const
{
    double x = t - sign * r;
    return 3 * get_Fd2(x) / pow(r, 3) + sign * 9. * get_Fd1(x)  / pow(r, 4) + 3. * get_F(x) / pow(r, 5);
}

double EppleyPacket::get_A_tot(double r, double t) const
{
    return get_A(r, t, 1) - get_A(r, t, -1);
}

double EppleyPacket::get_B(double r, double t, int sign) const
{
    double x = t - sign * r;
    return -1. *sign * get_Fd3(x) / (r*r)
           -3. * get_Fd2(x) / pow(r, 3)
           -6. * sign * get_Fd1(x) / pow(r, 4)
           -6. * get_F(x) / pow(r, 5);
}

double EppleyPacket::get_B_tot(double r, double t) const
{
    return get_B(r, t, 1) - get_B(r, t, -1);
}

double EppleyPacket::get_C(double r, double t, int sign) const
{
    double x = t - sign * r;
    return 0.25 * get_Fd4(x) / r
           + 0.5 * sign * get_Fd3(x) / (r*r)
           + 2.25 * get_Fd2(x) / pow(r, 3)
           + 5.25 * sign * get_Fd1(x) / pow(r, 4)
           + 5.25 * get_F(x) / pow(r, 5);
}

double EppleyPacket::get_C_tot(double r, double t) const
{
    return get_C(r, t, 1) - get_C(r, t, -1);
}

// m = 0 EppleyPacket

double EppleyPacketM0::get_gxx(double x, double y, double z, double r, double t) const
{
    return 1. + (-1. + 3.*y*y/(r*r) + 3.*x*x*z*z/pow(r,4)) * get_A_tot(r, t)
              - 6. * z*z * x*x * get_B_tot(r, t) / pow(r, 4)
              + 3. * (- y*y / (r*r) + x*x*z*z/pow(r,4)) * get_C_tot(r, t);
}

double EppleyPacketM0::get_gxy(double x, double y, double z, double r, double t) const
{
    return 3. * x * y * (-1.*get_A_tot(r, t)* (x*x + y*y) - 2*z*z*get_B_tot(r, t) + (r*r+z*z)*get_C_tot(r, t)) / pow(r, 4);
}

double EppleyPacketM0::get_gxz(double x, double y, double z, double r, double t) const
{
    return 3. * x * z * (z*z*get_A_tot(r, t) + (x*x + y*y - z*z)*get_B_tot(r, t) - (x*x+y*y)*get_C_tot(r, t)) / pow(r, 4);
}

double EppleyPacketM0::get_gyy(double x, double y, double z, double r, double t) const
{
    return 1. + (-1. + 3.*x*x/(r*r) + 3.*y*y*z*z/pow(r,4)) * get_A_tot(r, t)
              - 6. * z*z * y*y * get_B_tot(r, t) / pow(r, 4)
              + 3. * (- x*x / (r*r) + y*y*z*z/pow(r,4)) * get_C_tot(r, t);
}

double EppleyPacketM0::get_gyz(double x, double y, double z, double r, double t) const
{
    return 3. * y * z * (z*z*get_A_tot(r, t) + (x*x + y*y - z*z)*get_B_tot(r, t) - (x*x+y*y)*get_C_tot(r, t)) / pow(r, 4);
}

double EppleyPacketM0::get_gzz(double x, double y, double z, double r, double t) const
{
    return 1. + (-1. + 3.*pow(z, 4) / pow(r, 4)) * get_A_tot(r, t)
              + 6. * z*z * (x*x + y*y) * get_B_tot(r, t) / pow(r, 4)
                + 3. * (x*x + y*y)*(x*x + y*y) * get_C_tot(r, t) / pow(r, 4);
}

#endif /* EPPLEYPACKET_IMPL_HPP_ */
