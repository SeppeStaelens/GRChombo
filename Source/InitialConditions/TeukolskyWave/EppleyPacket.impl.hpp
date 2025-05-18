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

// Auxiliary functions. Note that we already imply the superposition here,
// meaning that we implement e.g. A_out  - A_in

double EppleyPacket::get_A(double r, double t) const
{
    return 6. * (get_Fd1(t - r) - get_Fd1(t + r)) / pow(r, 4);
}

double EppleyPacket::get_B(double r, double t) const
{
    return (get_Fd3(t+r) - get_Fd3(t-r)) / (r*r) + 
           6. * (get_Fd1(t+r) - get_Fd1(t-r)) / pow(r, 4);
}

double EppleyPacket::get_C(double r, double t) const
{
    return 0.5 * (get_Fd3(t-r) - get_Fd3(t+r)) / (r*r) +
           5.25 * (get_Fd1(t-r) - get_Fd1(t+r)) / pow(r, 4);
}

// m = 0 EppleyPacket

double EppleyPacketM0::get_gxx(double x, double y, double z, double r) const
{
    return 1. + (-1. + 3.*y*y/(r*r) + 3.*x*x*z*z/pow(r,4)) * get_A(r)
              - 6. * z*z * x*x * get_B(r) / pow(r, 4)
              + 3. * (- y*y / (r*r) + x*x*z*z/pow(r,4)) * get_C(r);
}

double EppleyPacketM0::get_gxy(double x, double y, double z, double r) const
{
    return 3. * x * y * (-1.*get_A(r) (x*x + y*y) - 2*z*z*get_B(r) + (r*r+z*z)*get_C(r)) / pow(r, 4);
}

double EppleyPacketM0::get_gxz(double x, double y, double z, double r) const
{
    return 3. * x * z * (z*z*get_A(r) + (x*x + y*y - z*z)*get_B(r) - (x*x+y*y)*get_C(r)) / pow(r, 4);
}

double EppleyPacketM0::get_gyy(double x, double y, double z, double r) const
{
    return 1. + (-1. + 3.*x*x/(r*r) + 3.*y*y*z*z/pow(r,4)) * get_A(r)
              - 6. * z*z * y*y * get_B(r) / pow(r, 4)
              + 3. * (- x*x / (r*r) + y*y*z*z/pow(r,4)) * get_C(r);
}

double EppleyPacketM0::get_gyz(double x, double y, double z, double r) const
{
    return 3. * y * z * (z*z*get_A(r) + (x*x + y*y - z*z)*get_B(r) - (x*x+y*y)*get_C(r)) / pow(r, 4);
}

double EppleyPacketM0::get_gzz(double x, double y, double z, double r) const
{
    return 1. + (-1. + 3.*pow(z, 4) / pow(r, 4)) * get_A(r)
              + 6. * z*z * (x*x + y*y) * get_B(r) / pow(r, 4)
                + 3. * (x*x + y*y)*(x*x + y*y) * get_C(r) / pow(r, 4);
}

#endif /* EPPLEYPACKET_IMPL_HPP_ */
