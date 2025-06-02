/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(TEUKOLSKYWAVE_HPP_)
#error "This file should only be included through TeukolskyWave.hpp"
#endif

#ifndef TEUKOLSKYWAVE_IMPL_HPP_
#define TEUKOLSKYWAVE_IMPL_HPP_

#include "EppleyPacket.hpp" //for EppleyPacket class

inline TeukolskyWave::TeukolskyWave(EppleyPacket_params_t a_params_eppley_packet, double a_dx)
    : m_dx(a_dx), m_params_eppley_packet(a_params_eppley_packet), 
      m_eppley_packet(EppleyPacketM0(a_params_eppley_packet))
{
}

// Compute the value of the initial vars on the grid
template <class data_t>
void TeukolskyWave::compute(Cell<data_t> current_cell) const
{
    // CHANGE TO VACUUM
    CCZ4Vars::VarsWithGauge<data_t> vars;

    current_cell.load_vars(vars);

    // Coordinates for centre of mass
    Coordinates<data_t> coords(current_cell, m_dx,
                               m_params_eppley_packet.wave_centre);

    // First star positioning
    double x = coords.x;
    double z = coords.z;
    double y = coords.y;
    double r = sqrt(x * x + y * y + z * z);
    double t = m_params_eppley_packet.time_offset;

    // get the metric componennts
    double gxx = m_eppley_packet.get_gxx(x, y, z, r, t);
    double gxy = m_eppley_packet.get_gxy(x, y, z, r, t);
    double gxz = m_eppley_packet.get_gxz(x, y, z, r, t);
    double gyy = m_eppley_packet.get_gyy(x, y, z, r, t);
    double gyz = m_eppley_packet.get_gyz(x, y, z, r, t);
    double gzz = m_eppley_packet.get_gzz(x, y, z, r, t);

    double g[3][3] = {{gxx, gxy, gxz}, {gxy, gyy, gyz}, {gxz, gyz, gzz}};

    // Define the conformal factor
    double det_g = g[0][0] * (g[1][1] * g[2][2] - g[1][2] * g[2][1]) -
                 g[0][1] * (g[1][0] * g[2][2] - g[1][2] * g[2][0]) +
                 g[0][2] * (g[1][0] * g[2][1] - g[1][1] * g[2][0]);
    double chi = pow(det_g, -1. / 3.);

    // Initialise conformal metric
    double h[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
    FOR2(i, j) h[i][j] = g[i][j] / chi;

    // Define initial conformal factor
    vars.chi += chi;

    // Define initial lapse
    vars.lapse += 1.0;

    // Define initial trace of K and A_ij
    FOR2(i, j) vars.h[i][j] = h[i][j];
    vars.K = 0.;
    FOR2(i, j) vars.A[i][j] = 0.;

    current_cell.store_vars(vars);
}

#endif /* TEUKOLSKYWAVE_IMPL_HPP_ */
