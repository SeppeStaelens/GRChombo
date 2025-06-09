/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(TSBOSONSTARTEUKOLSKYWAVE_HPP_)
#error "This file should only be included through TSBosonStarTeukolskyWave.hpp"
#endif

#ifndef TSBOSONSTARTEUKOLSKYWAVE_IMPL_HPP_
#define TSBOSONSTARTEUKOLSKYWAVE_IMPL_HPP_

#include "ThinShellSolution.hpp" //for BosonStarSolution class
#include "EppleyPacket.hpp" //for EppleyPacket class

template <class packet_t>
inline TSBosonStarTeukolskyWave<packet_t>::TSBosonStarTeukolskyWave(BosonStar_params_t a_params_BosonStar,
                                        Potential::params_t a_params_potential,
                                        EppleyPacket_params_t a_params_eppley_packet,
                                        double a_dx, double a_L)
    : m_dx(a_dx), m_params_BosonStar(a_params_BosonStar),
      m_params_potential(a_params_potential), m_params_eppley_packet(a_params_eppley_packet),
      m_eppley_packet(packet_t(a_params_eppley_packet))
{
    read_1d_solution(4. * a_L);
}

template <class packet_t>
void TSBosonStarTeukolskyWave<packet_t>::read_1d_solution(const double max_r)
{
    try
    {
        // set initial parameters and then read the solution
        m_1d_sol.set_initialcondition_params(m_params_BosonStar,
                                             m_params_potential, max_r);

        pout() << "Reading boson star functions from file" << endl;
        m_1d_sol.initialise_from_file();
        pout() << "Completed for star 1" << endl;
    }
    catch (std::exception &exception)
    {
        pout() << exception.what() << "\n";
    }
}

// Compute the value of the initial vars on the grid
template <class packet_t>
template <class data_t>
void TSBosonStarTeukolskyWave<packet_t>::compute(Cell<data_t> current_cell) const
{
    MatterCCZ4<ComplexScalarField<>>::Vars<data_t> vars;
    // Load variables (should be set to zero if this is a single BS)

    current_cell.load_vars(vars);
    // VarsTools::assign(vars, 0.); // Set only the non-zero components below

    // Coordinates for centre of mass
    Coordinates<data_t> coords(current_cell, m_dx,
                               m_params_BosonStar.star_centre);

    // Boson star positioning
    double x = coords.x;
    double z = coords.z;
    double y = coords.y;
    double r = sqrt(x * x + y * y + z * z);
    double areal_r = m_1d_sol.r_from_R_Spline(r);

    // auxiliary variables
    double Phi_ = m_1d_sol.PhiSpline(areal_r);
    double f_ = m_1d_sol.fSpline(areal_r);
    double A_ = m_1d_sol.ASpline(areal_r);
    double w_ = m_1d_sol.get_BSfrequency();

    // First star physical variables
    double lapse_bs = exp(Phi_);

    // Add on to evolution equations
    vars.phi_Re = A_;
    vars.phi_Im = 0.;
    vars.Pi_Re = 0.;
    vars.Pi_Im += -(1. / lapse_) * (w_ * A_);

    // Initialise conformal metric
    double h[3][3] = {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};

    double chi_bs = f_ * f_;

    double g_bs[3][3] = {{1. / chi_bs, 0., 0.}, {0., 1. / chi_bs, 0.}, {0., 0., 1. / chi_bs}};

    // ----- Add the Teukolsky wave packet -----
    double t = m_params_eppley_packet.time_offset;

    // get the metric componennts
    double gxx = m_eppley_packet.get_gxx(x, y, z, r, t);
    double gxy = m_eppley_packet.get_gxy(x, y, z, r, t);
    double gxz = m_eppley_packet.get_gxz(x, y, z, r, t);
    double gyy = m_eppley_packet.get_gyy(x, y, z, r, t);
    double gyz = m_eppley_packet.get_gyz(x, y, z, r, t);
    double gzz = m_eppley_packet.get_gzz(x, y, z, r, t);

    double g[3][3] = {{gxx, gxy, gxz}, {gxy, gyy, gyz}, {gxz, gyz, gzz}};

    // plain superposition of the two metrics
    g[0][0] += g_bs[0][0] - 1.0;
    g[1][1] += g_bs[1][1] - 1.0;
    g[2][2] += g_bs[2][2] - 1.0;

    // Define the conformal factor
    double det_g = g[0][0] * (g[1][1] * g[2][2] - g[1][2] * g[2][1]) -
                 g[0][1] * (g[1][0] * g[2][2] - g[1][2] * g[2][0]) +
                 g[0][2] * (g[1][0] * g[2][1] - g[1][1] * g[2][0]);
    double chi = pow(det_g, -1. / 3.);

    // Initialise conformal metric
    FOR2(i, j) h[i][j] = g[i][j] * chi;

     // Define initial conformal factor
    vars.chi += chi;

    // Define initial lapse
    vars.lapse += lapse_;

    // Define initial trace of K and A_ij
    FOR2(i, j) vars.h[i][j] = h[i][j];
    vars.K = 0.;
    FOR2(i, j) vars.A[i][j] = 0.;

    current_cell.store_vars(vars);
}

#endif /* TSBOSONSTARTEUKOLSKYWAVE_IMPL_HPP_ */
