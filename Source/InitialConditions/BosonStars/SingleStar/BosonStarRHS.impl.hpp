/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(BOSONSTARRHS_HPP_)
#error "This file should only be included through BosonStarRHS.hpp"
#endif

#ifndef BOSONSTARRHS_IMPL_HPP_
#define BOSONSTARRHS_IMPL_HPP_

BosonStarRHS::BosonStarRHS(Potential::params_t a_params_potential,
    double a_G_Newton)
    : m_rescaled_phi4_coeff(a_params_potential.phi4_coeff
        / (4.0 * M_PI * a_G_Newton * a_params_potential.scalar_mass *
        a_params_potential.scalar_mass)), m_G_Newton(a_G_Newton) {}

template <typename initial_state_t>
void BosonStarRHS::operator() (const initial_state_t &a_vars,
initial_state_t &rhs_out, const double &a_radius)
{
    //introduce references to make equations easier to read
    auto &f = a_vars[0]; //f = (1/2)*log(g_tt) - log(frequency/m)
    auto &g = a_vars[1]; //g = (1/2)*log(g_rr)
    auto &psi = a_vars[2]; //psi = sqrt(4*pi*G)*|phi|
    auto &Psi = a_vars[3]; //Psi = d(psi)/d(radius)

    auto &f_rhs = rhs_out[0];
    auto &g_rhs = rhs_out[1];
    auto &psi_rhs = rhs_out[2];
    auto &Psi_rhs = rhs_out[3];

    //the RHS contains removable singularities so need to handle the case
    //radius = 0 separately
    if( a_radius < std::numeric_limits<double>::epsilon() )
    {
        f_rhs = 0.0;
        g_rhs = 0.0;
        psi_rhs = 0.0;
        Psi_rhs = 0.0;
    }
    else
    {
        //first calculate some quantities which are used more than once
        auto exp_minus2f = std::exp(-2.0 * f);
        auto exp_plus2g = std::exp(2.0 * g);
        auto psi2 = psi * psi;
        auto Psi2 = Psi * Psi;
        auto radius_inv = 1.0 / a_radius;

        f_rhs = 0.5 * (a_radius * (exp_plus2g * psi2 * (exp_minus2f - 1.0
            - 0.5 * m_rescaled_phi4_coeff * psi2) + Psi2)
            + radius_inv * (exp_plus2g - 1.0));
        g_rhs = 0.5 * (a_radius * (exp_plus2g * psi2 * (exp_minus2f + 1.0
            + 0.5 * m_rescaled_phi4_coeff * psi2) + Psi2)
            - radius_inv * (exp_plus2g -1.0));
        psi_rhs = Psi;
        Psi_rhs = (exp_plus2g * a_radius * psi2 * (1.0
            + 0.5 * m_rescaled_phi4_coeff * psi2) - radius_inv * (exp_plus2g +
            1.0)) * Psi + exp_plus2g * (1.0 + m_rescaled_phi4_coeff * psi2 -
            exp_minus2f) * psi;
    }
}

#endif /* BOSONSTARRHS_IMPL_HPP_ */
