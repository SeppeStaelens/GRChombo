/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(BINARYBS_HPP_)
#error "This file should only be included through BinaryBS.impl.hpp"
#endif

#ifndef BINARYBS_IMPL_HPP_
#define BINARYBS_IMPL_HPP_

inline BinaryBS::BinaryBS(BosonStar_params_t a_boson_star1_params,
                          BosonStar_params_t a_boson_star2_params,
                          Potential::params_t a_params_potential,
                          double a_G_Newton, double a_dx, bool a_identical,
                          int a_verbosity)
    : m_dx(a_dx), m_boson_star1(a_boson_star1_params, a_params_potential,
        a_G_Newton, a_dx, a_verbosity), m_boson_star2(a_boson_star2_params,
        a_params_potential, a_G_Newton, a_dx, a_verbosity),
        m_identical(a_identical), m_verbosity(a_verbosity)
{
}

void BinaryBS::compute_profiles(const double a_L)
{
    if(m_verbosity)
    {
        pout() << "BinaryBS::compute_profiles: Computing boson star 1 profile"
               << std::endl;
    }
    m_boson_star1.compute_1d_solution(3.5 * a_L);

    // only need to compute profile for 2nd star if different from first
    if(!m_identical)
    {
        if(m_verbosity)
        {
            pout() << "BinaryBS::compute_profiles:"
                      " Computing boson star 2 profile" << std::endl;
        }
        m_boson_star2.compute_1d_solution(3.5 * a_L);
    }
    else
    {
        if(m_verbosity)
        {
            pout() << "BinaryBS::compute_profiles: Boson star 2 identical"
                      " to star 1; skipping profile computation" << std::endl;
        }

        // Copy boson_star1 into boson_star2 keeping phase and centre
        auto boson_star2_centre = m_boson_star2.m_params_BosonStar.star_centre;
        double boson_star2_phase = m_boson_star2.m_params_BosonStar.phase;
        m_boson_star2 = m_boson_star1;
        m_boson_star2.m_params_BosonStar.star_centre = boson_star2_centre;
        m_boson_star2.m_params_BosonStar.phase = boson_star2_phase;
    }
}

template <class data_t>
void BinaryBS::compute(Cell<data_t> current_cell) const
{
    MatterCCZ4<ComplexScalarField<>>::Vars<data_t> vars;
    // Load variables (should be set to zero before)
    current_cell.load_vars(vars);

    // Get two different Coordinate objects centred on each star
    Coordinates<data_t> coords1(current_cell, m_dx,
        m_boson_star1.m_params_BosonStar.star_centre);
    Coordinates<data_t> coords2(current_cell, m_dx,
        m_boson_star2.m_params_BosonStar.star_centre);

    // Get the distance from the centre of each star
    double r1 = coords1.get_radius();
    double r2 = coords2.get_radius();

    // CCZ4 variable superposition
    data_t chi1 = m_boson_star1.m_1d_sol.m_chi(r1);
    data_t chi2 = m_boson_star2.m_1d_sol.m_chi(r2);
    data_t lapse1 = m_boson_star1.m_1d_sol.m_lapse(r1);
    data_t lapse2 = m_boson_star2.m_1d_sol.m_lapse(r2);

    // 1/chi = 1/chi1 + 1/chi2 - 1
    vars.chi += (chi1 * chi2) / (chi1 + chi2 - chi1 * chi2);
    // lapse^2 = lapse1^2 + lapse2^2 - 1
    vars.lapse += sqrt(lapse1 * lapse1 + lapse2 * lapse2 - 1);

    // Conformal metric is flat
    FOR1(i) vars.h[i][i] += 1.;

    // Matter superposition
    double phase1 = m_boson_star1.m_params_BosonStar.phase;
    double phase2 = m_boson_star2.m_params_BosonStar.phase;
    double frequency1 = m_boson_star1.m_1d_sol.m_frequency_over_mass
                        * m_boson_star1.m_params_potential.scalar_mass;
    double frequency2 = m_boson_star2.m_1d_sol.m_frequency_over_mass
                        * m_boson_star2.m_params_potential.scalar_mass;
    data_t mod_phi1 = m_boson_star1.m_1d_sol.m_phi(r1);
    data_t mod_phi2 = m_boson_star2.m_1d_sol.m_phi(r2);
    vars.phi_Re += mod_phi1 * cos(phase1) + mod_phi2 * cos(phase2);
    vars.phi_Im += mod_phi1 * sin(phase1) + mod_phi2 * sin(phase2);
    vars.Pi_Re += frequency1 * mod_phi1 * sin(phase1) / lapse1
                + frequency2 * mod_phi2 * sin(phase2) / lapse2;
    vars.Pi_Im += -frequency1 * mod_phi1 * cos(phase1) / lapse1
                 - frequency2 * mod_phi2 * cos(phase2) / lapse2;

    // Store variables
    current_cell.store_vars(vars);
}

#endif /* BINARYBS_IMPL_HPP_ */
