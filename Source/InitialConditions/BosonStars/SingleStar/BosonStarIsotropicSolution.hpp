/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef BOSONSTARISOTROPICSOLUTION_HPP_
#define BOSONSTARISOTROPICSOLUTION_HPP_

#include "BosonStarParams.hpp"
#include "BosonStarSolution.hpp"
#include "ComplexPotential.hpp"
#include <cmath>
#include <algorithm>
#include "SplineInterpolator.hpp"
#include <boost/numeric/odeint.hpp>
#include <limits>
#include "parstream.H" //gives pout
//#include <iomanip> //for setprecision()

//! This class constructs interpolation functions for the boson star solution
//! in isotropic coordinates given a rescaled solution in polar-areal coordinates

template <template<typename...> class initial_data_t, typename initial_state_t>
class BosonStarIsotropicSolution
{
public:
    tools::spline<initial_data_t> m_chi;
    tools::spline<initial_data_t> m_lapse;
    tools::spline<initial_data_t> m_phi;
    double m_frequency_over_mass;

    //! New constructor which can be called before a polar areal solution is
    //! computed
    BosonStarIsotropicSolution(BosonStar_params_t a_params_BosonStar,
        Potential::params_t a_params_potential, double a_G_Newton,
        int a_verbosity);

    //! If the new constructor is used, this function must be called afterwards
    //! to construct the isotropic solution from a polar-areal solution.
    void makeFromPolarArealSolution(
        BosonStarSolution<initial_data_t, initial_state_t>
        &a_polar_areal_solution, const double a_max_radius);

private:
    BosonStar_params_t m_params_BosonStar;
    Potential::params_t m_params_potential;
    initial_data_t<double> m_isotropic_grid = {};
    initial_data_t<double> m_polar_areal_grid = {};
    double m_G_Newton;
    int m_verbosity;

    //! This function calculates the isotropic grid from the polar areal grid.
    void calculateIsotropicGrid(
        BosonStarSolution<initial_data_t, initial_state_t>
        &a_polar_areal_solution, const double a_max_radius);

    //! Adds the points for the conformal factor chi interpolation function
    void construct_chi(BosonStarSolution<initial_data_t, initial_state_t>
        &a_polar_areal_solution);

    //! Adds the points for the lapse and scalar field interpolation functions
    //! and gets the frequency/mass from the polar-areal solution.
    void construct_phi_and_lapse(
        BosonStarSolution<initial_data_t, initial_state_t>
        &a_polar_areal_solution);
};

#include "BosonStarIsotropicSolution.impl.hpp"

#endif /* BOSONSTARISOTROPICSOLUTION_HPP_ */
