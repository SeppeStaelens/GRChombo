/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef BOSONSTARINTEGRATOR_HPP_
#define BOSONSTARINTEGRATOR_HPP_

#include "BosonStarParams.hpp"
#include "BosonStarRHS.hpp"
#include "BosonStarSolutionObserver.hpp"
#include "BosonStarSolution.hpp"
#include <boost/numeric/odeint.hpp>
#include "parstream.H" //for pout


//! Class that uses odeint to integrate the equations for a static, spherically
//! symmetric boson star in polar-areal coordinates.
template <template<typename...> class initial_data_t, typename initial_state_t>
class BosonStarIntegrator
{
public:
    //! Constructor
    BosonStarIntegrator(BosonStar_params_t a_params_BosonStar,
        Potential::params_t a_params_potential, double a_G_Newton, int a_verbosity);

    //! Do integration
    void doIntegration(const double a_f_central);

    //! Returns the solution in a BosonStarSolution object.
    BosonStarSolution<initial_data_t, initial_state_t>& getSolution();

protected:
    BosonStar_params_t m_params_BosonStar;
    BosonStarRHS m_boson_star_rhs;
    BosonStarSolution<initial_data_t, initial_state_t> m_boson_star_solution;
    BosonStarSolutionObserver<initial_data_t, initial_state_t> m_sol_observer;
    int m_verbosity;
};

#include "BosonStarIntegrator.impl.hpp"

#endif /* BOSONSTARINTEGRATOR_HPP_ */
