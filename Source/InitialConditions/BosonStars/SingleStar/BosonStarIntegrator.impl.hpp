/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(BOSONSTARINTEGRATOR_HPP_)
#error "This file should only be included through BosonStarIntegrator.hpp"
#endif

#ifndef BOSONSTARINTEGRATOR_IMPL_HPP_
#define BOSONSTARINTEGRATOR_IMPL_HPP_

template <template<typename...> class initial_data_t, typename initial_state_t>
BosonStarIntegrator<initial_data_t, initial_state_t>::BosonStarIntegrator(
    BosonStar_params_t a_params_BosonStar,
    Potential::params_t a_params_potential, double a_G_Newton, int a_verbosity)
    : m_params_BosonStar(a_params_BosonStar),
    m_boson_star_rhs(a_params_potential, a_G_Newton),
    m_sol_observer(m_boson_star_solution), m_verbosity(a_verbosity) {}

template <template<typename...> class initial_data_t, typename initial_state_t>
void BosonStarIntegrator<initial_data_t, initial_state_t>
    ::doIntegration(const double a_f_central)
{
    //First clear arrays to make sure nothing has been left from the last
    //integration.
    m_boson_star_solution.clear();

    //identify fixed BCs
    const double g_central{0.0};
    const double Psi_central{0.0};
    const double central_radius{0.0}; /*bit of a weird name but makes it easier
                                        to understand later on */
    //Set central BCs
    initial_state_t central_vars{a_f_central, g_central,
        m_params_BosonStar.central_amplitude_CSF, Psi_central};

    using namespace boost::numeric::odeint;
    typedef runge_kutta_dopri5<initial_state_t> error_stepper_t;
    //Need to put this in a try block as the solution observer can throw
    //exceptions.
    try
    {
        integrate_adaptive(make_controlled<error_stepper_t>
            (m_params_BosonStar.abs_error, m_params_BosonStar.rel_error),
            m_boson_star_rhs, central_vars, central_radius,
            m_params_BosonStar.max_radius, m_params_BosonStar.initial_step_size,
            m_sol_observer);
    }
    catch (std::exception &exception)
    {
        if (m_verbosity >= 2)
        {
            pout() << exception.what() << " Max areal radius = " <<
                m_boson_star_solution.get_grid().back() << " / m." << "\n";
        }
    }
}

template <template<typename...> class initial_data_t, typename initial_state_t>
BosonStarSolution<initial_data_t, initial_state_t>&
    BosonStarIntegrator<initial_data_t, initial_state_t>
    ::getSolution()
{
    return m_boson_star_solution;
}


#endif /* BOSONSTARINTEGRATOR_IMPL_HPP_ */
