/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(BOSONSTARBINARYSEARCH_HPP_)
#error "This file should only be included through BosonStarBinarySearch.hpp"
#endif

#ifndef BOSONSTARBINARYSEARCH_IMPL_HPP_
#define BOSONSTARBINARYSEARCH_IMPL_HPP_

/*
template <template<typename...> class initial_data_t, typename initial_state_t>
BosonStarBinarySearch<initial_data_t, initial_state_t>::BosonStarBinarySearch(
    BosonStar_params_t a_params_BosonStar,
    Potential::params_t a_params_potential,
    BosonStarSolution<initial_data_t, initial_state_t> &a_sol_min,
    BosonStarSolution<initial_data_t, initial_state_t> &a_sol_max)
    : BosonStarIntegrator<initial_data_t, initial_state_t>
    (a_params_BosonStar, a_params_potential),
    m_f_central_min(a_sol_min.get_f()[0]),
    m_f_central_max(a_sol_max.get_f()[0]),
    m_sol_min(a_sol_min), m_sol_max(a_sol_max), m_sol_mid(a_sol_min)
{
    if( !checkValidInterval() )
    {
        throw std::runtime_error("The interval provided does not bound a"\
        " ground state soluion. Please try again!");
    }
}
*/

template <template<typename...> class initial_data_t, typename initial_state_t>
BosonStarBinarySearch<initial_data_t, initial_state_t>::BosonStarBinarySearch(
    BosonStar_params_t a_params_BosonStar,
    Potential::params_t a_params_potential, double a_G_Newton, int a_verbosity,
    const double a_f_central_guess)
    : BosonStarIntegrator<initial_data_t, initial_state_t>
    (a_params_BosonStar, a_params_potential, a_G_Newton, a_verbosity),
    m_f_central_max(a_f_central_guess)
{
    findInterval();
    if( !checkValidInterval() )
    {
        MayDay::Error("BosonStarBinarySearch: " \
        "The found interval does not bound a ground state solution.");
    }
}

template <template<typename...> class initial_data_t, typename initial_state_t>
bool BosonStarBinarySearch<initial_data_t, initial_state_t>
    ::checkValidInterval()
{
    return( m_sol_min.get_psi_at_max_radius() * m_sol_max.get_psi_at_max_radius()
        < 0.0 && m_sol_min.get_num_psi_roots() == 1
        && m_sol_max.get_num_psi_roots() == 0 );
}

template <template<typename...> class initial_data_t, typename initial_state_t>
void BosonStarBinarySearch<initial_data_t, initial_state_t>::findInterval()
{
    int num_psi_roots_max, num_psi_roots_min, num_psi_roots_mid;
    this->doIntegration(m_f_central_max);
    num_psi_roots_max = (this->getSolution()).get_num_psi_roots();

    //First find a value of f_central_max that's close enough to 0 so that
    //psi in sol_max has no roots.
    while(num_psi_roots_max != 0)
    {
        m_f_central_max *= 0.5;
        this->doIntegration(m_f_central_max);
        num_psi_roots_max = (this->getSolution()).get_num_psi_roots();
    }
    m_sol_max = this->getSolution();

    //Now find a value of f_central_min for which psi in sol_min has at
    //at least one root.
    double f_central_step{0.125};
    m_f_central_min = m_f_central_max - f_central_step;
    this->doIntegration(m_f_central_min);
    num_psi_roots_min = (this->getSolution()).get_num_psi_roots();

    while(num_psi_roots_min == 0)
    {
        m_f_central_min -= f_central_step;
        this->doIntegration(m_f_central_min);
        num_psi_roots_min = (this->getSolution()).get_num_psi_roots();
    }
    m_sol_min = this->getSolution();

    //Now find a value of f_central_min for which psi in sol_min has exactly
    //one root
    while(num_psi_roots_min != 1)
    {
        m_f_central_mid = 0.5 * (m_f_central_max + m_f_central_min);
        this->doIntegration(m_f_central_mid);
        num_psi_roots_mid = (this->getSolution()).get_num_psi_roots();
        if(num_psi_roots_mid == 0)
        {
            m_f_central_max = m_f_central_mid;
            m_sol_max = this->getSolution();
        }
        else
        {
            m_f_central_min = m_f_central_mid;
            m_sol_min = this->getSolution();
            num_psi_roots_min = num_psi_roots_mid;

        }
    }
}

//Note that for some reason, C++ won't let you access member variables of a
//template base class without derefrencing the implicit this pointer. If this
//needs to be modified and any member function/variable of BosonStarIntegrator
//need to be accessed, remember to derefrence them from this*.
template <template<typename...> class initial_data_t, typename initial_state_t>
void BosonStarBinarySearch<initial_data_t, initial_state_t>::shoot()
{
    int n_binary_search_iter{0};

    //keep bisecting until we are within the desired tolerance
    while(m_f_central_max - m_f_central_min
        > this->m_params_BosonStar.binary_search_tol &&
        n_binary_search_iter <= this->m_params_BosonStar.max_binary_search_iter)
    {
        //Set central BCs for this iteration
        m_f_central_mid =
            0.5 * ( m_f_central_min + m_f_central_max );
        this->doIntegration(m_f_central_mid);

        //Copy the storage arrays into a solution object
        m_sol_mid = this->getSolution();

        if(m_sol_min.get_psi_at_max_radius()
            * m_sol_mid.get_psi_at_max_radius() > 0.0)
        {
            m_sol_min = m_sol_mid;
            m_f_central_min = m_f_central_mid;
        }
        else
        {
            m_sol_max = m_sol_mid;
            m_f_central_max = m_f_central_mid;
        }
    }

    //Calculate final mid solution.
    m_f_central_mid = 0.5 * ( m_f_central_min + m_f_central_max );
    this->doIntegration(m_f_central_mid);

    //Copy the storage arrays into a solution object
    m_sol_mid = this->getSolution();
}

template <template<typename...> class initial_data_t, typename initial_state_t>
BosonStarSolution<initial_data_t, initial_state_t>&
    BosonStarBinarySearch<initial_data_t, initial_state_t>::getShootedSolution()
{
    return m_sol_mid;
}

#endif /* BOSONSTARBINARYSEARCH_IMPL_HPP_ */
