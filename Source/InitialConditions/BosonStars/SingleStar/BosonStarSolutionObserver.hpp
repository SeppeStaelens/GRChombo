/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef BOSONSTARSOLUTIONOBSERVER_HPP_
#define BOSONSTARSOLUTIONOBSERVER_HPP_

#include <cmath> //for std::abs
#include <stdexcept> //for std::exception
#include "BosonStarSolution.hpp" //BosonStarSolution Class

//! This is a class to store the solution found by odeint in storage arrays as
//! it is calculated.

template <template<typename...> class initial_data_t, typename initial_state_t>
class BosonStarSolutionObserver
{
public:
    BosonStarSolution<initial_data_t, initial_state_t> &m_boson_star_solution;
    int m_num_psi_roots; //!< this counts the number of roots in psi

    //! Constructor
    BosonStarSolutionObserver(
        BosonStarSolution<initial_data_t, initial_state_t>
        &a_boson_star_solution, int a_num_psi_roots = 0)
         : m_boson_star_solution(a_boson_star_solution),
         m_num_psi_roots(a_num_psi_roots) {}

    //! Overloaded () operator required for interface with boost odeint library
    void operator() (const initial_state_t &a_vars, double a_radius)
    {
        //don't do checks if this is the first step as the arrays will have
        //zero length.
        if (m_boson_star_solution.get_num_grid_points() > 0)
        {
            //first remind us which variable we're referring to
            auto f_new = a_vars[0];
            auto g_new = a_vars[1];
            auto psi_new = a_vars[2];
            auto Psi_new = a_vars[3];
            auto psi_old = m_boson_star_solution.get_psi_at_max_radius();
            auto psi_central = m_boson_star_solution.get_psi()[0];

            //if our solution has grown too large throw an exception
            if (std::abs(psi_new) > 2.0 * psi_central
            || std::abs(f_new) > 1.0e4 || std::abs(g_new) > 1.0e4 ||
            std::abs(Psi_new) > 1.0e4)
            {
                throw std::runtime_error("Solution blow up; "\
                "Integration stopped.");
            }

            //if psi has changed sign, increase our count of the number of roots
            if (psi_new * psi_old < 0.0)
            {
                ++m_num_psi_roots;
            }
        }
        //finally store the result of this step and carry on
        m_boson_star_solution.push_back(a_vars, a_radius);
    }
};

#endif /* BOSONSTARSOLUTIONOBSERVER_HPP_ */
