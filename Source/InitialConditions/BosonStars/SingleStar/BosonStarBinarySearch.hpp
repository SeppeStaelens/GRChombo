/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef BOSONSTARBINARYSEARCH_HPP_
#define BOSONSTARBINARYSEARCH_HPP_

#include "BosonStarSolution.hpp" //for BosonStarSolution class
#include "BosonStarIntegrator.hpp" //for inheritance from BosonStarIntegrator class
#include "ComplexPotential.hpp" //for Potential::params_t struct
#include "MayDay.H"

/*! Class that implements the binary search shooting algorithm to find the
static, spherically-symmetric, ground state boson star solutions. It uses a
standard interval bisection method to find the central value of
f = (1/2)log(g_tt) - log(frequency/m) for a given central scalar field
amplitude.
*/

template <template<typename...> class initial_data_t, typename initial_state_t>
class BosonStarBinarySearch : \
    public BosonStarIntegrator<initial_data_t, initial_state_t>
{
public:
    /*
    //! Constructor which requires pre-computed bounding solutions sol_min and
    //! sol_max - not used anymore
    BosonStarBinarySearch(BosonStar_params_t a_params_BosonStar,
        Potential::params_t a_params_potential,
        BosonStarSolution<initial_data_t, initial_state_t> &a_sol_min,
        BosonStarSolution<initial_data_t, initial_state_t> &a_sol_max);
    */

    //! Constructor which calls findInterval to find sol_min and sol_max
    BosonStarBinarySearch(BosonStar_params_t a_params_BosonStar,
        Potential::params_t a_params_potential, double a_G_Newton,
        int a_verbosity, const double a_f_central_guess = -0.5);

    /*! Function called by the constructor to check the values of
    f_central_min and f_central_max bound the desired value of
    f_central. This checks that the sign of psi at the edge of the
    computational domain differs between the two solutions and that they have
    the correct number of roots (0 for solution corresponding to
    f_central_max and 1 for solution corresponding to f_central min) so
    that we are getting the ground state
    */
    bool checkValidInterval();

    //! Function to find bounding solutions sol_min and sol_max
    void findInterval();

    /*! This function performs the binary search shooting algorithm to find a
    value of f that is within binary_search_tol of the true value. It
    stores the final computed solution in m_sol_mid.
    */
    void shoot();

    /*! This function returns the final "mid" solution to the caller.
    */
    BosonStarSolution<initial_data_t, initial_state_t>& getShootedSolution();

private:
    double m_f_central_min; //!< Central value of f for the "min" solution
    double m_f_central_max; //!< Central value of f for the "max" solution
    double m_f_central_mid; //!< Central value of f for the "mid" solution

    BosonStarSolution<initial_data_t, initial_state_t> m_sol_min; //!< Stores the "min" solution.
    BosonStarSolution<initial_data_t, initial_state_t> m_sol_max; //!< Stores the "max" solution.
    BosonStarSolution<initial_data_t, initial_state_t> m_sol_mid; //!< Stores the "mid" solution.
};

#include "BosonStarBinarySearch.impl.hpp"

#endif /* BOSONSTARBINARYSEARCH_HPP_ */
