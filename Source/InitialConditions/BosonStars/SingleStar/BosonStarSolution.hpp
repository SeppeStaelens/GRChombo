/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef BOSONSTARSOLUTION_HPP_
#define BOSONSTARSOLUTION_HPP_

#include <cmath> //for std::exp and std::abs
#include <limits> //for double epsilon
#include <algorithm> //for std::min_element and std::upper_bound
#include <iterator> //for std::iterator

//! Container class which stores the grid values of the rescaled Boson Star
//! solutions in polar areal coordinates.

template <template<typename...> class initial_data_t, typename initial_state_t>
class BosonStarSolution
{
public:
    //! Constructor which just initialises the member arrays to be empty
    BosonStarSolution();

    //! Push back function is called by BosonStarSolutionObserver to push
    //! back the variables at every integration step
    void push_back(const initial_state_t a_vars, const double a_radius);

    //! Clears all arrays
    void clear();

    //! Returns the number of grid points
    int get_num_grid_points() const;

    //! Returns the maximum radius of the solution
    double get_max_radius() const;

    //! Approximately calculates the index of an array at which there is an
    //! inflection point. Uses a crude first order approximation.
    int find_inflection_index(const initial_data_t<double> &a_array) const;

    //! Returns the value of psi at the maximum radius
    double get_psi_at_max_radius() const;

    //! Returns the number of roots in Psi
    int get_num_psi_roots();

    //! Returns the last valid index for g
    int get_last_good_g_index();

    //! Calculates the complex oscillation frequency (divided by the scalar
    //! mass) of the solution and stores it in m_frequency
    void calculate_frequency_over_mass();

    //! Calculates the ADM mass of the solution
    void calculate_ADM_mass();

    //! Returns the frequency divided by the scalar mass
    double get_frequency_over_mass();

    //! Returns the ADM mass in rescaled units (M_pl^2/m)
    double get_ADM_mass();

    //! Returns m_initial_grid
    initial_data_t<double>& get_grid();

    //! Returns the m_f_array
    initial_data_t<double>& get_f();

    //! Returns the m_g_array
    initial_data_t<double>& get_g();

    //! Returns the m_psi_array
    initial_data_t<double>& get_psi();

    //! Returns the m_Psi_array
    initial_data_t<double>& get_Psi();

private:
    initial_data_t<double> m_initial_grid; //!< array holding the grid radial coordinates
    initial_data_t<double> m_f_array = {}; //!< array to hold f grid values
    initial_data_t<double> m_g_array = {}; //!< array to hold g grid values
    initial_data_t<double> m_psi_array = {}; //!< array to hold psi grid values
    initial_data_t<double> m_Psi_array = {}; //!< array to hold Psi grid values
    int m_num_grid_points; //!< this stores the number of grid values per variable
    int m_num_psi_roots; //!< this stores the number of roots in psi
    int m_last_good_g_index; /*!< this stores the largest grid index for which
                                g is still valid (before growing modes take
                                over) */
    double m_frequency_over_mass; //!< this stores the complex oscillation frequency/m
    double m_ADM_mass; //!< this stores the ADM mass in units of M_pl^2/m

    //! Function to count the number of roots in psi and store in m_num_psi_roots
    void count_num_psi_roots();
};

#include "BosonStarSolution.impl.hpp"

#endif /* BOSONSTARSOLUTION_HPP_ */
