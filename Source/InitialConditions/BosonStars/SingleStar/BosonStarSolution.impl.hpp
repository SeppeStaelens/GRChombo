/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(BOSONSTARSOLUTION_HPP_)
#error "This file should only be included through BosonStarSolution.hpp"
#endif

#ifndef BOSONSTARSOLUTION_IMPL_HPP_
#define BOSONSTARSOLUTION_IMPL_HPP_

template <template<typename...> class initial_data_t, typename initial_state_t>
BosonStarSolution<initial_data_t, initial_state_t>::BosonStarSolution()
    : m_initial_grid {}, m_f_array {}, m_g_array {}, m_psi_array {},
    m_Psi_array {}, m_num_grid_points {0}, m_num_psi_roots {0},
    m_last_good_g_index {0}, m_frequency_over_mass {NAN}, m_ADM_mass {NAN} {}

template <template<typename...> class initial_data_t, typename initial_state_t>
void BosonStarSolution<initial_data_t, initial_state_t>::push_back(
    const initial_state_t a_vars, const double a_radius)
{
    m_initial_grid.push_back(a_radius);
    m_f_array.push_back(a_vars[0]);
    m_g_array.push_back(a_vars[1]);
    m_psi_array.push_back(a_vars[2]);
    m_Psi_array.push_back(a_vars[3]);
    ++m_num_grid_points;
}

template <template<typename...> class initial_data_t, typename initial_state_t>
void BosonStarSolution<initial_data_t, initial_state_t>::clear()
{
    m_initial_grid.clear();
    m_f_array.clear();
    m_g_array.clear();
    m_psi_array.clear();
    m_Psi_array.clear();
    m_num_grid_points = 0;
    m_num_psi_roots = 0;
    m_frequency_over_mass = NAN;
}

template <template<typename...> class initial_data_t, typename initial_state_t>
int BosonStarSolution<initial_data_t, initial_state_t>::get_num_grid_points()
const
{
    return m_num_grid_points;
}

template <template<typename...> class initial_data_t, typename initial_state_t>
double BosonStarSolution<initial_data_t, initial_state_t>::
    get_max_radius() const
{
    return m_initial_grid.back();
}

template <template<typename...> class initial_data_t, typename initial_state_t>
int BosonStarSolution<initial_data_t, initial_state_t>::
    find_inflection_index(const initial_data_t<double> &a_array) const
{
    //Start looking only in the outer 75% radius
    auto start_iterator = std::upper_bound(std::begin(m_initial_grid),
        std::end(m_initial_grid), 0.25 * m_initial_grid.back());
    const int start_index{static_cast<int>(std::distance(std::begin(
        m_initial_grid),start_iterator))};

    //first calculate the absolute approximate first order derivative
    initial_data_t<double> array_abs_diff(m_num_grid_points - start_index - 1);
    double temp_abs_diff;

    for(int i = start_index; i < m_num_grid_points - 1; ++i)
    {
        temp_abs_diff = std::abs((a_array[i+1] - a_array[i]) /
            ( m_initial_grid[i+1] - m_initial_grid[i] ) );
        array_abs_diff[i - start_index] =
            temp_abs_diff > std::numeric_limits<double>::epsilon() ?
            temp_abs_diff : HUGE_VAL;
    }

    //now get an iterator to the minimum element of this abs diff array
    auto iterator = std::min_element(std::begin(array_abs_diff),
        std::end(array_abs_diff));

    auto inflection_index =
        start_index + std::distance(std::begin(array_abs_diff), iterator);

    return inflection_index;
}

template <template<typename...> class initial_data_t, typename initial_state_t>
double BosonStarSolution<initial_data_t, initial_state_t>::
    get_psi_at_max_radius() const
{
    return m_psi_array.back();
}

template <template<typename...> class initial_data_t, typename initial_state_t>
void BosonStarSolution<initial_data_t, initial_state_t>::count_num_psi_roots()
{
    for (int i = 1; i < m_num_grid_points; ++i)
    {
        if(m_psi_array[i - 1] * m_psi_array[i] < 0)
        {
            ++m_num_psi_roots;
        }
    }
}

template <template<typename...> class initial_data_t, typename initial_state_t>
int BosonStarSolution<initial_data_t, initial_state_t>::get_num_psi_roots()
{
    if(m_num_psi_roots == 0)
        count_num_psi_roots();
    return m_num_psi_roots;
}

template <template<typename...> class initial_data_t, typename initial_state_t>
int BosonStarSolution<initial_data_t, initial_state_t>
    ::get_last_good_g_index()
{
    if(std::isnan(m_ADM_mass))
        calculate_ADM_mass();
    return m_last_good_g_index;
}

template <template<typename...> class initial_data_t, typename initial_state_t>
void BosonStarSolution<initial_data_t, initial_state_t>
    ::calculate_frequency_over_mass()
{
    //first calculate the grid values of the function for which the limit in the
    //far field is the frequency/m
    initial_data_t<double> frequency_limit_function;
    frequency_limit_function.reserve(m_num_grid_points);
    for(int i = 0; i < m_num_grid_points; ++i)
    {
        frequency_limit_function.push_back( std::exp( - m_f_array[i]
            - m_g_array[i] ) );
    }

    //Get the lower index between which the frequency should be evaluated
    int lower_eval_index = find_inflection_index(frequency_limit_function);

    //Use linear interpolation to evaluate the frequency limit function midway
    //between grid points at the point the derivative is smallest.
    m_frequency_over_mass = 0.5 * (frequency_limit_function[lower_eval_index]
        + frequency_limit_function[lower_eval_index + 1] );
}

template <template<typename...> class initial_data_t, typename initial_state_t>
void BosonStarSolution<initial_data_t, initial_state_t>::calculate_ADM_mass()
{
    //first calculate the grid values of the mass aspect function
    initial_data_t<double> mass_aspect_function;
    mass_aspect_function.reserve(m_num_grid_points);
    for(int i = 0; i < m_num_grid_points; ++i)
    {
        mass_aspect_function.push_back( 0.5 * m_initial_grid[i] *
            (1.0 - std::exp(-2.0 * m_g_array[i])) );
    }

    //Get the lower index between which the ADM mass should be evaluated
    m_last_good_g_index
        = find_inflection_index(mass_aspect_function);

    //Use linear interpolation to get the ADM mass
    m_ADM_mass = 0.5 * (mass_aspect_function[m_last_good_g_index]
        + mass_aspect_function[m_last_good_g_index + 1] );
}

template <template<typename...> class initial_data_t, typename initial_state_t>
double BosonStarSolution<initial_data_t, initial_state_t>
    ::get_frequency_over_mass()
{
    if(std::isnan(m_frequency_over_mass))
        calculate_frequency_over_mass();
    return m_frequency_over_mass;
}

template <template<typename...> class initial_data_t, typename initial_state_t>
double BosonStarSolution<initial_data_t, initial_state_t>::get_ADM_mass()
{
    if(std::isnan(m_ADM_mass))
        calculate_ADM_mass();
    return m_ADM_mass;
}

template <template<typename...> class initial_data_t, typename initial_state_t>
initial_data_t<double>& BosonStarSolution<initial_data_t, initial_state_t>::
    get_grid()
{
    return m_initial_grid;
}

template <template<typename...> class initial_data_t, typename initial_state_t>
initial_data_t<double>& BosonStarSolution<initial_data_t, initial_state_t>::
    get_f()
{
    return m_f_array;
}

template <template<typename...> class initial_data_t, typename initial_state_t>
initial_data_t<double>& BosonStarSolution<initial_data_t, initial_state_t>::
    get_g()
{
    return m_g_array;
}

template <template<typename...> class initial_data_t, typename initial_state_t>
initial_data_t<double>& BosonStarSolution<initial_data_t, initial_state_t>::
    get_psi()
{
    return m_psi_array;
}

template <template<typename...> class initial_data_t, typename initial_state_t>
initial_data_t<double>& BosonStarSolution<initial_data_t, initial_state_t>::
    get_Psi()
{
    return m_Psi_array;
}

#endif /* BOSONSTARSOLUTION_IMPL_HPP_ */
