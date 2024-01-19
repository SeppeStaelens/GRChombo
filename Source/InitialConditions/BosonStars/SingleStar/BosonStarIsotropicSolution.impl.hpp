/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(BOSONSTARISOTROPICSOLUTION_HPP_)
#error "This file should only be included through BosonStarIsotropicSolution.hpp"
#endif

#ifndef BOSONSTARISOTROPICSOLUTION_IMPL_HPP_
#define BOSONSTARISOTROPICSOLUTION_IMPL_HPP_

template <template<typename...> class initial_data_t, typename initial_state_t>
BosonStarIsotropicSolution<initial_data_t, initial_state_t>
    ::BosonStarIsotropicSolution(BosonStar_params_t a_params_BosonStar,
    Potential::params_t a_params_potential, double a_G_Newton, int a_verbosity)
    : m_params_BosonStar(a_params_BosonStar),
    m_params_potential(a_params_potential), m_G_Newton(a_G_Newton),
    m_verbosity(a_verbosity) {}

template <template<typename...> class initial_data_t, typename initial_state_t>
void BosonStarIsotropicSolution<initial_data_t, initial_state_t>
    ::makeFromPolarArealSolution(
    BosonStarSolution<initial_data_t, initial_state_t>
    &a_polar_areal_solution, const double a_max_radius)
{
    //pout() << setprecision(12);
    calculateIsotropicGrid(a_polar_areal_solution, a_max_radius);
    construct_chi(a_polar_areal_solution);
    construct_phi_and_lapse(a_polar_areal_solution);
}

template <template<typename...> class initial_data_t, typename initial_state_t>
void BosonStarIsotropicSolution<initial_data_t, initial_state_t>
    ::calculateIsotropicGrid(
    BosonStarSolution<initial_data_t, initial_state_t>
    &a_polar_areal_solution, const double a_max_radius)
{
    //Copy grid and calculate exp(g)
    initial_data_t<double> exp_g_array(
        a_polar_areal_solution.get_num_grid_points());
    initial_data_t<double> extended_grid = a_polar_areal_solution.get_grid();

    for(int i = 0; i < a_polar_areal_solution.get_num_grid_points(); ++i)
    {
        exp_g_array[i] = std::exp(a_polar_areal_solution.get_g()[i]);
    }

    const double rho_match{a_polar_areal_solution.get_grid()[
        a_polar_areal_solution.get_last_good_g_index()]};

    if (m_verbosity >= 2)
    {
        pout() << "Switched to asymptotics at areal radius = "
            << rho_match / m_params_potential.scalar_mass << "\n";
    }

    //if we want a solution up to a larger radius than is possible from the
    //polar areal solution, we can use asymptotics.
    if(rho_match < m_params_potential.scalar_mass * a_max_radius)
    {
        //Remove elements where the solution is no longer valid
        exp_g_array.resize(
            a_polar_areal_solution.get_last_good_g_index());
        extended_grid.resize(a_polar_areal_solution.get_last_good_g_index());

        //Now add in extra grid-points using the asymptotic solution
        const double max_computed_radius{extended_grid.back()};

        //Note the division by scalar_mass as a_max_radius is not rescaled
        //but max_computed_radius is rescaled
        int num_extra_points = std::ceil((m_params_potential.scalar_mass
            * a_max_radius - max_computed_radius)
            / m_params_BosonStar.initial_step_size);
        extended_grid.reserve(static_cast<int>(extended_grid.size())
            + num_extra_points);
        exp_g_array.reserve(extended_grid.size());

        double radius, exp_g_value;
        for(int i = 1; i <= num_extra_points; ++i)
        {
            radius = max_computed_radius
                + m_params_BosonStar.initial_step_size * i;
            //asymptotics = e^(-2g) ~ 1 - 2M/rho
            exp_g_value = std::pow(1.0 -
                2.0 * a_polar_areal_solution.get_ADM_mass() / radius, -0.5);
            extended_grid.push_back(radius);
            exp_g_array.push_back(exp_g_value);
        }
    }

    //Now construct interpolation function
    tools::spline<initial_data_t> exp_g_interpolated;
    exp_g_interpolated.set_points(extended_grid, exp_g_array);

    //RHS for converting to isotropic coordinates
    auto rhs_lambda = [&](const initial_state_t &R , initial_state_t &dRdr ,
        double r )
    {
        dRdr[0] = exp_g_interpolated(m_params_potential.scalar_mass * r)
            * R[0]/r;
    };

    //outer boundary condition
    initial_state_t R_max {0.25 * a_max_radius * (1.0 /
        exp_g_interpolated(m_params_potential.scalar_mass * a_max_radius)
        + 1.0) * (1.0 / exp_g_interpolated(m_params_potential.scalar_mass
        * a_max_radius) + 1.0)};

    //Since the ODE involves division by r, we can only integrate up to some
    //r > 0
    const double r_inner_limit {1.0e-10};

    using namespace boost::numeric::odeint;
    typedef runge_kutta_dopri5<initial_state_t> error_stepper_t;
    integrate_adaptive(make_controlled<error_stepper_t>
        (m_params_BosonStar.abs_error, m_params_BosonStar.rel_error),
        rhs_lambda, R_max, a_max_radius, r_inner_limit,
        -m_params_BosonStar.initial_step_size, //need to reverse sign as integrating inwards
        [&] (const initial_state_t &a_R, double &a_r)
        {
            m_isotropic_grid.push_back(a_R[0]);
            m_polar_areal_grid.push_back(a_r);
        });

    //reverse both arrays as they currently start from the outer limit
    std::reverse(m_isotropic_grid.begin(), m_isotropic_grid.end());
    std::reverse(m_polar_areal_grid.begin(), m_polar_areal_grid.end());

}

template <template<typename...> class initial_data_t, typename initial_state_t>
void BosonStarIsotropicSolution<initial_data_t, initial_state_t>
    ::construct_chi(BosonStarSolution<initial_data_t, initial_state_t>
    &a_polar_areal_solution)
{
    //first make an array to hold the isotropic grid values of chi.
    initial_data_t<double> chi_array(m_polar_areal_grid.size());

    for(int i = 0; i < static_cast<int>(m_polar_areal_grid.size()); i++)
    {
        chi_array[i] = m_isotropic_grid[i] * m_isotropic_grid[i] /
            ( m_polar_areal_grid[i] * m_polar_areal_grid[i] );
    }

    /*
    //For debugging
    pout() << setw(16) << "R" << "\t" << setw(16) << "chi" << "\n";
    for(int i = 0; i < static_cast<int>(m_isotropic_grid.size()); i++)
    {
        pout() << setw(16) << m_isotropic_grid[i] << "\t" << setw(16)
        << chi_array[i] << "\n";
    }
    */

    //construct interpolation function
    m_chi.set_points(m_isotropic_grid, chi_array);
}

template <template<typename...> class initial_data_t, typename initial_state_t>
void BosonStarIsotropicSolution<initial_data_t, initial_state_t>
    ::construct_phi_and_lapse(BosonStarSolution<initial_data_t, initial_state_t>
    &a_polar_areal_solution)
{
    //first construct interpolation function of f and psi = sqrt(4 pi G) * phi
    tools::spline<initial_data_t> f_interp;
    f_interp.set_points(a_polar_areal_solution.get_grid(),
        a_polar_areal_solution.get_f());
    tools::spline<initial_data_t> psi_interp;
    psi_interp.set_points(a_polar_areal_solution.get_grid(),
        a_polar_areal_solution.get_psi());

    //as for the isotropic grid, we'll want to switch to asymptotics for large radii
    const double rho_match{a_polar_areal_solution.get_grid()[
        a_polar_areal_solution.get_last_good_g_index()]};
    const double r_match{rho_match / m_params_potential.scalar_mass};


    //make arrays to hold the isotropic grid values of the
    //lapse = (frequency/m) * e^(f) and phi
    initial_data_t<double> lapse_array(m_polar_areal_grid.size());
    initial_data_t<double> phi_array(m_polar_areal_grid.size());

    double r, rho;
    double psi_factor = 1.0 / std::sqrt(4.0 * M_PI * m_G_Newton);
    m_frequency_over_mass = a_polar_areal_solution.get_frequency_over_mass();
    double phi_decay_rate
            = std::sqrt(1.0 - m_frequency_over_mass * m_frequency_over_mass);
    for(int i = 0; i < static_cast<int>(m_polar_areal_grid.size()); i++)
    {
        r = m_polar_areal_grid[i];
        rho = m_params_potential.scalar_mass * r;
        if(m_polar_areal_grid[i] <= r_match)
        {
            lapse_array[i] = m_frequency_over_mass * std::exp(f_interp(rho));
            phi_array[i] = psi_factor * psi_interp(rho);
        }
        else
        {
            //for large radii, lapse ~ sqrt(1 - 2M/rho)
            lapse_array[i] = std::sqrt(1.0 -
                2.0 * a_polar_areal_solution.get_ADM_mass() / rho);
            //for large radii, psi ~ K/rho * exp( -rho * sqrt(1 - frequency^2/m^2))
            phi_array[i] = psi_factor * psi_interp(rho_match) * std::exp(
                - (rho - rho_match) * phi_decay_rate) * rho_match / rho;
        }
    }

    /*
    //For debugging
    pout() << setw(16) << "R" << "\t" << setw(16) << "lapse" << "\t" << setw(16)
    << "phi" << "\n";
    for(int i = 0; i < static_cast<int>(m_isotropic_grid.size()); i++)
    {
        pout() << setw(16) << m_isotropic_grid[i] << "\t" << setw(16)
        << lapse_array[i] << "\t" << setw(16) << phi_array[i] << "\n";
    }
    */

    //construct interpolation functions
    m_lapse.set_points(m_isotropic_grid, lapse_array);
    m_phi.set_points(m_isotropic_grid, phi_array);
}

#endif /* BOSONSTARISOTROPICSOLUTION_IMPL_HPP_ */
