/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

// General includes
#include "GRParmParse.hpp"
#include "SimulationParametersBase.hpp"

// Problem specific includes:
#include "ComplexPotential.hpp"
#include "BosonStarParams.hpp"

class SimulationParameters : public SimulationParametersBase
{
public:
    SimulationParameters(GRParmParse &pp) : SimulationParametersBase(pp)
    {
        // Read the problem specific params
        readParams(pp);
    }

    void readParams(GRParmParse &pp)
    {
        // for regridding
        pp.load("regrid_threshold_phi", regrid_threshold_phi);
        pp.load("regrid_threshold_chi", regrid_threshold_chi);

        // Gravitional constant
        pp.load("G_Newton", G_Newton, 1.0);

        // Boson Star initial data params
        pp.load("central_amplitude_CSF",
                bosonstar_params.central_amplitude_CSF);
        pp.load("phase", bosonstar_params.phase, 0.0);
        pp.load("abs_error", bosonstar_params.abs_error, 1.0e-14);
        pp.load("rel_error", bosonstar_params.rel_error, 1.0e-14);
        pp.load("initial_step_size", bosonstar_params.initial_step_size,
                0.015625);
        pp.load("max_radius", bosonstar_params.max_radius, 256.);
        pp.load("binary_search_tol", bosonstar_params.binary_search_tol,
                1.0e-15);
        pp.load("max_binary_search_iter",
                bosonstar_params.max_binary_search_iter, 1000);
        pp.load("star_centre", bosonstar_params.star_centre,
                {0.5 * L, 0.5 * L, 0.5 * L});

        // Potential params
        pp.load("scalar_mass", potential_params.scalar_mass, 1.0);
        pp.load("phi4_coeff", potential_params.phi4_coeff, 0.0);

        // Mass extraction
        pp.load("activate_mass_extraction", activate_mass_extraction, 0);
        pp.load("num_mass_extraction_radii",
                mass_extraction_params.num_extraction_radii, 1);
        pp.load("mass_extraction_levels",
                mass_extraction_params.extraction_levels,
                mass_extraction_params.num_extraction_radii, 0);
        pp.load("mass_extraction_radii",
                mass_extraction_params.extraction_radii,
                mass_extraction_params.num_extraction_radii, 0.1);
        pp.load("num_points_phi_mass", mass_extraction_params.num_points_phi,
                2);
        pp.load("num_points_theta_mass",
                mass_extraction_params.num_points_theta, 4);
        pp.load("mass_extraction_center",
                mass_extraction_params.extraction_center,
                {0.5 * L, 0.5 * L, 0.5 * L});

        // Work out the minimum extraction level
        auto min_extraction_level_it =
            std::min_element(mass_extraction_params.extraction_levels.begin(),
                             mass_extraction_params.extraction_levels.end());
        mass_extraction_params.min_extraction_level = *(min_extraction_level_it);

        // Do we cant to calculate L2 norms of constraint violations
        pp.load("calculate_constraint_violations",
                calculate_constraint_violations, false);

        // Do we want to calculate and write the Noether Charge to a file
        pp.load("calculate_noether_charge", calculate_noether_charge, false);

        // Variables for outputting to plot files
        pp.load("num_plot_vars", num_plot_vars, 0);
        pp.load("plot_vars", plot_vars, num_plot_vars, 0);

        // Variables for outputting inf-norm
        pp.load("num_vars_inf_norm", num_vars_inf_norm, 0);
        pp.load("vars_inf_norm", vars_inf_norm, num_vars_inf_norm, 0);
    }

    // Tagging thresholds
    Real regrid_threshold_phi, regrid_threshold_chi;

    // Initial data for matter and potential
    double G_Newton;
    BosonStar_params_t bosonstar_params;
    Potential::params_t potential_params;

    // Mass extraction
    int activate_mass_extraction;
    extraction_params_t mass_extraction_params;

    // Do we want to write a file with the L2 norms of contraints?
    bool calculate_constraint_violations;

    // Do we want to write the Noether Charge to a file
    bool calculate_noether_charge;

    // Vars for outputting in plot files
    int num_plot_vars;
    std::vector<int> plot_vars;

    // Vars for outputting inf-norms
    int num_vars_inf_norm;
    std::vector<int> vars_inf_norm;
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
