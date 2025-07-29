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
#include "BosonStarParams.hpp"
#include "ComplexPotential.hpp"

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

        // ######################################
        //  Single Boson Star Solver Parameters
        // ######################################

        // Boson Star initial data params
        pp.load("central_amplitude_CSF",
                bosonstar_params.central_amplitude_CSF);
        pp.load("phase", bosonstar_params.phase, 0.0);
        pp.load("gridpoints", bosonstar_params.gridpoints, 1000000);
        pp.load("BS_solver_verbosity", bosonstar_params.BS_solver_verbosity,
                false);

        pp.load("star_centre", bosonstar_params.star_centre, center);

        // Potential params
        pp.load("scalar_mass", potential_params.scalar_mass, 1.0);
        pp.load("phi4_coeff", potential_params.phi4_coeff, 0.0);
        pp.load("solitonic", potential_params.solitonic, false);
        pp.load("sigma_solitonic", potential_params.sigma_soliton, 0.2);

        positionA[0] = bosonstar_params.star_centre[0];
        positionA[1] = bosonstar_params.star_centre[1];
        positionA[2] = bosonstar_params.star_centre[2];

        pout() << "Star A is at x-position " << positionA[0] << endl;
        pout() << "Star A is at y-position " << positionA[1] << endl;
        pout() << "Star A is at z-position " << positionA[2] << endl;

        // Star Tracking
        pp.load("do_star_track", do_star_track, false);
        pp.load("number_of_stars", number_of_stars, 1);
        pp.load("star_points", star_points, 81);
        pp.load("star_track_width_A", star_track_width_A, 4.);
        pp.load("star_track_width_B", star_track_width_B, 4.);
        pp.load("direction_of_motion", star_track_direction_of_motion);
        pp.load("star_track_level", star_track_level, 5);

#ifdef USE_AHFINDER
        pp.load("AH_1_initial_guess", AH_1_initial_guess,
                0.5 * bosonstar_params.mass);
        pp.load("AH_set_origins_to_punctures", AH_set_origins_to_punctures,
                false);
#endif

        // Mass extraction
        pp.load("activate_mass_extraction", activate_mass_extraction, 0);
        pp.load("mass_write_extraction",
                mass_extraction_params.write_extraction, false);
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

        // Weyl extraction
        pp.load("activate_gw_extraction", activate_weyl_extraction, 0);

        // Effective potential extraction
        pp.load("activate_effective_potential_extraction",
                activate_effective_potential_extraction, 0);
        pp.load("effective_potential_write_extraction",
                effective_potential_extraction_params.write_extraction, false);
        std::string effective_potential_extraction_path;
        if (pp.contains("effective_potential_extraction_subpath"))
        {
        pp.load("effective_potential_extraction_subpath", effective_potential_extraction_path);
        if (!effective_potential_extraction_path.empty() && effective_potential_extraction_path.back() != '/')
                effective_potential_extraction_path += "/";
        if (output_path != "./" && !output_path.empty())
                effective_potential_extraction_path = output_path + effective_potential_extraction_path;
        }
        else
                effective_potential_extraction_path = data_path;
        pp.load("effective_potential_extraction_file_prefix",
                effective_potential_extraction_params.extraction_file_prefix,
                std::string("Veff"));
        pp.load("num_effective_potential_extraction_radii",
                effective_potential_extraction_params.num_extraction_radii, 2);
        double min_r, max_r;
        int effective_potential_extraction_level;
        pp.load("effective_potential_extraction_level",
                effective_potential_extraction_level, 0);
        pp.load("effective_potential_min_r", min_r, 10.);
        pp.load("effective_potential_max_r", max_r, 11.);
        std::vector<double> radii(
            effective_potential_extraction_params.num_extraction_radii);
        std::vector<int> levels(
            effective_potential_extraction_params.num_extraction_radii);
        for (int i = 0;
             i < effective_potential_extraction_params.num_extraction_radii;
             ++i)
        {
            radii[i] = min_r + i * (max_r - min_r) /
                                   (effective_potential_extraction_params
                                        .num_extraction_radii -
                                    1);
            levels[i] = effective_potential_extraction_level;
        }
        effective_potential_extraction_params.extraction_radii = radii;
        effective_potential_extraction_params.extraction_levels = levels;
        pp.load("num_points_phi_effective_potential",
                effective_potential_extraction_params.num_points_phi, 2);
        pp.load("num_points_theta_effective_potential",
                effective_potential_extraction_params.num_points_theta, 4);
        pp.load("effective_potential_extraction_center",
                effective_potential_extraction_params.extraction_center,
                {0.5 * L, 0.5 * L, 0.5 * L});
	
	// Phi extraction
	pp.load("activate_phi_extraction", activate_phi_extraction, false);
	
	if (activate_phi_extraction)
        {
	    phi_extraction_params.extraction_radii = effective_potential_extraction_params.extraction_radii;
	    phi_extraction_params.extraction_levels = effective_potential_extraction_params.extraction_levels;
            phi_extraction_params.write_extraction = effective_potential_extraction_params.write_extraction;
	    phi_extraction_params.center = effective_potential_extraction_params.center;
	    pp.load("num_points_phi_Phi", phi_extraction_params.num_points_phi, 2);
            pp.load("num_points_theta_Phi", phi_extraction_params.num_points_theta, 5);
            if (extraction_params.num_points_theta % 2 == 0)
            {
                extraction_params.num_points_theta += 1;
                pout() << "Parameter: num_points_theta incompatible with "
                          "Simpson's "
                       << "rule so increased by 1.\n";
            }

            if (pp.contains("modes"))
            {
                pp.load("num_modes", extraction_params.num_modes);
                std::vector<int> extraction_modes_vect(
                    2 * extraction_params.num_modes);
                pp.load("modes", extraction_modes_vect,
                        2 * extraction_params.num_modes);
                extraction_params.modes.resize(extraction_params.num_modes);
                for (int i = 0; i < extraction_params.num_modes; ++i)
                {
                    extraction_params.modes[i].first =
                        extraction_modes_vect[2 * i];
                    extraction_params.modes[i].second =
                        extraction_modes_vect[2 * i + 1];
                }
            }
            else
            {
                // by default extraction (l,m) = (2,0), (2,1) and (2,2)
                extraction_params.num_modes = 3;
                extraction_params.modes.resize(3);
                for (int i = 0; i < 3; ++i)
                {
                    extraction_params.modes[i].first = 2;
                    extraction_params.modes[i].second = i;
                }
            }

            std::string phi_extraction_path;
            if (pp.contains("phi_extraction_subpath"))
            {
                pp.load("phi_extraction_subpath", phi_extraction_path);
                if (!phi_extraction_path.empty() && phi_extraction_path.back() != '/')
                    phi_extraction_path += "/";
                if (output_path != "./" && !output_path.empty())
                    phi_extraction_path = output_path + phi_extraction_path;
            }
            else
                phi_extraction_path = data_path;

            phi_extraction_params.data_path = data_path;
            phi_extraction_params.extraction_path = phi_extraction_path;

            // default names to Weyl extraction
            pp.load("phi_extraction_file_prefix",
                    phi_extraction_params.extraction_file_prefix,
                    std::string("Phi_extraction_"));
            pp.load("phi_integral_file_prefix",
                    phi_extraction_params.integral_file_prefix,
                    std::string("Phi_mode_"));
        }

        // Do we cant to calculate L2 norms of constraint violations
        pp.load("calculate_constraint_violations",
                calculate_constraint_violations, false);

        // Do we want to calculate and write the Noether Charge to a file
        pp.load("calculate_noether_charge", calculate_noether_charge, false);

        // Variables for outputting inf-norm
        pp.load("num_vars_inf_norm", num_vars_inf_norm, 0);
        pp.load("vars_inf_norm", vars_inf_norm, num_vars_inf_norm, 0);
    }

    // Tagging thresholds
    Real regrid_threshold_phi, regrid_threshold_chi, regrid_threshold_rho;
    Real tag_radius_A, tag_radius_B, tag_buffer;

    std::array<int, 2> tag_punctures_max_levels;
    std::array<int, 2> tag_horizons_max_levels;

    // Initial data for matter and potential
    double G_Newton;
    bool identical; // whether or not the 2 boson stars have the same profile

    BosonStar_params_t bosonstar_params;
    Potential::params_t potential_params;

    // Mass extraction
    int activate_mass_extraction;
    extraction_params_t mass_extraction_params;

    int activate_weyl_extraction;

    int activate_effective_potential_extraction;
    extraction_params_t effective_potential_extraction_params;

    bool activate_phi_extraction;
    spherical_extraction_params_t phi_extraction_params;
    
    // Do we want to write a file with the L2 norms of contraints?
    bool calculate_constraint_violations;

    // Do we want to write the Noether Charge to a file
    bool calculate_noether_charge;

    // Vars for outputting in plot files
    // int num_plot_vars;
    // std::vector<int> plot_vars;

    // Vars for outputting inf-norms
    int num_vars_inf_norm;
    std::vector<int> vars_inf_norm;

    bool do_star_track;
    int number_of_stars;
    int star_points;
    double star_track_width_A;
    double star_track_width_B;
    std::string star_track_direction_of_motion;
    int star_track_level;

    std::array<double, CH_SPACEDIM> positionA, positionB;

#ifdef USE_AHFINDER
    double AH_1_initial_guess;
    bool AH_set_origins_to_punctures;
#endif

    //     int flux_extraction_level; // specifies times (level) to do angmom
    //     flux extraction bool do_flux_integration;
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
