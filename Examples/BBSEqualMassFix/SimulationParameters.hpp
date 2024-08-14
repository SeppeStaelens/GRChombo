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
        pp.load("BS_solver_psc", bosonstar_params.PSC, 2.0);
        pp.load("BS_solver_omc", bosonstar_params.OMC, 0.5);
        pp.load("BS_solver_verbosity", bosonstar_params.BS_solver_verbosity,
                false);

        pp.load("star_centre", bosonstar_params.star_centre, center);

        // Potential params
        pp.load("scalar_mass", potential_params.scalar_mass, 1.0);
        pp.load("phi4_coeff", potential_params.phi4_coeff, 0.0);
        pp.load("solitonic", potential_params.solitonic, false);
        pp.load("sigma_solitonic", potential_params.sigma_solitonic, 0.2);

        // ######################################
        //  Binary Boson Star Parameters
        // ######################################

        pp.load("antiboson", bosonstar_params.antiboson, false);
        pp.load("BS_rapidity", bosonstar_params.BS_rapidity, 0.0);
        pp.load("BS_separation", bosonstar_params.BS_separation, 0.0);
        pp.load("BS_mass", bosonstar_params.mass, 1.0);
        pp.load("BS_impact_parameter", bosonstar_params.BS_impact_parameter,
                0.0);

        // Initialize values for bosonstar2_params to same as bosonstar_params
        // and then assign that ones that should differ below
        bosonstar2_params = bosonstar_params;

        // Boson Star 2 parameters
        pp.load("central_amplitude_CSF2",
                bosonstar2_params.central_amplitude_CSF);
        pp.load("BS_rapidity2", bosonstar2_params.BS_rapidity);
        pp.load("BS_mass2", bosonstar2_params.mass, 1.0);

        // std::array<double, CH_SPACEDIM> positionA, positionB;

        positionA[0] =
            (bosonstar_params.star_centre[0] +
             bosonstar_params.mass_ratio * bosonstar_params.BS_separation /
                 (bosonstar_params.mass_ratio + 1.));
        positionA[1] = bosonstar_params.star_centre[1] -
                       bosonstar_params.mass_ratio *
                           bosonstar_params.BS_impact_parameter /
                           (bosonstar_params.mass_ratio + 1.);
        positionA[2] = bosonstar_params.star_centre[2];

        positionB[0] = (bosonstar_params.star_centre[0] -
                        bosonstar_params.BS_separation /
                            (bosonstar_params.mass_ratio + 1.));
        positionB[1] = bosonstar_params.star_centre[1] +
                       bosonstar_params.BS_impact_parameter /
                           (bosonstar_params.mass_ratio + 1.);
        positionB[2] = bosonstar_params.star_centre[2];

        pout() << "Star A is at x-position " << positionA[0] << endl;
        pout() << "Star A is at y-position " << positionA[1] << endl;
        pout() << "Star A is at z-position " << positionA[2] << endl;

        pout() << "Star B is at x-position " << positionB[0] << endl;
        pout() << "Star B is at y-position " << positionB[1] << endl;
        pout() << "Star B is at z-position " << positionB[2] << endl;

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
        pp.load("AH_2_initial_guess", AH_2_initial_guess,
                0.5 * bosonstar2_params.mass);
        pp.load("AH_set_origins_to_punctures", AH_set_origins_to_punctures,
                false);
#endif

        // Tagging
        pp.load("regrid_threshold_phi", regrid_threshold_phi);
        pp.load("regrid_threshold_rho", regrid_threshold_rho, 0.5);
        pp.load("regrid_threshold_chi", regrid_threshold_chi);

        pp.load("tag_radius_A", tag_radius_A, 4.);
        pp.load("tag_radius_B", tag_radius_B, 4.);
        pp.load("tag_buffer", tag_buffer, 0.5);
        pp.load("tag_punctures_max_levels", tag_punctures_max_levels,
                {max_level, max_level});
        pp.load("tag_horizons_max_levels", tag_horizons_max_levels,
                {max_level, max_level});

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

        // Work out the minimum extraction level
        auto min_extraction_level_it =
            mass_extraction_params.min_extraction_level();

        // Do we cant to calculate L2 norms of constraint violations
        pp.load("calculate_constraint_violations",
                calculate_constraint_violations, false);

        // Do we want to calculate and write the Noether Charge to a file
        pp.load("calculate_noether_charge", calculate_noether_charge, false);
    }

    // Tagging thresholds
    Real regrid_threshold_phi, regrid_threshold_chi, regrid_threshold_rho;
    Real tag_radius_A, tag_radius_B, tag_buffer;

    std::array<int, 2> tag_punctures_max_levels;
    std::array<int, 2> tag_horizons_max_levels;

    // Initial data for matter and potential
    double G_Newton;

    BosonStar_params_t bosonstar_params;
    BosonStar_params_t bosonstar2_params;
    Potential::params_t potential_params;

    // Mass extraction
    int activate_mass_extraction;
    extraction_params_t mass_extraction_params;

    int activate_weyl_extraction;

    // Do we want to write a file with the L2 norms of contraints?
    bool calculate_constraint_violations;

    // Do we want to write the Noether Charge to a file
    bool calculate_noether_charge;

    // For tracking
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
    double AH_2_initial_guess;
    bool AH_set_origins_to_punctures;
#endif
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
