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
#include "EppleyPacketParams.hpp"

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
        pp.load("regrid_threshold_phi", regrid_threshold_phi, 0.3);
        pp.load("regrid_threshold_chi", regrid_threshold_chi, 0.3);

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

        // Eppley Packet initial data params
        pp.load("Eppley_amplitude", eppley_packet_params.amplitude, 0.01);
        pp.load("Eppley_sigma", eppley_packet_params.sigma, 0.1);
        pp.load("magnetic", eppley_packet_params.magnetic, 0);
        pp.load("time_offset", eppley_packet_params.time_offset, 0.);
        pp.load("wave_centre", eppley_packet_params.wave_centre, center);

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

        // Effective potential extraction
        pp.load("activate_effective_potential_extraction",
                activate_effective_potential_extraction, 0);
        pp.load("effective_potential_write_extraction",
                effective_potential_extraction_params.write_extraction, false);
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

        // Do we cant to calculate L2 norms of constraint violations
        pp.load("calculate_constraint_violations",
                calculate_constraint_violations, false);

        // Do we want to calculate and write the Noether Charge to a file
        pp.load("calculate_noether_charge", calculate_noether_charge, false);
    }

    // Tagging thresholds
    Real regrid_threshold_phi, regrid_threshold_chi;

    // Initial data for matter and potential
    double G_Newton;

    // Boson star parameters
    BosonStar_params_t bosonstar_params;
    Potential::params_t potential_params;
    
    // Teukolsky wave parameters
    EppleyPacket_params_t eppley_packet_params;

    // Mass extraction
    int activate_mass_extraction;
    extraction_params_t mass_extraction_params;

    // Effective potential extraction
    int activate_effective_potential_extraction;
    extraction_params_t effective_potential_extraction_params;

    // Do we want to write a file with the L2 norms of contraints?
    bool calculate_constraint_violations;

    // Do we want to write the Noether Charge to a file
    bool calculate_noether_charge;
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
