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
        // Gravitional constant
        pp.load("G_Newton", G_Newton, 1.0);

        // Eppley Packete initial data params
        pp.load("Eppley_amplitude", eppley_packet_params.amplitude, 0.01);
        pp.load("Eppley_sigma", eppley_packet_params.sigma, 0.1);
        pp.load("magnetic", eppley_packet_params.magnetic, 0);
	      pp.load("time_offset", eppley_packet_params.time_offset, 0.);
        pp.load("wave_centre", eppley_packet_params.wave_centre, center);

        // Do we cant to calculate L2 norms of constraint violations
        pp.load("calculate_constraint_violations",
                calculate_constraint_violations, false);
    }

    // Tagging thresholds
    // TO DO

    EppleyPacket_params_t eppley_packet_params;

    // Initial data for matter and potential
    double G_Newton;

    // Do we want to write a file with the L2 norms of contraints?
    bool calculate_constraint_violations;
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
