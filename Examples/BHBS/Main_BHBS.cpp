/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "parstream.H" //Gives us pout()
#include <iostream>
#include <chrono>

#include "DefaultLevelFactory.hpp"
#include "GRAMR.hpp"
#include "GRParmParse.hpp"
#include "MultiLevelTask.hpp"
#include "SetupFunctions.hpp"
#include "SimulationParameters.hpp"

// Problem specific includes:
#include "BHBSLevel.hpp"

// BHBSAMR object
#include "BHBSAMR.hpp"

int runGRChombo(int argc, char *argv[])
{
    // Load the parameter file and construct the SimulationParameter class
    // To add more parameters edit the SimulationParameters file.
    char *in_file = argv[1];
    GRParmParse pp(argc - 2, argv + 2, NULL, in_file);
    SimulationParameters sim_params(pp);

    pout() << "#################################" << endl;
    pout() << "# Parameters succesfully loaded #" << endl;
    pout() << "#################################" << endl;

    // The line below selects the problem that is simulated
    // (To simulate a different problem, define a new child of AMRLevel
    // and an associated LevelFactory)
    BHBSAMR bhbs_amr;

    #ifdef USE_TWOPUNCTURES
    /* Gareth's addition */
    bhbs_amr.set_two_punctures_parameters(sim_params.tp_params);
    // Run TwoPunctures solver if id_choice is appropriate
	if (sim_params.binary_params.id_choice == 6 &&  sim_params.restart_from_checkpoint == false)
        pout() << "Running TwoPunctures solver..." << endl;
        bhbs_amr.m_two_punctures.Run();
        pout() << "TwoPunctures solver finished." << endl;
    #endif

    // !!!! Seems like this assumes the existence of two boson stars...
    bhbs_amr.m_star_tracker.initialise_star_tracking(
        sim_params.do_star_track, sim_params.number_of_stars, 
        {sim_params.position_BS, sim_params.position_BH}, sim_params.star_points, 
        sim_params.star_track_width_BS, sim_params.star_track_width_BH, 
        sim_params.star_track_direction_of_motion);
    DefaultLevelFactory<BHBSLevel> bh_bs_level_fact(bhbs_amr, sim_params);

    setupAMRObject(bhbs_amr, bh_bs_level_fact);

    // Instantiate AMR interpolator for mass/GW extraction
    AMRInterpolator<Lagrange<4>> interpolator(
        bhbs_amr, sim_params.origin, sim_params.dx, sim_params.boundary_params,
        sim_params.verbosity);
    bhbs_amr.set_interpolator(&interpolator);

    #ifdef USE_AHFINDER
    if (sim_params.AH_activate)
    {
        AHSurfaceGeometry sphere(sim_params.position_BH);

        bhbs_amr.m_ah_finder.add_ah(sphere, sim_params.AH_initial_guess,
                                  sim_params.AH_params);
    }
    #endif

    using Clock = std::chrono::steady_clock;
    using Minutes = std::chrono::duration<double, std::ratio<60, 1>>;

    std::chrono::time_point<Clock> start_time = Clock::now();

    // Add a scheduler to call specificPostTimeStep on every AMRLevel at t=0
    auto task = [](GRAMRLevel *level)
    {
        if (level->time() == 0.)
            level->specificPostTimeStep();
    };
    // call 'now' really now
    MultiLevelTaskPtr<> call_task(task);
    call_task.execute(bhbs_amr);

    // Engage! Run the evolution.
    bhbs_amr.run(sim_params.stop_time, sim_params.max_steps);

    auto now = Clock::now();
    auto duration = std::chrono::duration_cast<Minutes>(now - start_time);
    pout() << "Total simulation time (mins): " << duration.count() << ".\n";

    bhbs_amr.conclude();

    // Write Chombo timer report
    CH_TIMER_REPORT();

    return 0;
}

int main(int argc, char *argv[])
{
    mainSetup(argc, argv);

    int status = runGRChombo(argc, argv);

    if (status == 0)
        pout() << "GRChombo finished." << std::endl;
    else
        pout() << "GRChombo failed with return code " << status << std::endl;

    mainFinalize();
    return status;
}
