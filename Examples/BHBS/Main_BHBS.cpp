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
#include "SetupFunctions.hpp"
#include "SimulationParameters.hpp"
#include "CallDoAnalysis.hpp"

// Problem specific includes:
#include "BHBSLevel.hpp"
#include "BHBSAMR.hpp"

int runGRChombo(int argc, char *argv[])
{
    // Load the parameter file and construct the SimulationParameter class
    // To add more parameters edit the SimulationParameters file.
    char *in_file = argv[1];
    GRParmParse pp(argc - 2, argv + 2, NULL, in_file);
    SimulationParameters sim_params(pp);

    // The line below selects the problem that is simulated
    // (To simulate a different problem, define a new child of AMRLevel
    // and an associated LevelFactory)
    BHBSAMR bhbs_amr;

    #ifdef USE_TWOPUNCTURES
    bhbs_amr.set_two_punctures_parameters(sim_params.tp_params);
    // Run TwoPunctures solver if id_choice is appropriate
	if (sim_params.bosonstar_params.id_choice == 6)
    {
        pout() << "Running TwoPunctures solver..." << endl;
        bhbs_amr.m_two_punctures.Run();	
        pout() << "TwoPunctures solver finished." << endl;
    }
    #endif
    
    // !!!! Seems like this assumes the existence of two boson stars...
    bhbs_amr.m_star_tracker.initial_setup(sim_params.do_star_track,
        sim_params.number_of_stars, {sim_params.positionA, sim_params.positionB},
        sim_params.star_points, sim_params.star_track_width_A, sim_params.star_track_width_B, sim_params.star_track_direction_of_motion);
    DefaultLevelFactory<BHBSLevel> bh_bs_level_fact(bhbs_amr, sim_params);
    setupAMRObject(bhbs_amr, bh_bs_level_fact);

    // Instantiate AMR interpolator for mass/GW extraction
    AMRInterpolator<Lagrange<4>> interpolator(
        bhbs_amr, sim_params.origin, sim_params.dx, sim_params.boundary_params,
        sim_params.verbosity);

    // Add a scheduler to GRAMR which just calls doAnalysis on every AMRLevel
    // at time 0. It is called later in postTimeStep
    RefCountedPtr<CallDoAnalysis> call_do_analysis_ptr(new CallDoAnalysis);
    RefCountedPtr<Scheduler> scheduler_ptr(new Scheduler);
    scheduler_ptr->schedule(call_do_analysis_ptr, sim_params.max_steps);
    bhbs_amr.schedule(scheduler_ptr);


    bhbs_amr.set_interpolator(&interpolator);

    using Clock = std::chrono::steady_clock;
    using Minutes = std::chrono::duration<double, std::ratio<60, 1>>;

    std::chrono::time_point<Clock> start_time = Clock::now();

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
