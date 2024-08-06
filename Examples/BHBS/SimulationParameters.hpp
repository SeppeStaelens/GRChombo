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
#include "BHBSBinaryParams.hpp"
#include "AngMomFluxParams.hpp"
#include "BoostedBH.hpp"

#ifded USE_TWOPUNCTURES
#include "TP_Parameters.hpp"
#endif

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
	pp.load("regrid_threshold_rho", regrid_threshold_rho);
        pp.load("regrid_threshold_chi", regrid_threshold_chi);

        // Gravitional constant
        pp.load("G_Newton", G_Newton, 1.0);

        // Boson Star params
        pp.load("BS_mass", bosonstar_params.mass, 1.0);
        pp.load("BS_rapidity", bosonstar_params.rapidity, 0.0);

        pp.load("gridpoints",bosonstar_params.gridpoints, 400000);
        pp.load("central_amplitude_CSF",
                bosonstar_params.central_amplitude_CSF, 0.2);
        pp.load("phase", bosonstar_params.phase, 0.0);
        pp.load("eigen", bosonstar_params.eigen, 0);
        pp.load("antiboson", bosonstar_params.antiboson, false);

        pp.load("BS_radius_width", bosonstar_params.radius_width, 10.);
        pp.load("BS_bump_radius", bosonstar_params.bump_radius, 10.0);

        // Black Hole parameters
        pp.load("BH_mass", blackhole_params.mass, 1.0);
	pp.load("BH_rapidity", blackhole_params.rapidity, 0.0);

        pp.load("BH_radius_width", blackhole_params.radius_width, 10.0);
        pp.load("BH_bump_radius", blackhole_params.bump_radius, 10.0);   

        // Binary parameters
        pp.load("centre_of_mass", binary_params.centre_of_mass, center);
        pp.load("binary_separation", binary_params.separation, 16.0);
        pp.load("impact_parameter", binary_params.impact_parameter, 0.0);
        // pp.load("mass_ratio", binary_params.mass_ratio, 1.0);
        binary_params.mass_ratio = blackhole_params.mass / bosonstar_params.mass;
        pp.load ("do_rotation", bosonstar_params.do_rotation, false);
        
        pout() << "The BH/BS mass ratio is " << binary_params.mass_ratio << endl;

        pp.load("id_choice", binary_params.id_choice, 0);
        pp.load("epsilon", binary_params.epsilon, 0.1);
        pp.load("weight_function_choice", binary_params.weight_function_choice, 1);
	pp.load("weight_function_order", binary_params.weight_function_order, 4);

        pp.load("conformal_factor_power", binary_params.conformal_factor_power, -4);
        bosonstar_params.Newtons_constant = G_Newton;

        // Potential params
        pp.load("scalar_mass", potential_params.scalar_mass, 1.0);
        pp.load("phi4_coeff", potential_params.phi4_coeff, 0.0);
        pp.load("solitonic", potential_params.solitonic, false);
        pp.load("sigma_soliton", potential_params.sigma_soliton, 0.02);
        
	//std::array<double, CH_SPACEDIM> positionA, positionB;

	position_BS[0] = binary_params.centre_of_mass[0] + binary_params.mass_ratio * binary_params.separation / (binary_params.mass_ratio + 1.);
	position_BS[1] = binary_params.centre_of_mass[1] - binary_params.mass_ratio * binary_params.impact_parameter / (binary_params.mass_ratio + 1.);
	position_BS[2] = binary_params.centre_of_mass[2];

	position_BH[0] = binary_params.centre_of_mass[0] - binary_params.separation / (binary_params.mass_ratio + 1.);
	position_BH[1] = binary_params.centre_of_mass[1] + binary_params.impact_parameter / (binary_params.mass_ratio + 1.);
        position_BH[2] = binary_params.centre_of_mass[2];

	pout() << "Boson star is at x-position " << position_BS[0] << endl;
        pout() << "Boson star is at y-position " << position_BS[1] << endl;
        pout() << "Boson star is at z-position " << position_BS[2] << endl;	

	pout() << "Black hole is at x-position " << position_BH[0] << endl;
        pout() << "Black hole is at y-position " << position_BH[1] << endl;	
	pout() << "Black hole is at z-position " << position_BH[2] << endl;
	
	// Star Tracking
        pp.load("do_star_track", do_star_track, false);
        pp.load("number_of_stars", number_of_stars, 1);
        pp.load("star_points", star_points, 81);
        pp.load("star_track_width_BS", star_track_width_BS, 4.);
	pp.load("star_track_width_BH", star_track_width_BH, 4.);
        pp.load("direction_of_motion", star_track_direction_of_motion);
        pp.load("star_track_level", star_track_level, 5);
        
        //Tagging
        pp.load("tag_radius_BS", tag_radius_BS, 4.);
	pp.load("tag_radius_BH", tag_radius_BH, 4.);
        pp.load("tag_buffer", tag_buffer, 0.5);
	pp.load("tag_punctures_max_levels", tag_punctures_max_levels,
                {max_level, max_level});
        pp.load("tag_horizons_max_levels", tag_horizons_max_levels,
                {max_level, max_level});

        // Mass extraction
        pp.load("activate_mass_extraction", activate_mass_extraction, 0);
        pp.load("mass_write_extraction",
                mass_extraction_params.write_extraction,
                false);
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
        auto min_extraction_level_it = mass_extraction_params.min_extraction_level();

        // Do we cant to calculate L2 norms of constraint violations
        pp.load("calculate_constraint_violations",
                calculate_constraint_violations, false);

        // Do we want to calculate and write the Noether Charge to a file
        pp.load("calculate_noether_charge", calculate_noether_charge, false);

        // Variables for outputting to plot files
        //pp.load("num_plot_vars", num_plot_vars, 0);
        //pp.load("plot_vars", plot_vars, num_plot_vars, 0);

        // Variables for outputting inf-norm
        pp.load("num_vars_inf_norm", num_vars_inf_norm, 0);
        pp.load("vars_inf_norm", vars_inf_norm, num_vars_inf_norm, 0);


        pp.load("flux_extraction_level", flux_extraction_level, 0);
        /*pp.load("flux_number_of_radii", angmomflux_params.number_radii,1);
        pp.load("flux_do", angmomflux_params.do_flux_integration,false);
        pp.load("flux_extraction_level", angmomflux_params.extraction_level,0);
        pp.load("flux_num_theta", angmomflux_params.num_theta,10);
        pp.load("flux_num_phi", angmomflux_params.num_phi,10);
        pp.load("flux_extraction_centre", angmomflux_params.centre,
                                                {0.5 * L, 0.5 * L, 0.5 * L});

        angmomflux_params.radii.resize(angmomflux_params.number_radii);
        pp.load("flux_extraction_radii", angmomflux_params.radii,
                                                angmomflux_params.number_radii);*/
        #ifdef USE_AHFINDER
        pp.load("AH_initial_guess", AH_initial_guess, 0.5 * blackhole_params.mass);
        #endif
	
        #ifdef USE_TWOPUNCTURES

    	tp_params.verbose = (verbosity > 0);	

        double v_BS = -tanh(bosonstar_params.rapidity);
	double v_BH = tanh(blackhole_params.rapidity);

	double gamma_BS = 1 / sqrt ( 1 - v_BS * v_BS);
	double gamma_BH = 1 / sqrt ( 1 - v_BH * v_BH);

        double q = binary_params.mass_ratio;
        double d = binary_params.separation;
	double b = binary_params.impact_parameter;

        //amount the original configuration has been  rotated relative to x-axis parallel config. assumes impact param/ separation > 0. Needed to rotate BH momenta
        double rotation_angle = asin(b / sqrt(d * d + b * b ));
	 
	bool calculate_target_masses;
        pp.load("TP_calculate_target_masses", calculate_target_masses, false);
        tp_params.give_bare_mass = !calculate_target_masses;
	
	// masses
        if (calculate_target_masses)
        {
            // IDENTIFICATION: (+, BH1, BS) vs (-, BH2, BH)
            pp.load("TP_target_mass_plus", tp_params.target_M_plus, bosonstar_params.mass);
            pp.load("TP_target_mass_minus", tp_params.target_M_minus, blackhole_params.mass);
            pp.load("TP_adm_tol", tp_params.adm_tol, 1e-10);
            pout() << "The black holes have target ADM masses of "
                   << tp_params.target_M_plus << " and "
                   << tp_params.target_M_minus << "\n";
            bh1_params.mass = tp_params.target_M_plus;
            bh2_params.mass = tp_params.target_M_minus;
	}
	else
	{
            pp.load("TP_mass_plus", tp_params.par_m_plus);
            pp.load("TP_mass_minus", tp_params.par_m_minus);
            bh1_params.mass = tp_params.par_m_plus;
            bh2_params.mass = tp_params.par_m_minus;
            pout() << "The black holes have bare masses of "
                   << std::setprecision(16) << tp_params.par_m_plus << " and "
                   << tp_params.par_m_minus << "\n";
            // reset precision
            pout() << std::setprecision(6);
        }

        //magnitudes of momenta (in x-axis dir in final but not TP coord system)
        double p1 = bosonstar_params.mass * gamma_BS * v_BS;
        double p2 = blackhole_params.mass * gamma_BH * v_BH;
        std::array<double, CH_SPACEDIM> momentum_BS{p1, 0, 0};
	std::array<double, CH_SPACEDIM> momentum_BH{p2, 0, 0};
        if (bosonstar_params.do_rotation)
	{
	    pout() << "Doing TP coordinate rotation by angle " << rotation_angle << endl;
	    
	    //rotate into frame in which black holes are on x-axis as required for TwoPunctures; we'll undo when constructing initial data
	    momentum_BS[0] = cos(rotation_angle) *p1;
       	    momentum_BS[1] = sin(rotation_angle) *p1;

	    momentum_BH[0] = cos(rotation_angle) *p2;
            momentum_BH[1] = sin(rotation_angle) *p2;
	}

        // BH spin and momenta
        std::array<double, CH_SPACEDIM> spin_minus, spin_plus;
        pp.load("TP_momentum_plus", bh1_params.momentum, momentum_BS);
        pp.load("TP_momentum_minus", bh2_params.momentum, momentum_BH);
        pp.load("TP_spin_plus", spin_plus);
        pp.load("TP_spin_minus", spin_minus);
        FOR(i)
        {
            tp_params.par_P_plus[i] = bh1_params.momentum[i];
            tp_params.par_P_minus[i] = bh2_params.momentum[i];
            tp_params.par_S_minus[i] = spin_minus[i];
            tp_params.par_S_plus[i] = spin_plus[i];
        }

        pout() << "The corresponding momenta are:";
        pout() << "\nP_plus = ";
        FOR(i) { pout() << tp_params.par_P_plus[i] << " "; }
        pout() << "\nP_minus = ";
        FOR(i) { pout() << tp_params.par_P_minus[i] << " "; }

        pout() << "\nThe corresponding spins are:";
        pout() << "\nS_plus = ";
        FOR(i) { pout() << tp_params.par_S_plus[i] << " "; }
        pout() << "\nS_minus = ";
        FOR(i) { pout() << tp_params.par_S_minus[i] << " "; }
        pout() << "\n";

	// interpolation type
        bool use_spectral_interpolation;
        pp.load("TP_use_spectral_interpolation", use_spectral_interpolation,
                false);
        tp_params.grid_setup_method =
            (use_spectral_interpolation) ? "evaluation" : "Taylor expansion";

        // initial_lapse (default to psi^n)
        pp.load("TP_initial_lapse", tp_params.initial_lapse,
                std::string("psi^n"));
        if (tp_params.initial_lapse != "twopunctures-antisymmetric" &&
            tp_params.initial_lapse != "twopunctures-averaged" &&
            tp_params.initial_lapse != "psi^n" &&
            tp_params.initial_lapse != "brownsville")
        {
            std::string message = "Parameter: TP_initial_lapse: ";
            message += tp_params.initial_lapse;
            message += " invalid";
            MayDay::Error(message.c_str());
        }
        if (tp_params.initial_lapse == "psi^n")
        {
            pp.load("TP_initial_lapse_psi_exponent",
                    tp_params.initial_lapse_psi_exponent, -2.0);
        }

        // Spectral grid parameters
        pp.load("TP_npoints_A", tp_params.npoints_A, 30);
        pp.load("TP_npoints_B", tp_params.npoints_B, 30);
        pp.load("TP_npoints_phi", tp_params.npoints_phi, 16);
        if (tp_params.npoints_phi % 4 != 0)
        {
            MayDay::Error("TP_npoints_phi must be a multiple of 4");
        }

        // Solver parameters and tolerances
        pp.load("TP_Newton_tol", tp_params.Newton_tol, 1e-10);
        pp.load("TP_Newton_maxit", tp_params.Newton_maxit, 5);
        pp.load("TP_epsilon", tp_params.TP_epsilon, 1e-6);
        pp.load("TP_Tiny", tp_params.TP_Tiny, 0.0);
        pp.load("TP_Extend_Radius", tp_params.TP_Extend_Radius, 0.0);

        //total distance between BH and BS
	double total_sep = sqrt(d*d + b*b);

        // BH positions
        pp.load("TP_offset_minus", tp_offset_minus, - total_sep / (q + 1.));
        pp.load("TP_offset_plus", tp_offset_plus, total_sep * q / (q + 1.));
        bh1_params.center = center;
        bh2_params.center = center;
        bh1_params.center[0] += tp_offset_minus;
        bh2_params.center[0] += tp_offset_plus;
        double center_offset_x = 0.5 * (tp_offset_plus + tp_offset_minus);
        tp_params.center_offset[0] = center_offset_x;
        // par_b is half the distance between BH_minus and BH_plus
        tp_params.par_b = 0.5 * (tp_offset_plus - tp_offset_minus);
        pp.load("TP_swap_xz", tp_params.swap_xz, false);

        // Debug output
        pp.load("TP_do_residuum_debug_output",
                tp_params.do_residuum_debug_output, false);
        pp.load("TP_do_initial_debug_output", tp_params.do_initial_debug_output,
                false);

        // Irrelevant parameters set to default value
        tp_params.keep_u_around = false;
        tp_params.use_sources = false;
        tp_params.rescale_sources = true;
        tp_params.use_external_initial_guess = false;
        tp_params.multiply_old_lapse = false;
        tp_params.schedule_in_ADMBase_InitialData = true;
        tp_params.solve_momentum_constraint = false;
        tp_params.metric_type = "something else";
        tp_params.conformal_storage = "not conformal at all";
        tp_params.conformal_state = 0;
        tp_params.mp = 0;
        tp_params.mm = 0;
        tp_params.mp_adm = 0;
        tp_params.mm_adm = 0;
	
	#endif
    }

    // Tagging thresholds
    Real regrid_threshold_phi, regrid_threshold_chi, regrid_threshold_rho;
    Real tag_radius_BS, tag_radius_BH, tag_buffer;

    std::array<int, 2> tag_punctures_max_levels;
    std::array<int, 2> tag_horizons_max_levels;

    // Initial data for matter and potential
    double G_Newton;

    BosonStar_params_t bosonstar_params;
    BlackHole_params_t blackhole_params;
    Binary_params_t binary_params;
    Potential::params_t potential_params;

    // Mass extraction
    int activate_mass_extraction;
    extraction_params_t mass_extraction_params;

    int activate_weyl_extraction;

    // Do we want to write a file with the L2 norms of contraints?
    bool calculate_constraint_violations;

    // Do we want to write the Noether Charge to a file
    bool calculate_noether_charge;

    // Vars for outputting in plot files
    //int num_plot_vars;
    //std::vector<int> plot_vars;

    // Vars for outputting inf-norms
    int num_vars_inf_norm;
    std::vector<int> vars_inf_norm;

    bool do_star_track;
    int number_of_stars;
    int star_points;
    double star_track_width_BS;
    double star_track_width_BH;
    std::string star_track_direction_of_motion;
    int star_track_level;

    std::array<double, CH_SPACEDIM> position_BS, position_BH;

    int flux_extraction_level; // specifies times (level) to do angmom flux extraction

    #ifdef USE_TWOPUNCTURES
    double tp_offset_plus, tp_offset_minus;
   
    //param sets for TP data and each boosted BH
    TP::Parameters tp_params;
    BoostedBH::params_t bh2_params;
    BoostedBH::params_t bh1_params;
    #endif

    #ifdef USE_AHFINDER
    double AH_initial_guess;
    #endif
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
