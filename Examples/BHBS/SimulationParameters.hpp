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
#include "AngMomFluxParams.hpp"
#include  "BoostedBH.hpp"
#include "TP_Parameters.hpp"

class SimulationParameters : public SimulationParametersBase
{
public:
    SimulationParameters(GRParmParse &pp) : SimulationParametersBase(pp)
    {
        // Read the problem specific params
        readParams(pp);
    }

     #ifdef USE_TWOPUNCTURES
    double tp_offset_plus, tp_offset_minus;
   
    //param sets for TP data and each boosted BH
    TP::Parameters tp_params;
    BoostedBH::params_t bh2_params;
    BoostedBH::params_t bh1_params;
    #endif

    void readParams(GRParmParse &pp)
    {
        // for regridding
        pp.load("regrid_threshold_phi", regrid_threshold_phi);
	pp.load("regrid_threshold_rho", regrid_threshold_rho);
        pp.load("regrid_threshold_chi", regrid_threshold_chi);

        // Gravitional constant
        pp.load("G_Newton", G_Newton, 1.0);

        // Boson Star initial data params
        pp.load("central_amplitude_CSF",
                bosonstar_params.central_amplitude_CSF);
        pp.load("phase", bosonstar_params.phase, 0.0);
        pp.load("eigen", bosonstar_params.eigen, 0);
        pp.load("gridpoints",bosonstar_params.gridpoints,400000);

        pp.load("star_centre", bosonstar_params.star_centre, center);

        // Potential params
        pp.load("scalar_mass", potential_params.scalar_mass, 1.0);
        pp.load("phi4_coeff", potential_params.phi4_coeff, 0.0);
        pp.load("solitonic", potential_params.solitonic, false);
        pp.load("sigma_soliton", potential_params.sigma_soliton, 0.02);
        pp.load("antiboson", bosonstar_params.antiboson, false);
        pp.load("BS_rapidity", bosonstar_params.BS_rapidity, 0.0);
        pp.load("binary_separation", bosonstar_params.binary_separation, 0.0);
        pp.load("BS_mass", bosonstar_params.mass);
        pp.load("BS_impact_parameter", bosonstar_params.BS_impact_parameter, 0.0);
        pp.load("id_choice", bosonstar_params.id_choice, 2);
        pp.load("mass_ratio", bosonstar_params.mass_ratio, 1.0);
        pp.load("radius_width1", bosonstar_params.radius_width1, 10.);
        pp.load("radius_width2", bosonstar_params.radius_width2, 20.);
        pp.load("conformal_factor_power", bosonstar_params.conformal_factor_power, -4);
        pp.load("G_Newton", bosonstar_params.Newtons_constant, 1.0);
	pp.load("BS_bump_radius", bosonstar_params.BS_bump_radius, 10.0);
	pp.load("BH_bump_radius", bosonstar_params.BH_bump_radius, 10.0);   
        
        // BH parameters
        pp.load("BlackHoleMass", blackhole_params.BlackHoleMass, 0.);
	pp.load("BH_rapidity", blackhole_params.BH_rapidity, 0.0);
	pp.load("weight_function_order", blackhole_params.weight_function_order, 4);

        pp.load("epsilon", bosonstar_params.epsilon, 0.1);
        pp.load("weight_function_choice", bosonstar_params.weight_function_choice, 1);
        
	//std::array<double, CH_SPACEDIM> positionA, positionB;

	positionA[0] = (bosonstar_params.star_centre[0] + bosonstar_params.mass_ratio * bosonstar_params.binary_separation / (bosonstar_params.mass_ratio + 1.));
	positionA[1] = bosonstar_params.star_centre[1] - bosonstar_params.mass_ratio * bosonstar_params.BS_impact_parameter / (bosonstar_params.mass_ratio + 1.);
	positionA[2] = bosonstar_params.star_centre[2];

	positionB[0] = (bosonstar_params.star_centre[0] - bosonstar_params.binary_separation / (bosonstar_params.mass_ratio + 1.));
	positionB[1] = bosonstar_params.star_centre[1] + bosonstar_params.BS_impact_parameter / (bosonstar_params.mass_ratio + 1.);
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
        
        //Tagging
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
	
                                                #ifdef USE_TWOPUNCTURES

    	tp_params.verbose = (verbosity > 0);	
	 
	bool calculate_target_masses;
        pp.load("TP_calculate_target_masses", calculate_target_masses, false);
        tp_params.give_bare_mass = !calculate_target_masses;
	
	// masses
        if (calculate_target_masses)
        {
            pp.load("TP_target_mass_plus", tp_params.target_M_plus);
            pp.load("TP_target_mass_minus", tp_params.target_M_minus);
            pp.load("TP_adm_tol", tp_params.adm_tol, 1e-10);
            pout() << "The black holes have target ADM masses of "
                   << tp_params.target_M_plus << " and "
                   << tp_params.target_M_minus << "\n";
            bh1_params.mass = tp_params.target_M_minus;
            bh2_params.mass = tp_params.target_M_plus;
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

	// BH spin and momenta
        std::array<double, CH_SPACEDIM> spin_minus, spin_plus;
        pp.load("TP_momentum_minus", bh1_params.momentum);
        pp.load("TP_momentum_plus", bh2_params.momentum);
        pp.load("TP_spin_plus", spin_plus);
        pp.load("TP_spin_minus", spin_minus);
        FOR(i)
        {
            tp_params.par_P_minus[i] = bh1_params.momentum[i];
            tp_params.par_P_plus[i] = bh2_params.momentum[i];
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

        // BH positions
        pp.load("TP_offset_minus", tp_offset_minus);
        pp.load("TP_offset_plus", tp_offset_plus);
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
    Real tag_radius_A, tag_radius_B, tag_buffer;

    std::array<int, 2> tag_punctures_max_levels;
    std::array<int, 2> tag_horizons_max_levels;

    // Initial data for matter and potential
    double G_Newton;

    BosonStar_params_t bosonstar_params;
    BlackHole_params_t blackhole_params;
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
    double star_track_width_A;
    double star_track_width_B;
    std::string star_track_direction_of_motion;
    int star_track_level;

    std::array<double, CH_SPACEDIM> positionA, positionB;

    int flux_extraction_level; // specifies times (level) to do angmom flux extraction
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
