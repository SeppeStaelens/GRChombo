/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// General includes common to most GR problems
#include "BinaryBSLevel.hpp"
#include "BoxLoops.hpp"
#include "NanCheck.hpp"
#include "PositiveChiAndAlpha.hpp"
#include "TraceARemoval.hpp"

// For RHS update
#include "MatterCCZ4.hpp"

// For constraints calculation
#include "MatterConstraints.hpp"
#include "ConstraintViolations.hpp"

// For tag cells
#include "ComplexPhiAndChiExtractionTaggingCriterion.hpp"

// Problem specific includes
#include "ComputePack.hpp"
#include "ComplexPotential.hpp"
#include "BinaryBS.hpp"
#include "ComplexScalarField.hpp"
#include "SetValue.hpp"

// For mass extraction
#include "ADMMass.hpp"
#include "MassExtraction.hpp"

// For GW extraction
#include "Weyl4.hpp"
#include "WeylExtraction.hpp"

// For Noether Charge calculation
#include "SmallDataIO.hpp"
#include "NoetherCharge.hpp"

// Things to do at each advance step, after the RK4 is calculated
void BinaryBSLevel::specificAdvance()
{
    // Enforce trace free A_ij and positive chi and alpha
    BoxLoops::loop(make_compute_pack(TraceARemoval(), PositiveChiAndAlpha()),
                   m_state_new, m_state_new, INCLUDE_GHOST_CELLS);

    // Check for nan's
    if (m_p.nan_check)
        BoxLoops::loop(NanCheck(), m_state_new, m_state_new, EXCLUDE_GHOST_CELLS,
                       disable_simd());
}

// Initial data for field and metric variables
void BinaryBSLevel::initialData()
{
    CH_TIME("BinaryBSLevel::initialData");
    if (m_verbosity)
        pout() << "BinaryBSLevel::initialData " << m_level << endl;

    // First initalise a BinaryBS object
    BinaryBS binary_bs(m_p.boson_star1_params, m_p.boson_star2_params,
                       m_p.potential_params, m_p.G_Newton, m_dx,
                       m_p.identical, m_verbosity);

    // Compute the profiles for each star
    binary_bs.compute_profiles(m_p.L);

    // First set everything to zero. Note that the BinaryBS compute class
    // assumes this is done
    BoxLoops::loop(make_compute_pack(SetValue(0.0), binary_bs),
                   m_state_new, m_state_new, INCLUDE_GHOST_CELLS,
                   disable_simd());
}

// Things to do before outputting a checkpoint file
void BinaryBSLevel::preCheckpointLevel()
{
    fillAllGhosts();
    Potential potential(m_p.potential_params);
    ComplexScalarFieldWithPotential complex_scalar_field(potential);
    BoxLoops::loop(make_compute_pack(
                    MatterConstraints<ComplexScalarFieldWithPotential>(
                    complex_scalar_field, m_dx, m_p.G_Newton), NoetherCharge()),
                   m_state_new, m_state_new, EXCLUDE_GHOST_CELLS);

}

// Things to do before outputting a plot file
void BinaryBSLevel::prePlotLevel()
{
    fillAllGhosts();
    Potential potential(m_p.potential_params);
    ComplexScalarFieldWithPotential complex_scalar_field(potential);
    BoxLoops::loop(make_compute_pack(
                    MatterConstraints<ComplexScalarFieldWithPotential>(
                    complex_scalar_field, m_dx, m_p.G_Newton), NoetherCharge()),
                   m_state_new, m_state_new, EXCLUDE_GHOST_CELLS);

}

// Things to do in RHS update, at each RK4 step
void BinaryBSLevel::specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                       const double a_time)
{
    // Enforce trace free A_ij and positive chi and alpha
    BoxLoops::loop(
        make_compute_pack(TraceARemoval(), PositiveChiAndAlpha()), a_soln,
        a_soln, INCLUDE_GHOST_CELLS);

    // Calculate MatterCCZ4 right hand side with matter_t = ComplexScalarField
    // We don't want undefined values floating around in the constraints so
    // zero these
    Potential potential(m_p.potential_params);
    ComplexScalarFieldWithPotential complex_scalar_field(potential);
    MatterCCZ4<ComplexScalarFieldWithPotential> my_ccz4_matter(
        complex_scalar_field, m_p.ccz4_params, m_dx, m_p.sigma, m_p.formulation,
        m_p.G_Newton);
    SetValue set_analysis_vars_zero(0.0, Interval(c_Pi_Im + 1, NUM_VARS - 1));
    auto compute_pack =
        make_compute_pack(my_ccz4_matter, set_analysis_vars_zero);
    BoxLoops::loop(compute_pack, a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
}

// Things to do at ODE update, after soln + rhs
void BinaryBSLevel::specificUpdateODE(GRLevelData &a_soln,
                                         const GRLevelData &a_rhs, Real a_dt)
{
    // Enforce trace free A_ij
    BoxLoops::loop(TraceARemoval(), a_soln, a_soln, INCLUDE_GHOST_CELLS);
}

// Things to do for analysis after each timestep and at the start
void BinaryBSLevel::doAnalysis()
{
    CH_TIME("BinaryBSLevel::specificPostTimeStep");

    fillAllGhosts();
    Potential potential(m_p.potential_params);
    ComplexScalarFieldWithPotential complex_scalar_field(potential);
    auto analysis_compute_pack =
        make_compute_pack(Weyl4(m_p.extraction_params.extraction_center, m_dx),
                          ADMMass(m_p.L, m_dx),
                          MatterConstraints<ComplexScalarFieldWithPotential>(
                          complex_scalar_field, m_dx, m_p.G_Newton),
                          NoetherCharge());
    BoxLoops::loop(analysis_compute_pack, m_state_new, m_state_new,
                   EXCLUDE_GHOST_CELLS);

    bool first_step = (m_dt == m_time);

    // Do the extraction on the min extraction level
    if (m_p.activate_gw_extraction == 1 &&
        m_level == m_p.extraction_params.min_extraction_level)
    {
        if (m_verbosity)
        {
            pout() << "BinaryBSLevel::specificPostTimeStep:"
                      " Extracting gravitational waves." << endl;
        }

        // Refresh the interpolator and do the interpolation
        m_gr_amr.m_interpolator->refresh();
        WeylExtraction gw_extraction(m_p.extraction_params, m_dt, m_time,
                                     first_step, m_restart_time);
        gw_extraction.execute_query(m_gr_amr.m_interpolator);
    }

    if (m_p.activate_mass_extraction == 1 &&
        m_level == m_p.mass_extraction_params.min_extraction_level)
    {
        if (m_verbosity)
        {
            pout() << "BinaryBSLevel::specificPostTimeStep:"
                      " Extracting mass." << endl;
        }

        // Now refresh the interpolator and do the interpolation
        m_gr_amr.m_interpolator->refresh();
        MassExtraction mass_extraction(m_p.mass_extraction_params, m_dt,
                                    m_time, m_restart_time,
                                    first_step);
        mass_extraction.execute_query(m_gr_amr.m_interpolator);
    }

    if (m_level == 0)
    {
        if (m_p.calculate_constraint_violations)
        {
            // Write constraint violations to file
            ConstraintViolations constraint_violations(c_Ham,
                Interval(c_Mom1, c_Mom3), &m_gr_amr, m_p.coarsest_dx, m_dt,
                m_time, m_restart_time, "ConstraintViolations.dat",
                first_step);
            constraint_violations.execute();
        }

        if (m_p.calculate_noether_charge)
        {
            // compute integrated volume weighted noether charge integral
            double noether_charge = m_gr_amr.compute_sum(c_N, m_dx);
            SmallDataIO noether_charge_file("NoetherCharge.dat", m_dt, m_time,
                                            m_restart_time,
                                            SmallDataIO::APPEND,
                                            first_step);
            noether_charge_file.remove_duplicate_time_data();
            if (m_time == 0.)
            {
                noether_charge_file.write_header_line({"Noether Charge"});
            }
            noether_charge_file.write_time_data_line({noether_charge});
        }

        // Compute the maximum of mod_phi and write it to a file
        double mod_phi_max = m_gr_amr.compute_max(
                                Interval(c_mod_phi, c_mod_phi));
        SmallDataIO mod_phi_max_file("mod_phi_max.dat", m_dt, m_time,
                                     m_restart_time,
                                     SmallDataIO::APPEND,
                                     first_step);
        mod_phi_max_file.remove_duplicate_time_data();
        if (m_time == 0.)
        {
            mod_phi_max_file.write_header_line({"max mod phi"});
        }
        mod_phi_max_file.write_time_data_line({mod_phi_max});

        // Calculate the infinity-norm of all variables specified in params file
        // and output them
        if (m_p.num_vars_inf_norm > 0)
        {
            pout() << "Variable infinity norms:\n";
            for (int icomp : m_p.vars_inf_norm)
            {
                std::string var_name = UserVariables::variable_names[icomp];
                double var_norm = m_gr_amr.compute_norm(Interval(icomp, icomp),
                                    0., m_p.coarsest_dx);
                pout() << var_name << ": " << var_norm << "\t";
            }
            pout() << std::endl;
        }
    }
}

// Specify if you want any plot files to be written, with which vars
void BinaryBSLevel::specificWritePlotHeader(
    std::vector<int> &plot_states) const
{
    plot_states = m_p.plot_vars;
}

void BinaryBSLevel::computeTaggingCriterion(FArrayBox &tagging_criterion,
                                               const FArrayBox &current_state)
{
    BoxLoops::loop(ComplexPhiAndChiExtractionTaggingCriterion(m_dx, m_level,
                   m_p.mass_extraction_params, m_p.regrid_threshold_phi,
                   m_p.regrid_threshold_chi), current_state, tagging_criterion);
}
