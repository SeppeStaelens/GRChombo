/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// General includes common to most GR problems
#include "TeukolskyWaveLevel.hpp"
#include "BoxLoops.hpp"
#include "GammaCalculator.hpp"
#include "NanCheck.hpp"
#include "PositiveChiAndAlpha.hpp"
#include "TraceARemoval.hpp"
#include "ComputePack.hpp"
#include "SetValue.hpp"

// Tagging criterion
#include "ChiTaggingCriterion.hpp"

// For RHS update
#include "IntegratedMovingPunctureGauge.hpp"
#include "CCZ4RHS.hpp"

// For constraints calculation
#include "NewConstraints.hpp"

// Problem specific includes
#include "TeukolskyWave.hpp"

// For GW extraction
#include "WeylExtraction.hpp"
#include "Weyl4.hpp"

// DataIO
#include "SmallDataIO.hpp"

// for chombo grid Functions
#include "AMRReductions.hpp"
#include "InterpolationQuery.hpp"

// Things to do at each advance step, after the RK4 is calculated
void TeukolskyWaveLevel::specificAdvance()
{
    // Enforce trace free A_ij and positive chi and alpha
    BoxLoops::loop(make_compute_pack(TraceARemoval(), PositiveChiAndAlpha()),
                   m_state_new, m_state_new, INCLUDE_GHOST_CELLS);

    // Check for nan's
    if (m_p.nan_check)
        BoxLoops::loop(
            NanCheck(m_dx, m_p.center, "NaNCheck in specific Advance"),
            m_state_new, m_state_new, EXCLUDE_GHOST_CELLS, disable_simd());
}

// Initial data for field and metric variables
void TeukolskyWaveLevel::initialData()
{
    CH_TIME("TeukolskyWaveLevel::initialData");
    if (m_verbosity)
        pout() << "TeukolskyWaveLevel::initialData " << m_level << endl;

    // First initalise a BosonStar object
    TeukolskyWave teukolsky_wave(m_p.eppley_packet_params, m_dx);

    // First set everything to zero ... we don't want undefined values in
    // constraints etc, then  initial conditions for Boson Star
    BoxLoops::loop(make_compute_pack(SetValue(0.0), teukolsky_wave), m_state_new,
                   m_state_new, INCLUDE_GHOST_CELLS, disable_simd());

    BoxLoops::loop(GammaCalculator(m_dx), m_state_new, m_state_new,
                   EXCLUDE_GHOST_CELLS, disable_simd());

    fillAllGhosts();
}

// Things to do before outputting a checkpoint file
void TeukolskyWaveLevel::preCheckpointLevel()
{
    CH_TIME("TeukolskyWaveLevel::preCheckpointLevel");

    fillAllGhosts();
    BoxLoops::loop(Constraints(m_dx, c_Ham, Interval(c_Mom1, c_Mom3)),
        m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
}

// Things to do before outputting a plot file
void TeukolskyWaveLevel::prePlotLevel()
{
    CH_TIME("TeukolskyWaveLevel::prePlotLevel");

    fillAllGhosts();
    BoxLoops::loop(Constraints(m_dx, c_Ham, Interval(c_Mom1, c_Mom3)),
        m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
}

// Things to do in RHS update, at each RK4 step
void TeukolskyWaveLevel::specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                     const double a_time)
{
    // Enforce trace free A_ij and positive chi and alpha
    BoxLoops::loop(make_compute_pack(TraceARemoval(), PositiveChiAndAlpha()),
                   a_soln, a_soln, INCLUDE_GHOST_CELLS);

    BoxLoops::loop(
            CCZ4RHS<IntegratedMovingPunctureGauge, FourthOrderDerivatives>(
                m_p.ccz4_params, m_dx, m_p.sigma, m_p.formulation),
            a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
}

// Things to do at ODE update, after soln + rhs
void TeukolskyWaveLevel::specificUpdateODE(GRLevelData &a_soln,
                                       const GRLevelData &a_rhs, Real a_dt)
{
    // Enforce trace free A_ij
    BoxLoops::loop(TraceARemoval(), a_soln, a_soln, INCLUDE_GHOST_CELLS);
}

void TeukolskyWaveLevel::specificPostTimeStep()
{
    CH_TIME("TeukolskyWaveLevel::specificPostTimeStep");

    bool first_step = (m_time == 0.0);

    if (m_p.activate_weyl_extraction == 1)
    {
        int min_level = m_p.extraction_params.min_extraction_level();
        bool calculate_weyl = at_level_timestep_multiple(min_level);
        if (calculate_weyl)
        {
            // Populate the Weyl Scalar values on the grid
            fillAllGhosts();
            BoxLoops::loop(
                Weyl4(m_p.extraction_params.center, m_dx, m_p.formulation),
                m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);

            // Do the extraction on the min extraction level
            if (m_level == min_level)
            {
                CH_TIME("WeylExtraction");
                // Now refresh the interpolator and do the interpolation
                // fill ghosts manually to minimise communication
                bool fill_ghosts = false;
                m_gr_amr.m_interpolator->refresh(fill_ghosts);
                m_gr_amr.fill_multilevel_ghosts(
                    VariableType::diagnostic, Interval(c_Weyl4_Re, c_Weyl4_Im),
                    min_level);
                WeylExtraction my_extraction(m_p.extraction_params, m_dt,
                                             m_time, first_step,
                                             m_restart_time);
                my_extraction.execute_query(m_gr_amr.m_interpolator);
            }
        }
    }

    if (m_p.calculate_constraint_violations)
    {
        fillAllGhosts();
        BoxLoops::loop(Constraints(m_dx, c_Ham, Interval(c_Mom1, c_Mom3)),
                       m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
        if (m_level == 0)
        {
            AMRReductions<VariableType::diagnostic> amr_reductions(m_gr_amr);

        // constraints calculated pre check and pre plot so done here already
        double L2_Ham = amr_reductions.norm(c_Ham, 2, true);
        double L2_Mom = amr_reductions.norm(Interval(c_Mom1, c_Mom3), 2, true);
        double L1_Ham = amr_reductions.norm(c_Ham, 1, true);
        double L1_Mom = amr_reductions.norm(Interval(c_Mom1, c_Mom3), 1, true);

        std::string constraints_filename = m_p.data_path + "constraint_norms";
        SmallDataIO constraints_file(constraints_filename, m_dt, m_time,
                                     m_restart_time, SmallDataIO::APPEND,
                                     first_step);
        constraints_file.remove_duplicate_time_data();
        if (first_step)
        {
            constraints_file.write_header_line({"L^2_Ham", "L^2_Mom", "L^1_Ham",
                                                "L^1_Mom"});
        }
        constraints_file.write_time_data_line(
            {L2_Ham, L2_Mom, L1_Ham, L1_Mom});
        }
    }
}

void TeukolskyWaveLevel::computeTaggingCriterion(
    FArrayBox &tagging_criterion, const FArrayBox &current_state,
    const FArrayBox &current_state_diagnostics)
{
    BoxLoops::loop(ChiTaggingCriterion(m_dx), current_state, tagging_criterion);
}
