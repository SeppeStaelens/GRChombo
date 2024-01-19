/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef BINARYBSLEVEL_HPP_
#define BINARYBSLEVEL_HPP_

#include "DefaultLevelFactory.hpp"
#include "GRAMRLevel.hpp"
// Problem specific includes
#include "ComplexPotential.hpp"
#include "ComplexScalarField.hpp"

//!  A class for the evolution of a binary boson star system
/*!
     The class calculates the profiles of two spherically symmetric boson star
     models, superposes them and then evolves using the CCZ4 equations.
*/
class BinaryBSLevel : public GRAMRLevel
{
    friend class DefaultLevelFactory<BinaryBSLevel>;
    // Inherit the contructors from GRAMRLevel
    using GRAMRLevel::GRAMRLevel;

    // Typedef for scalar field
    typedef ComplexScalarField<Potential> ComplexScalarFieldWithPotential;

    //! Things to do at the end of the advance step, after RK4 calculation
    virtual void specificAdvance() override;

    //! Initialize data for the field and metric variables
    virtual void initialData() override;

    //! routines to do before outputing checkpoint file
    virtual void preCheckpointLevel() override;

    //! routines to do before outputing plot file
    virtual void prePlotLevel() override;

    //! RHS routines used at each RK4 step
    virtual void specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                 const double a_time) override;

    //! Things to do in UpdateODE step, after soln + rhs update
    virtual void specificUpdateODE(GRLevelData &a_soln,
                                   const GRLevelData &a_rhs,
                                   Real a_dt) override;

    //! Specify which variables to write at plot intervals
    virtual void specificWritePlotHeader(std::vector<int> &plot_states) const;

    //! Tell Chombo how to tag cells for regridding
    virtual void computeTaggingCriterion(FArrayBox &tagging_criterion,
                                         const FArrayBox &current_state)
                                         override;

    //! Things to do for analysis after each timestep and at the start
    virtual void doAnalysis() override;
};

#endif /* BINARYBSLEVEL_HPP_ */
