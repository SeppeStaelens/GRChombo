/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef CALL_DO_ANALYSIS_
#define CALL_DO_ANALYSIS_

#include "Scheduler.H"
#include "AMR.H"
#include "GRAMRLevel.hpp"

//! This is just an interface for the AMR scheduler to call doAnalysis on
//! every AMRLevel
class CallDoAnalysis : public Scheduler::PeriodicFunction
{
public:
    // Use default condstructor of PeriodicFunction
    using PeriodicFunction::PeriodicFunction;

    AMR* m_amr_ptr; //! pointer to AMR object

    virtual void setUp(AMR &a_AMR, int a_interval) override
    {
        m_amr_ptr = &a_AMR;
    }

    virtual void operator()(int a_step, Real a_time) override
    {
        auto amr_level_ptrs = m_amr_ptr->getAMRLevels().stdVector();
        // need to reverse this vector so that doAnalysis is called in order of
        // finest level to coarsest.
        std::reverse(std::begin(amr_level_ptrs), std::end(amr_level_ptrs));
        for (AMRLevel* amr_level_ptr : amr_level_ptrs)
        {
            GRAMRLevel::gr_cast(amr_level_ptr)->doAnalysis();
        }
    }
};

#endif /* CALL_DO_ANALYSIS_ */
