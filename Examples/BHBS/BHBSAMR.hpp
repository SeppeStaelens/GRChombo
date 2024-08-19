/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef BHBSAMR_HPP_
#define BHBSAMR_HPP_

// Even if USE_TWOPUNCTURES is not defined, this file will include BHAMR.hpp
#include "BHAMR.hpp"

#ifdef USE_TWOPUNCTURES
#include "TwoPunctures.hpp"
#endif

#include "StarTracker.hpp"

/// A descendent of Chombo's AMR class to interface with tools which require
/// access to the whole AMR hierarchy, and those of GRAMR
/**
 * This object inherits from BHAMR and combines functionalities from TPAMR and STAMR
 */
class BHBSAMR : public BHAMR // maybe inherit just from GRAMR?
{
  public:

    StarTracker m_star_tracker;
    #ifdef USE_TWOPUNCTURES
    TP::TwoPunctures m_two_punctures;
    #endif

    BHBSAMR() {}

    void set_interpolator(AMRInterpolator<Lagrange<4>> *a_interpolator) override
    {
        BHAMR::set_interpolator(a_interpolator);
        m_star_tracker.set_interpolator(a_interpolator);
    }

    #ifdef USE_TWOPUNCTURES
    void set_two_punctures_parameters(const TP::Parameters &params)
    {
        // explicitly invoke copy constructor of base Parameters class
        m_two_punctures.Parameters::operator=(params);
    }
    #endif
};

#endif /* BHBSAMR_HPP_ */
