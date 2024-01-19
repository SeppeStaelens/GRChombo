/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef BOSONSTARRHS_HPP_
#define BOSONSTARRHS_HPP_

#include "ComplexPotential.hpp"
#include <cmath>
#include <limits>

//! Class which contains the rescaled RHS of the initial data equations for a
//! spherically symmetric boson star

class BosonStarRHS
{
public:
    //! Constructor
    BosonStarRHS(Potential::params_t a_params_potential, double a_G_Newton);

    //! operator() which computes the RHS of the rescaled equations
    template <typename initial_state_t>
    void operator() (
        const initial_state_t &a_vars,  //!< the values of the initial variables
        initial_state_t &rhs_out,       //!< where the RHS is stored
        const double &a_radius);     //!< rescaled radial coordinate = mr

private:
    const double m_rescaled_phi4_coeff = 0.0;
    const double m_G_Newton;
};

#include "BosonStarRHS.impl.hpp"

#endif /* BOSONSTARRHS_HPP_ */
