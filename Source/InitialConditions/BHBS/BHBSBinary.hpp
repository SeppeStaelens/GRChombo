/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef BHBSBINARY_HPP_
#define BHBSBINARY_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "MatterCCZ4.hpp"
#include "ComplexScalarField.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total no. components
#include "VarsTools.hpp"
#include "simd.hpp"
#include "ComplexPotential.hpp"
#include "BosonStarParams.hpp"
#include "BosonStarSolution.hpp"
#include "WeightFunction.hpp"
#include <vector>
#include "parstream.H" //gives pout

//! Class which solves for the initial data for a spherically symmetric boson
//! star with phi^4 coupling
class BHBSBinary
{

public:
    //! The constructor
    BHBSBinary(BosonStar_params_t a_params_BosonStar, BlackHole_params_t a_params_BlackHole,
        Potential::params_t a_params_potential, double a_G_Newton, double a_dx, int a_verbosity);

    //! Computes the 1d solution and stores in m_1d_sol
    void compute_1d_solution(const double max_r);

    //! Function to compute the value of all the initial vars on the grid
    template <class data_t>
    void compute(Cell<data_t> current_cell) const;

    //The object that stores the solution found by the 1d ODE integrator */
    BosonStarSolution m_1d_sol;

protected:
    double m_dx;
    double m_G_Newton;
    BosonStar_params_t m_params_BosonStar;
    BlackHole_params_t m_params_BlackHole;
    Potential::params_t m_params_potential; //!< The potential params
    int m_verbosity;
};

#include "BHBSBinary.impl.hpp"

#endif /* BHBS_HPP_ */