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
#include "BHBSBinaryParams.hpp"
#include "BosonStarSolution.hpp"
#include "WeightFunction.hpp"
#include <vector>
#include "parstream.H" //gives pout

#ifdef USE_TWOPUNCTURES
#include "TwoPunctures.hpp"
#endif

//! Class which solves for the initial data for a spherically symmetric boson
//! star with phi^4 coupling
class BHBSBinary
{

public:
    //! The constructor
    BHBSBinary(BosonStar_params_t a_params_BosonStar, BlackHole_params_t a_params_BlackHole,
                Binary_params_t a_params_Binary, Potential::params_t a_params_potential, 
                double a_G_Newton, double a_dx, int a_verbosity);
    
    #ifdef USE_TWOPUNCTURES
    //! Constructor in case of Two Punctures being used
    BHBSBinary(BosonStar_params_t a_params_BosonStar, BlackHole_params_t a_params_BlackHole,
                Binary_params_t a_params_Binary, Potential::params_t a_params_potential, 
                double a_G_Newton, double a_dx, int a_verbosity, TP::TwoPunctures &a_two_punctures);
    #endif

    //!  This function computes the 1d solution for the BS in the binary
    void compute_1d_BS_solution(const double max_r);

    //! Function to compute the value of all the initial vars on the grid
    template <class data_t>
    void compute(Cell<data_t> current_cell) const;

    //!The object that stores the solution found by the 1d ODE integrator
    BosonStarSolution m_1d_sol;

    #ifdef USE_TWOPUNCTURES
    //! The TwoPunctures object
    TP::TwoPunctures *m_two_punctures;
    #endif

protected:
    double m_dx;
    double m_G_Newton;
    BosonStar_params_t m_params_BosonStar;
    BlackHole_params_t m_params_BlackHole;
    Binary_params_t m_params_Binary;
    Potential::params_t m_params_potential; //!< The potential params
    int m_verbosity;
};

#include "BHBSBinary.impl.hpp"

#endif /* BHBS_HPP_ */
