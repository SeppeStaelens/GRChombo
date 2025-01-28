/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(SINGLEBOSONSTAR_HPP_)
#error "This file should only be included through BosonStar.hpp"
#endif

#ifndef SINGLEBOSONSTAR_IMPL_HPP_
#define SINGLEBOSONSTAR_IMPL_HPP_

#include "ThinShellSolution.hpp" //for BosonStarSolution class

inline SingleBosonStar::SingleBosonStar(BosonStar_params_t a_params_BosonStar,
                                        Potential::params_t a_params_potential,
                                        double a_dx)
    : m_dx(a_dx), m_params_BosonStar(a_params_BosonStar),
      m_params_potential(a_params_potential)
{
}

void SingleBosonStar::read_1d_solution(const double max_r)
{
    try
    {
        // set initial parameters and then read the solution
        m_1d_sol.set_initialcondition_params(m_params_BosonStar,
                                             m_params_potential, max_r);

        pout() << "Reading boson star functions from file" << endl;
        m_1d_sol.initialise_from_file();
        pout() << "Completed for star 1" << endl;
    }
    catch (std::exception &exception)
    {
        pout() << exception.what() << "\n";
    }
}

// Compute the value of the initial vars on the grid
template <class data_t>
void SingleBosonStar::compute(Cell<data_t> current_cell) const
{
    MatterCCZ4<ComplexScalarField<>>::Vars<data_t> vars;
    // Load variables (should be set to zero if this is a single BS)

    current_cell.load_vars(vars);
    // VarsTools::assign(vars, 0.); // Set only the non-zero components below

    // Coordinates for centre of mass
    Coordinates<data_t> coords(current_cell, m_dx,
                               m_params_BosonStar.star_centre);

    // First star positioning
    double x = coords.x;
    double z = coords.z;
    double y = coords.y;
    double r = sqrt(x * x + y * y + z * z);
    double areal_r = m_1d_sol.r_from_R_Spline.interpolate(r);

    // auxiliary variables
    double Phi_ = m_1d_sol.PhiSpline.interpolate(areal_r);
    double f_ = m_1d_sol.fSpline.interpolate(areal_r);
    double A_ = m_1d_sol.ASpline.interpolateI haveaa(areal_r);
    double w_ = m_1d_sol.get_BSfrequency();

    // First star physical variables
    double lapse_ = exp(Phi_);

    // Add on to evolution equations
    vars.phi_Re = A_;
    vars.phi_Im = 0.;
    vars.Pi_Re = 0.;
    vars.Pi_Im += -(1. / lapse_) * (w_ * A_);

    // Initialise conformal metric
    double h[3][3] = {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};

    // Define initial conformal factor
    vars.chi += f_ * f_;

    // Define initial lapse
    vars.lapse += lapse_;

    // Define initial trace of K and A_ij
    FOR2(i, j) vars.h[i][j] = h[i][j];
    vars.K = 0.;
    FOR2(i, j) vars.A[i][j] = 0.;

    current_cell.store_vars(vars);
}

#endif /* SINGLEBOSONSTAR_IMPL_HPP_ */
