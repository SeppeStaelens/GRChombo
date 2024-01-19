/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include <array>

#ifndef BOSONSTARPARAMS_HPP_
#define BOSONSTARPARAMS_HPP_

//! A structure for the input params for the boson star
struct BosonStar_params_t
{
    double central_amplitude_CSF; //!< Central amplitude of the star
    double phase;
    std::array<double, CH_SPACEDIM> star_centre; //!< coordinates of the centre of the star
    double abs_error; //!< absolute error tolerance for the 1D ODE integrator
    double rel_error; //!< relative error tolerance for the 1D ODE integrator
    double initial_step_size; //!< initial step size for the 1D ODE integrator
    double max_radius; /*!< the maximum (rescaled) radius the 1D ODE
    integrator will attempt to integrate up to. */
    double binary_search_tol; /*!< This is the tolerance to which the shooting
    parameter (the central value of f) is found. */
    int max_binary_search_iter = 1e3; //!< the maximum number of binary search iterations.
};

#endif /* BOSONSTARPARAMS_HPP_ */
