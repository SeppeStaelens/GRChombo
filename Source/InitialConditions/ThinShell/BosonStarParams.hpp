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
    int gridpoints;
    bool BS_solver_verbosity;
    double central_amplitude_CSF;
    double scalar_mass;
    double phi4_coeff;
    bool solitonic;
    double sigma_solitonic;
    double phase;
    double mass;
    int n_power;
    int conformal_factor_power;
    bool antiboson;
    std::array<double, CH_SPACEDIM>
        star_centre; //!< coordinates of the centre of the star
};

#endif /* BOSONSTARPARAMS_HPP_ */
