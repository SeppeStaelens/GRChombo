/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include <array>

#ifndef BHBSPARAMS_HPP_
#define BHBSPARAMS_HPP_

//! A structure for the input params for the boson star
struct BosonStar_params_t
{
    int gridpoints; // numer of gridpoints used to create boson star
    double central_amplitude_CSF; //!< Central amplitude of the star
    double phase;
    int eigen; // radial eigenstate of the boson star (0=ground)
    double binary_separation;
    double mass;
    double BS_impact_parameter;
    double Newtons_constant;
    double BlackHoleMass;
    double BS_rapidity;
    double mass_ratio;
    int n_power;
    int id_choice; // initial data choice: 0 - plain superposition, 1 - Thomas' trick, 2 - fixing conformal factor method
    double radius_width;
    int conformal_factor_power;
    bool antiboson;
    std::array<double, CH_SPACEDIM> star_centre; //!< coordinates of the centre of the star
};

#endif /* BOSONSTARPARAMS_HPP_ */
