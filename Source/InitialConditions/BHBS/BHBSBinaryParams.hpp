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
    // dynamical    
    double mass;                    //! mass of the boson star
    double rapidity;                //! rapidity of the boson star

    // BS solution construction
    int gridpoints;                 //! numer of gridpoints used to create boson star
    double central_amplitude_CSF;   //!< Central scalar field amplitude of the star
    double phase;                   //! Phase of the scalar field at the centre of the star
    int eigen;                      //! radial eigenstate of the boson star (0=ground)
    bool antiboson;                 //! whether to use an antiboson star
    double Newtons_constant;        //! Newton's constant
    
    // initial condition parameters
    double radius_width;            //! for initial conditions
    double bump_radius;             //! parameter determining the width of the bump functions in BS_BH initial data methods
};

//! A structure for the input params for the black hole
struct BlackHole_params_t
{
    // dynamical
    double mass;                    //! mass of the black hole
    double rapidity;                //! rapidity of the black hole

    // initial condition parameters -- can these be combined into one?
    double radius_width;            //! for initial conditions -- specify which
    double bump_radius;             //! parameter determining the width of the bump functions in BS_BH initial data methods
};

struct Binary_params_t
{
    // dynamical
    std::array<double, CH_SPACEDIM> centre_of_mass; //!< coordinates of the centre of mass
    double separation;              //! separation of the binary components (along x-axis)
    double impact_parameter;        //! impact parameter of the binary (along y-axis)
    double mass_ratio;              //! mass of the black hole divided by the mass of the boson star

    // initial condition parameters
    int id_choice; // initial data choice: 0 - plain superposition, 1 - Thomas' trick, 2 - fixing conformal factor method
    double epsilon;
    int weight_function_choice;
    int weight_function_order;

    // other
    int conformal_factor_power;
};

#endif /* BHBSPARAMS_HPP_ */
