/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include <array>

#ifndef EPPLEYPACKETPARAMS_HPP_
#define EPPLEYPACKETPARAMS_HPP_

//! A structure for the input params for the boson star
struct EppleyPacket_params_t
{
    double sigma; //!< width of the packet
    double amplitude; //!< amplitude of the packet
    double magnetic; //!< magnetic quantum number
    std::array<double, CH_SPACEDIM>
        wave_centre; //!< coordinates of the centre of the star
};

#endif /* EPPLEYPACKETPARAMS_HPP_ */
