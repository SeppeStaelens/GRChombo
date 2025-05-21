/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef TEUKOLSKYWAVE_HPP_
#define TEUKOLSKYWAVE_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "EppleyPacket.hpp"
#include "EppleyPacketParams.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total no. components
#include "VarsTools.hpp"
#include "parstream.H" //gives pout
#include "simd.hpp"

//! Class which initialises an Eppley packet Teukolsky wave
class TeukolskyWave
{

  public:
    //! The constructor
    TeukolskyWave(EppleyPacket_params_t a_params_eppley_packet, double a_dx);

    //! Function to compute the value of all the initial vars on the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const;

    // The object that stores the metric functions
    EppleyPacket m_eppley_packet; //!< The Eppley packet object

  protected:
    double m_dx;
    EppleyPacket_params_t m_params_eppley_packet;  //!< The params for the Eppley packet
};

#include "TeukolskyWave.impl.hpp"

#endif /* TEUKOLSKYWAVE_HPP_ */
