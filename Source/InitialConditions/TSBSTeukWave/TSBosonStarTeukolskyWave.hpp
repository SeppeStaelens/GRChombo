/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef TSBOSONSTARTEUKOLSKYWAVE_HPP_
#define TSBOSONSTARTEUKOLSKYWAVE_HPP_

#include "BosonStarParams.hpp"
#include "Cell.hpp"
#include "ComplexPotential.hpp"
#include "ComplexScalarField.hpp"
#include "Coordinates.hpp"
#include "EppleyPacket.hpp"
#include "EppleyPacketParams.hpp"
#include "MatterCCZ4.hpp"
#include "Tensor.hpp"
#include "ThinShellSolution.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total no. components
#include "VarsTools.hpp"
#include "parstream.H" //gives pout
#include "simd.hpp"

//! Class which reads the initial data for a spherically symmetric thin shell
//! boson star solution and initialises the Eppley packet Teukolsky wave around it
template <class packet_t = EppleyPacketM0>
class TSBosonStarTeukolskyWave
{

  public:
    //! The constructor
    TSBosonStarTeukolskyWave(BosonStar_params_t a_params_BosonStar,
                    Potential::params_t a_params_potential, 
                    EppleyPacket_params_t a_params_eppley_packet,
                    double a_dx, double a_L);

    //! Computes the 1d solution and stores in m_1d_sol
    void read_1d_solution(const double max_r);

    //! Function to compute the value of all the initial vars on the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const;

    // The object that stores the metric functions
    packet_t m_eppley_packet; //!< The Eppley packet object

    // The object that stores the interpolators
    ThinShellSolution m_1d_sol;

  protected:
    double m_dx;
    BosonStar_params_t m_params_BosonStar;  //!< The complex scalar field params
    Potential::params_t m_params_potential; //!< The potential params
    EppleyPacket_params_t m_params_eppley_packet;  //!< The params for the Eppley packet
};

#include "TSBosonStarTeukolskyWave.impl.hpp"

#endif /* TSBOSONSTARTEUKOLSKYWAVE_HPP_ */
