/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef EFFECTIVEPOTENTIAL_HPP_
#define EFFECTIVEPOTENTIAL_HPP_

#include "CCZ4Vars"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "UserVariables.hpp" //This files needs c_NUM - total number of components

/*!
   The class allows the user to extract data from the grid for the
   effective potential for null geodesics over spherical shells at specified
   radii. The values may then be written to an output file, or integrated across
   the surfaces.
*/

class EffectivePotential
{
  public:
    template <class data_t> using Vars = CCZ4Vars::VarsWithGauge<data_t>;

    //! The constructor
    EffectivePotential(const std::array<double, CH_SPACEDIM> a_center,
                       const double a_dx)
        : m_dx(a_dx), m_center(a_center)
    {
    }

    //! The compute member which calculates the effective potential at each
    //! point on the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        const Coordinates<data_t> coords(current_cell, m_dx, m_center);
        const double R = coords.get_radius();

        const auto vars = current_cell.template load_vars<Vars>();
        // const double V_eff = vars.lapse / (vars.chi * R);
        const double V_eff = 1.0;

        current_cell.store_vars(V_eff, c_V_eff);
    }

  protected:
    const std::array<double, CH_SPACEDIM> m_center; //!< The grid center
    const double m_dx;                              //!< the grid spacing
}

#endif /* EFFECTIVEPOTENTIAL_HPP_ */
