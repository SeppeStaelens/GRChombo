/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef EFFECTIVEPOTENTIAL_HPP_
#define EFFECTIVEPOTENTIAL_HPP_

#include "CCZ4Vars.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "UserVariables.hpp" //This files needs c_NUM - total number of components

/*!
   The class allows the user to extract data from the grid to construct a proxy
   for the effective potential for null geodesics over spherical shells at
   specified radii. The values may then be written to an output file, or
   integrated across the surfaces. Formulas are based on Cunha et al. 2023, PRL
   130.6
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
        data_t R = coords.get_radius();
        data_t x = coords.x;
        data_t y = coords.y;
        data_t z = coords.z;
        data_t sint = sqrt(1. - z * z / (R * R));
        data_t cosp = x / (R * sint);
        data_t sinp = y / (R * sint);

        const auto vars = current_cell.template load_vars<Vars>();

        // The numerator of the effective potential is calculated by averaging
        // the lapse and shift squared over a sphere. The square root gets taken
        // in the Extraction
        const auto alpha_squared = vars.lapse * vars.lapse;

	auto norm_shift_squared = vars.chi;
	FOR2(i, j)
        norm_shift_squared +=
            vars.h[i][j] * vars.shift[i] * vars.shift[j] / vars.chi;
	norm_shift_squared -= vars.chi;

        // To account from the (numerical) departure of spherical symmetry, we
        // have to take the off-diagonal terms of the metric into account. The
        // following calculates the determinant of the restriction of the
        // spatial metric to the sphere. We can then calculate the proper
        // surface area of the sphere by integrating this over the domain of
        // theta and phi.
        const auto g_thth =
            z * z * (cosp * cosp * vars.h[0][0] + sinp * sinp * vars.h[1][1]) +
            R * R * sint * sint * vars.h[2][2] +
            2 * z *
                (cosp * vars.h[0][1] * z * sinp -
                 R * sint * (cosp * vars.h[0][2] + sinp * vars.h[1][2]));
        const auto g_phph = y * y * vars.h[0][0] + x * x * vars.h[1][1] -
                            2 * x * y * vars.h[0][1];
        const auto g_thph = -y * z * cosp * vars.h[0][0] +
                            x * z * sinp * vars.h[1][1] +
                            (x * cosp - y * sinp) * z * vars.h[0][1] +
                            R * sint * (y * vars.h[0][2] - x * vars.h[1][2]);
        const auto root_det_g =
            sqrt(g_thth * g_phph - g_thph * g_thph) / vars.chi;

        // As the volume element in the Spherical Extraction is hardcoded to be
        // r^2 sin theta, we need to divide by this.
        const auto V_eff_denominator = root_det_g / (R * R * sint);

        current_cell.store_vars(alpha_squared, c_alpha2);
        current_cell.store_vars(norm_shift_squared, c_beta2);
        current_cell.store_vars(V_eff_denominator, c_V_eff_d);
    }

  protected:
    const std::array<double, CH_SPACEDIM> m_center; //!< The grid center
    const double m_dx;                              //!< the grid spacing
};

#endif /* EFFECTIVEPOTENTIAL_HPP_ */
