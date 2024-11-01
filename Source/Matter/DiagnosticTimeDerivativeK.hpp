/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef NOETHERCHARGE_HPP_
#define NOETHERCHARGE_HPP_

#include "ADMConformalVars.hpp" // needed for CCz4 and matter variables
#include "CCZ4Geometry.hpp"
#include "Cell.hpp"
#include "ComplexScalarField.hpp"
#include "Coordinates.hpp"
#include "DimensionDefinitions.hpp"
#include "EMTensor.hpp"
#include "FourthOrderDerivatives.hpp"
#include "Interval.H"
#include "UserVariables.hpp"
#include "simd.hpp"

//! Calculates the Noether Charge integrand values and the modulus of the
//! complex scalar field on the grid
template <class deriv_t = FourthOrderDerivatives>
class DiagnosticTimeDerivativeK
{
    // Need matter variables and chi
    template <class data_t> using CCZ4Vars = CCZ4Vars::VarsNoGauge<data_t>;
    template <class data_t>
    using MatterVars = ComplexScalarField<>::Vars<data_t>;
    template <class data_t>
    using CCZ4VarsWithGauge = CCZ4Vars::VarsWithGauge<data_t>;
    template <class data_t>
    using Diff2Vars = CCZ4Vars::Diff2VarsWithGauge<data_t>;

  protected:
    const double m_G_Newton;
    const EMTensor<ComplexScalarField<>> m_emtensor;
    const deriv_t m_deriv;
    const double m_kappa1;
    const double m_kappa2;
    const m_formulation = USE_CCZ4;

  public:
    DiagnosticTimeDerivativeK(const double a_G_Newton = 1.0,
                              const EMTensor<ComplexScalarField<>> a_emtensor,
                              const double a_dx, const double a_kappa1,
                              const double a_kappa2)
        : m_G_Newton(a_G_Newton), m_emtensor(a_emtensor), m_deriv(a_dx),
          m_kappa1(a_kappa1), m_kappa2(a_kappa2)
    {
    }

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // load vars locally
        const auto ccz4_vars = current_cell.template load_vars<CCZ4Vars>();
        const auto matter_vars = current_cell.template load_vars<MatterVars>();
        const auto d1 = m_deriv.template diff1<CCZ4VarsWithGauge>(current_cell);
        const auto d2 = m_deriv.template diff2<Diff2Vars>(current_cell);

        const auto advec = this->m_deriv.template advection<CCZ4Vars>(
            current_cell, ccz4_vars.shift);

        data_t kappa1_times_lapse;
        if (m_params.covariantZ4)
            kappa1_times_lapse = m_kappa1;
        else
            kappa1_times_lapse = m_kappa1 * ccz4_vars.lapse;

        using namespace TensorAlgebra;

        auto h_UU = compute_inverse_sym(ccz4_vars.h);
        auto chris = compute_christoffel(d1.h, h_UU);

        Tensor<1, data_t> Z_over_chi;
        Tensor<1, data_t> Z;

        if (m_formulation == USE_BSSN)
        {
            FOR(i) Z_over_chi[i] = 0.0;
        }
        else
        {
            FOR(i)
            Z_over_chi[i] = 0.5 * (ccz4_vars.Gamma[i] - chris.contracted[i]);
        }
        FOR(i) Z[i] = ccz4_vars.chi * Z_over_chi[i];

        auto ricci = CCZ4Geometry::compute_ricci_Z(ccz4_vars, d1, d2, h_UU,
                                                   chris, Z_over_chi);

        data_t divshift = compute_trace(d1.shift);
        data_t Z_dot_d1lapse = compute_dot_product(Z, d1.lapse);
        data_t dlapse_dot_dchi = compute_dot_product(d1.lapse, d1.chi, h_UU);

        Tensor<2, data_t> covdtilde2lapse;
        Tensor<2, data_t> covd2lapse;
        FOR(k, l)
        {
            covdtilde2lapse[k][l] = d2.lapse[k][l];
            FOR(m)
            {
                covdtilde2lapse[k][l] -= chris.ULL[m][k][l] * d1.lapse[m];
            }
            covd2lapse[k][l] =
                ccz4_vars.chi * covdtilde2lapse[k][l] +
                0.5 * (d1.lapse[k] * d1.chi[l] + d1.chi[k] * d1.lapse[l] -
                       ccz4_vars.h[k][l] * dlapse_dot_dchi);
        }

        data_t tr_covd2lapse = -(GR_SPACEDIM / 2.0) * dlapse_dot_dchi;
        FOR(i)
        {
            tr_covd2lapse -= ccz4_vars.chi * chris.contracted[i] * d1.lapse[i];
            FOR(j)
            {
                tr_covd2lapse += h_UU[i][j] * (ccz4_vars.chi * d2.lapse[i][j] +
                                               d1.lapse[i] * d1.chi[j]);
            }
        }

        data_t dtK = advec.K +
                     ccz4_vars.lapse *
                         (ricci.scalar +
                          ccz4_vars.K * (ccz4_vars.K - 2 * ccz4_vars.Theta)) -
                     kappa1_times_lapse * GR_SPACEDIM * (1 + m_kappa2) *
                         ccz4_vars.Theta -
                     tr_covd2lapse;
        // ignoring the bit with cosmologivcal constant for now
        dtK += 4.0 * M_PI * m_G_Newton * ccz4_vars.lapse *
               (m_emtensor.S - 3 * m_emtensor.rho);

        // store the RHS of dtK as diagnostic variable
        current_cell.store_vars(dtK, c_dtK);
    }
};

#endif /* NOETHERCHARGE_HPP_ */
