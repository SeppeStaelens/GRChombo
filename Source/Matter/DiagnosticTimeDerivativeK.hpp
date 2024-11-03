/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef DIAGNOSTICTIMEDERIVATIVEK_HPP_
#define DIAGNOSTICTIMEDERIVATIVEK_HPP_

#include "ADMConformalVars.hpp" // needed for CCz4 and matter variables
#include "CCZ4Geometry.hpp"
#include "CCZ4RHS.hpp"
#include "Cell.hpp"
#include "ComplexScalarField.hpp"
#include "Coordinates.hpp"
#include "DimensionDefinitions.hpp"
#include "EMTensor.hpp"
#include "FourthOrderDerivatives.hpp"
#include "Interval.H"
#include "UserVariables.hpp"
#include "simd.hpp"
#include "CCZ4Vars.hpp"

//! Calculates the Noether Charge integrand values and the modulus of the
//! complex scalar field on the grid
template <class matter_t, class gauge_t = MovingPunctureGauge, class deriv_t = FourthOrderDerivatives>
class DiagnosticTimeDerivativeK
{
  public:
    template <class data_t> using MatterCCZ4Vars = typename MatterCCZ4<matter_t>::template Vars<data_t>;	
    using params_t = CCZ4_params_t<typename gauge_t::params_t>;
    // Need matter variables and chi
    template <class data_t>
    using CCZ4VarsWithGauge = CCZ4Vars::VarsWithGauge<data_t>;
    template <class data_t>
    using Diff2Vars = CCZ4Vars::Diff2VarsWithGauge<data_t>;
    template <class data_t> using MatterCCZ4RHSVars = typename MatterCCZ4RHS<matter_t>::template Vars<data_t>;

    enum{ USE_CCZ4, USE_BSSN};


  protected:
    const double m_G_Newton;
    const matter_t &m_matter;
    deriv_t m_deriv;
    double m_kappa1 = 0.1;
    double m_kappa2 = 0.;
    bool m_covariantZ4 = 1;
    const int m_formulation = USE_CCZ4;
    const int m_c_rho;
    const Interval m_c_Si;
    const Interval m_c_Sij;
    const int m_c_dtK;

  public:
    DiagnosticTimeDerivativeK(const double a_G_Newton,
                              const matter_t &a_matter,
                              const double a_dx, params_t a_params,
                              const int a_c_rho = -1, const Interval a_c_Si = Interval(), 
			      const Interval a_c_Sij = Interval(), const int a_c_dtK = -1)
        : m_G_Newton(a_G_Newton), m_matter(a_matter), m_deriv(a_dx), m_c_rho(a_c_rho), m_c_Si(a_c_Si), m_c_Sij(a_c_Sij), m_c_dtK(a_c_dtK)
    {
    m_kappa1 = a_params.kappa1;
    m_kappa2 = a_params.kappa2;
    m_covariantZ4 = a_params.covariantZ4;
    // add asserts from EMTensor
    }

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // load vars locally
        const auto ccz4_vars = current_cell.template load_vars<CCZ4VarsWithGauge>();
        const auto matter_vars = current_cell.template load_vars<MatterCCZ4Vars>();
        
	const auto d1 = m_deriv.template diff1<CCZ4VarsWithGauge>(current_cell);
	const auto d1_matter = m_deriv.template diff1<MatterCCZ4Vars>(current_cell);
        const auto d2 = m_deriv.template diff2<Diff2Vars>(current_cell);

        const auto advec = m_deriv.template advection<MatterCCZ4RHSVars>(
            current_cell, ccz4_vars.shift);

        data_t kappa1_times_lapse;
        if (m_covariantZ4)
            kappa1_times_lapse = m_kappa1;
        else
            kappa1_times_lapse = m_kappa1 * ccz4_vars.lapse;

        using namespace TensorAlgebra;

        auto h_UU = compute_inverse_sym(ccz4_vars.h);
        auto chris = compute_christoffel(d1.h, h_UU);

	const auto emtensor =
        m_matter.compute_emtensor(matter_vars, d1_matter, h_UU, chris.ULL);

	if (m_c_rho >= 0)
	{
	    current_cell.store_vars(emtensor.rho, m_c_rho);
	}

	if (m_c_Si.size() > 0)
    {   
#if DEFAULT_TENSOR_DIM == 3
        FOR(i) { current_cell.store_vars(emtensor.Si[i], m_c_Si.begin() + i); }                
#endif  
    }                          
                               
    if (m_c_Sij.size() > 0)    
    {                      
#if DEFAULT_TENSOR_DIM == 3
        current_cell.store_vars(emtensor.Sij[0][0], m_c_Sij.begin());
        current_cell.store_vars(emtensor.Sij[0][1], m_c_Sij.begin() + 1);                      
        current_cell.store_vars(emtensor.Sij[0][2], m_c_Sij.begin() + 2);                      
        current_cell.store_vars(emtensor.Sij[1][1], m_c_Sij.begin() + 3);                      
        current_cell.store_vars(emtensor.Sij[1][2], m_c_Sij.begin() + 4);                      
        current_cell.store_vars(emtensor.Sij[2][2], m_c_Sij.begin() + 5);                      
#endif  
    }   

	if (m_c_dtK >= 0)
	{
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
               (emtensor.S - 3 * emtensor.rho);

        // store the RHS of dtK as diagnostic variable
        current_cell.store_vars(dtK, c_dtK);
        }
    }
};

#endif /* DIAGNOSTICTIMEDERIVATIVEK_HPP_ */
