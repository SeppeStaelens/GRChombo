/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef NOETHERCHARGE_HPP_
#define NOETHERCHARGE_HPP_

#include "ADMConformalVars.hpp" // needed for CCz4 and matter variables
#include "Cell.hpp"
#include "ComplexScalarField.hpp"
#include "Coordinates.hpp"
#include "DimensionDefinitions.hpp"
#include "FourthOrderDerivatives.hpp"
#include "Interval.H"
#include "UserVariables.hpp"
#include "simd.hpp"

//! Calculates the Noether Charge integrand values and the modulus of the
//! complex scalar field on the grid
class NoetherCharge
{
    // Need matter variables and chi
    template <class data_t>
    using ADMVars = ADMConformalVars::VarsNoGauge<data_t>;
    template <class data_t>
    using MatterVars = ComplexScalarField<>::Vars<data_t>;

    template <class data_t> using RHSVars = MatterCCZ4RHS<>::Vars<data_t>;

    //   protected:
    //     FourthOrderDerivatives m_deriv;
    //     const int m_c_rho;      // var enum for the energy density
    //     const Interval m_c_Si;  // Interval of var enums for the momentum
    //     density const Interval m_c_Sij; // Interval of var enums for the
    //     spatial
    // stress-energy density

  public:
    // //! Constructor
    // NoetherCharge(const int a_c_rho = -1, const Interval a_c_Si = Interval(),
    //               const Interval a_c_Sij = Interval())
    //     : m_c_rho(a_c_rho), m_c_Si(a_c_Si), m_c_Sij(a_c_Sij)
    // {
    //     if (m_c_Si.size() != 0)
    //     {
    //         // Si is a vector
    //         CH_assert(m_c_Si.size() == DEFAULT_TENSOR_DIM);
    //     }

    //     if (m_c_Sij.size() != 0)
    //     {
    //         // Sij is a symmetric tensor
    //         CH_assert(m_c_Sij.size() ==
    //                   DEFAULT_TENSOR_DIM * (DEFAULT_TENSOR_DIM + 1) / 2);
    //     }
    // }

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // load vars locally
        const auto adm_vars = current_cell.template load_vars<ADMVars>();
        const auto matter_vars = current_cell.template load_vars<MatterVars>();
        const auto matter_rhs = current_cell.template load_vars<RHSVars>();

        // calculate Noether charge
        data_t N =
            pow(adm_vars.chi, -1.5) * (matter_vars.phi_Im * matter_vars.Pi_Re -
                                       matter_vars.phi_Re * matter_vars.Pi_Im);

        data_t mod_phi = sqrt(matter_vars.phi_Re * matter_vars.phi_Re +
                              matter_vars.phi_Im * matter_vars.phi_Im);

        current_cell.store_vars(N, c_N);
        current_cell.store_vars(mod_phi, c_mod_phi);

        /////// ADDITIONS

        // store the RHS of dtK as diagnostic variable
        current_cell.store_vars(matter_rhs.K, c_dtK);

        //         /////// EMTENSOR ADDITIONS

        //         const auto d1 = m_deriv.template diff1<Vars>(current_cell);

        //         using namespace TensorAlgebra;

        //         const auto h_UU = compute_inverse_sym(vars.h);
        //         const auto chris = compute_christoffel(d1.h, h_UU);

        //         const auto emtensor =
        //             m_matter.compute_emtensor(matter_vars, d1, h_UU,
        //             chris.ULL);

        //         if (m_c_rho >= 0)
        //         {
        //             current_cell.store_vars(emtensor.rho, m_c_rho);
        //         }

        //         if (m_c_Si.size() > 0)
        //         {
        // #if DEFAULT_TENSOR_DIM == 3
        //             FOR1(i)
        //             {
        //                 current_cell.store_vars(emtensor.Si[i],
        //                 m_c_Si.begin() + i);
        //             }
        // #endif
        //         }

        //         if (m_c_Sij.size() > 0)
        //         {
        // #if DEFAULT_TENSOR_DIM == 3
        //             current_cell.store_vars(emtensor.Sij[0][0],
        //             m_c_Sij.begin());
        //             current_cell.store_vars(emtensor.Sij[0][1],
        //             m_c_Sij.begin() + 1);
        //             current_cell.store_vars(emtensor.Sij[0][2],
        //             m_c_Sij.begin() + 2);
        //             current_cell.store_vars(emtensor.Sij[1][1],
        //             m_c_Sij.begin() + 3);
        //             current_cell.store_vars(emtensor.Sij[1][2],
        //             m_c_Sij.begin() + 4);
        //             current_cell.store_vars(emtensor.Sij[2][2],
        //             m_c_Sij.begin() + 5);
        // #endif
        //         }
    }
};

#endif /* NOETHERCHARGE_HPP_ */
