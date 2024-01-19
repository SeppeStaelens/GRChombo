
/* GRChombo
* Copyright 2012 The GRChombo collaboration.
* Please refer to LICENSE in GRChombo's root directory.
*/

#if !defined(DENSITY_HPP_)
#error "This file should only be included through Density.hpp"
#endif

#ifndef DENSITY_IMPL_HPP_
#define DENSITY_IMPL_HPP_
#include "DimensionDefinitions.hpp"

template <class matter_t>
Density<matter_t>::Density(const matter_t a_matter,
                                              double dx, double G_Newton)
   : Constraints(dx, 0.0 /*No cosmological constant*/), my_matter(a_matter),
     m_G_Newton(G_Newton)
{
}

template <class matter_t>
template <class data_t>
void Density<matter_t>::compute(Cell<data_t> current_cell) const
{
   // Load local vars and calculate derivs
   const auto vars = current_cell.template load_vars<Vars>();
   const auto d1 = m_deriv.template diff1<Vars>(current_cell);
   const auto d2 = m_deriv.template diff2<Vars>(current_cell);

   // Inverse metric and Christoffel symbol
   const auto h_UU = TensorAlgebra::compute_inverse_sym(vars.h);
   const auto chris = TensorAlgebra::compute_christoffel(d1.h, h_UU);

   // Energy Momentum Tensor
   const auto emtensor = my_matter.compute_emtensor(vars, d1, h_UU, chris.ULL);

   // Write the rhs into the output FArrayBox
   current_cell.store_vars(emtensor.rho, c_rho);
   current_cell.store_vars(emtensor.Si, GRInterval<c_s1, c_s3>());
}

#endif /* DENSITY_IMPL_HPP_ */
