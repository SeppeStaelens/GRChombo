/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(BHBSBINARY_HPP_)
#error "This file should only be included through BosonStar.hpp"
#endif

#ifndef BHBSBINARY_IMPL_HPP_
#define BHBSBINARY_IMPL_HPP_

#include "BosonStarSolution.hpp" //for BosonStarSolution class
#include "WeightFunction.hpp"
#include "DebuggingTools.hpp"
#include "Max.hpp"

inline BHBSBinary::BHBSBinary(BosonStar_params_t a_params_BosonStar, BlackHole_params_t a_params_BlackHole,
                    Potential::params_t a_params_potential, double a_G_Newton,
                    double a_dx, int a_verbosity)
    :m_dx(a_dx), m_G_Newton(a_G_Newton), m_params_BosonStar(a_params_BosonStar), m_params_BlackHole(a_params_BlackHole),
    m_params_potential(a_params_potential), m_verbosity(a_verbosity)
{
}

void BHBSBinary::compute_1d_solution(const double max_r)
/** This function computes the 1d solution for both BSs in the binary
*/
{
    try
    {   
        //set initial parameters and then run the solver (didnt put it in the constructor)
        
        pout() << "Setting initial conditions for Star 1" << endl;
        m_1d_sol.set_initialcondition_params(m_params_BosonStar,m_params_potential,max_r);
        int initial_data_construction = m_params_BosonStar.id_choice;
        pout() << "I am running initial data choice No: " << initial_data_construction << endl;
        pout() << "Running the solver for Star 1" << endl;
        m_1d_sol.main();
        pout() << "Completed for star 1" << endl;
    }
    catch (std::exception &exception)
    {
        pout() << exception.what() << "\n";
    }
}

// Compute the value of the initial vars on the grid
template <class data_t>
void BHBSBinary::compute(Cell<data_t> current_cell) const
{
    MatterCCZ4<ComplexScalarField<>>::Vars<data_t> vars;
    // Load variables (should be set to zero if this is a single BS)
    
    current_cell.load_vars(vars);
    //VarsTools::assign(vars, 0.); // Set only the non-zero components below
    
    // Coordinates for centre of mass
    Coordinates<data_t> coords(current_cell, m_dx,
        m_params_BosonStar.star_centre);

    // Import BS, BH parameters
    double rapidity = m_params_BosonStar.BS_rapidity;
    double rapidity2 = m_params_BlackHole.BH_rapidity;
    bool antiboson = m_params_BosonStar.antiboson;
    double M = m_params_BlackHole.BlackHoleMass;
    double separation = m_params_BosonStar.binary_separation;
    double impact_parameter = m_params_BosonStar.BS_impact_parameter;
    double q = m_params_BosonStar.mass_ratio;
    double radius_width1 = m_params_BosonStar.radius_width1;
    double radius_width2 = m_params_BosonStar.radius_width2;
    int conformal_power = m_params_BosonStar.conformal_factor_power;
    int initial_data_choice = m_params_BosonStar.id_choice;

    // Define boosts and coordinate objects, suppose star 1 is on the left of the centre of mass 
    // and star 2 is on the right of centre of mass

    //e.g. taking the centre of mass to be the origin, then STAR2 -------- (origin) -------- STAR1
    
    // First star positioning
    double c_ = cosh(rapidity);
    double s_ = sinh(rapidity);
    double v_ = tanh(rapidity);
    double t = (coords.x - q * separation / (q + 1.)) * s_; //set /tilde{t} to zero
    double x = (coords.x - q * separation / (q + 1.)) * c_;
    double z = coords.z; //set /tilde{t} to zero
    double y = coords.y + q * impact_parameter / (q + 1.);
    double r = sqrt(x * x + y * y + z * z);

    // Save relative coordinates to the star
    double x_star {x};
    double y_star {y};
    double z_star {z};
    double r_star {r};

    // First star physical variables
    double p_ = m_1d_sol.get_p_interp(r);
    double dp_ = m_1d_sol.get_dp_interp(r);
    double omega_ = m_1d_sol.get_lapse_interp(r);
    double omega_prime_ = m_1d_sol.get_dlapse_interp(r);
    double psi_ = m_1d_sol.get_psi_interp(r);
    double psi_prime_ = m_1d_sol.get_dpsi_interp(r);

    // Get scalar field modulus, conformal factor, lapse and their gradients
    double pc_os = psi_ * psi_ * c_ * c_ - omega_ * omega_ * s_ * s_;
    double lapse_1 = omega_ * psi_ / (sqrt(pc_os));
    double lapse_2 = 1.;
    double w_ = m_1d_sol.get_w();

    //Write in phase, shift, metric components of star 1 and initialise metric components of star 2
    double phase_ = m_params_BosonStar.phase * M_PI + w_ * t;
    double beta_x = s_ * c_ * (psi_ * psi_ - omega_ * omega_) / (pc_os);
    vars.shift[0] += beta_x;
    double g_zz_1 = psi_ * psi_;
    double g_yy_1 = psi_ * psi_;
    double g_xx_1 = pc_os;

    double g_xx_2 = 0., g_yy_2 = 0., g_zz_2 = 0., g_xx, g_yy, g_zz;

    //Add on to evolution equations
    vars.phi_Re += p_ * cos(phase_);
    vars.phi_Im += p_ * sin(phase_);
    vars.Pi_Re += -(1. / lapse_1) * ((x / r) * (s_ - beta_x * c_) * dp_ * cos(phase_) - w_ * (c_ - beta_x * s_) * p_ * sin(phase_));
    vars.Pi_Im += -(1. / lapse_1) * ((x / r) * (s_ - beta_x * c_) * dp_ * sin(phase_) + w_ * (c_ - beta_x * s_) * p_ * cos(phase_));

    //Initialise extrinsic curvature and metric with upper indices
    double KLL_1[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};
    double KLL_2[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};
    double KLL[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};
    double gammaLL[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};
    double gammaUU[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};
    double gammaUU_1[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};
    double gammaUU_2[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};
    double K1 = 0., K2 = 0.;

    // Fill them in
    gammaUU_1[0][0] = 1. / g_xx_1;
    gammaUU_1[1][1] = 1. / g_yy_1;
    gammaUU_1[2][2] = 1. / g_zz_1;

    KLL_1[2][2] = -lapse_1 * s_ * x * psi_prime_ / (r * psi_);
    KLL_1[1][1] = KLL_1[2][2];
    KLL_1[0][1] = lapse_1 * c_ * s_ * (y / r) * (psi_prime_ / psi_ - omega_prime_ / omega_);
    KLL_1[0][2] = lapse_1 * c_ * s_ * (z / r) * (psi_prime_ / psi_ - omega_prime_ / omega_);
    KLL_1[1][0] = KLL_1[0][1];
    KLL_1[2][0] = KLL_1[0][2];
    KLL_1[2][1] = 0.;
    KLL_1[1][2] = 0.;
    KLL_1[0][0] = lapse_1 * (x / r) * s_ * c_ * c_ * (psi_prime_ / psi_ - 2. * omega_prime_ / omega_ + v_ * v_ * omega_ * omega_prime_ * pow(psi_, -2));
    FOR2(i,j) K1 += gammaUU_1[i][j] * KLL_1[i][j];

    // Here we use Thomas Helfer's trick and find the corresponding fixed values to be substracted in the initial guess
    double helferLL[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};
    double helferLL2[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};
    // Note that for equal mass helferLL = helferLL2
    
    // Here is the conformal factor that will be differently calculated depending on the choice of the initial data
    double chi_;
    double chi_plain;

    // This is the effect of object 1 on object 2 and hence represents the value to be substracted in the initial data from the position of object 2
    // In the second block, stuff gets calculated at the position of the second star that got calculated before already. 
    double t_p = (-separation) * s_; //set /tilde{t} to zero
    double x_p = (-separation) * c_;
    double z_p = 0.; //set /tilde{t} to zero
    double y_p = impact_parameter;
    double r_p = sqrt(x_p * x_p + y_p * y_p + z_p * z_p);

    double p_p = m_1d_sol.get_p_interp(r_p);
    double dp_p = m_1d_sol.get_dp_interp(r_p);
    double omega_p = m_1d_sol.get_lapse_interp(r_p);
    double omega_prime_p = m_1d_sol.get_dlapse_interp(r_p);
    double psi_p = m_1d_sol.get_psi_interp(r_p);
    double psi_prime_p = m_1d_sol.get_dpsi_interp(r_p);

    double pc_os_p = psi_p * psi_p * c_ * c_ - omega_p * omega_p * s_ * s_;
           
    // compare this to g_ll_1 above
    helferLL[1][1] = psi_p * psi_p;
    helferLL[2][2] = psi_p * psi_p;
    helferLL[0][0] = pc_os_p;

    double chi_inf = pow((2. - helferLL[0][0]) * (2. - helferLL[1][1]) *
    (2. - helferLL[2][2]), -1./3.);
    double h00_inf = (2. - helferLL[0][0]) * chi_inf; 
    double h11_inf = (2. - helferLL[1][1]) * chi_inf;
    double h22_inf = (2. - helferLL[2][2]) * chi_inf;
    /*if (r<3){
    std::cout << "h00 = " << h00_inf << ", h11 = " << h11_inf
                        << ", h22 = " << h22_inf << ", chi inf = " <<
                        chi_inf << std::endl;}*/


    // BH positioning
    c_ = cosh(-rapidity2);
    s_ = sinh(-rapidity2);
    v_ = tanh(-rapidity2);
    t = (coords.x + separation / (q + 1.)) * s_; //set /tilde{t} to zero
    x = (coords.x + separation / (q + 1.)) * c_;
    z = coords.z;
    y = coords.y - impact_parameter / (q + 1.);
    r = sqrt(x * x + y * y + z * z);

    // Save relative coordinates to the BH
    double x_BH {x};
    double y_BH {y};
    double z_BH {z};
    double r_BH {r};

    // BH contributions
    double r_tilde;
    r_tilde = sqrt(r  *r + 10e-10);

    omega_ = (2. - M / r_tilde) / (2. + M / r_tilde);
    omega_prime_ = 4. * M / pow(2. * r_tilde + M, 2);
    psi_ = pow(1. + M/ (2. * r_tilde), 2);
    psi_prime_ = -(M / (r_tilde * r_tilde)) * (1. + M / (2. * r_tilde));

    pc_os = psi_ * psi_ * c_ *c_ - omega_ * omega_ * s_ * s_;
    lapse_2 = omega_ * psi_ / (sqrt(pc_os));

    if (antiboson)
    {
        w_ = - m_1d_sol.get_w();
    }
    else
    {
        w_ = m_1d_sol.get_w();
    }

    phase_ = w_ * t;
    beta_x = s_ * c_ * (psi_ * psi_ - omega_ * omega_) / (pc_os);
    vars.shift[0] += beta_x;
    g_zz_2 = psi_ * psi_;
    g_yy_2 = psi_ * psi_;
    g_xx_2 = pc_os;
    gammaUU_2[0][0] = 1. / g_xx_2;
    gammaUU_2[1][1] = 1. / g_yy_2;
    gammaUU_2[2][2] = 1. / g_zz_2;

    //do not need to modify scalar field or momentum if we have a black hole
   
    KLL_2[2][2] = -lapse_2 * s_ * x * psi_prime_ / (r * psi_);
    KLL_2[1][1] = KLL_2[2][2];
    KLL_2[0][1] = lapse_2 * c_ * s_ * (y / r) * (psi_prime_ / psi_ - omega_prime_ / omega_);
    KLL_2[0][2] = lapse_2 * c_ * s_ * (z / r) * (psi_prime_ / psi_ - omega_prime_ / omega_);
    KLL_2[1][0] = KLL_2[0][1];
    KLL_2[2][0] = KLL_2[0][2];
    KLL_2[2][1] = 0.;
    KLL_2[1][2] = 0.;
    KLL_2[0][0] = lapse_2 * (x / r) * s_ * c_ * c_ * (psi_prime_ / psi_ - 2. * omega_prime_ / omega_ + v_ * v_ * omega_ * omega_prime_ * pow(psi_, -2));
    FOR2(i,j) K2 += gammaUU_2[i][j] * KLL_2[i][j];

    //Plain superposition 
    if (initial_data_choice == 0)
    {
        g_xx = g_xx_1 + g_xx_2 - 1.0;
        g_yy = g_yy_1 + g_yy_2 - 1.0;
        g_zz = g_zz_1 + g_zz_2 - 1.0;

        //Now, compute upper and lower components
        gammaLL[0][0] = g_xx;
        gammaLL[1][1] = g_yy;
        gammaLL[2][2] = g_zz;
        gammaUU[0][0] = 1. / g_xx;
        gammaUU[1][1] = 1. / g_yy;
        gammaUU[2][2] = 1. / g_zz;
    
        double test_chi = g_xx * g_yy * g_zz;
    
        // Define initial conformal factor
        chi_ = pow(g_xx * g_yy * g_zz, -1. / 3.);

        vars.chi = chi_;

        // Define initial lapse
        vars.lapse += sqrt(vars.chi);

        // Define initial trace of K and A_ij
        double one_third = 1./3.;
        FOR2(i,j) vars.h[i][j] = vars.chi * gammaLL[i][j];
        FOR4(i,j,k,l) KLL[i][j] += gammaLL[i][l] * (gammaUU_1[l][k] * KLL_1[k][j] + gammaUU_2[l][k] * KLL_2[k][j]);
        FOR2(i,j) vars.K += KLL[i][j] * gammaUU[i][j];
        FOR2(i,j) vars.A[i][j] = vars.chi * (KLL[i][j] - one_third * vars.K * gammaLL[i][j]);

        current_cell.store_vars(vars);
    }
    
    //Thomas' trick
    if (initial_data_choice == 1)
    {
        g_xx = g_xx_1 + g_xx_2 - helferLL[0][0];
        g_yy = g_yy_1 + g_yy_2 - helferLL[1][1];
        g_zz = g_zz_1 + g_zz_2 - helferLL[2][2];

        //Now, compute upper and lower components
        gammaLL[0][0] = g_xx;
        gammaLL[1][1] = g_yy;
        gammaLL[2][2] = g_zz;
        gammaUU[0][0] = 1. / g_xx;
        gammaUU[1][1] = 1. / g_yy;
        gammaUU[2][2] = 1. / g_zz;
        
        // Define initial conformal factor
        chi_ = pow(g_xx * g_yy * g_zz, -1. / 3.);

        vars.chi = chi_;

        // Define initial lapse
        vars.lapse += sqrt(vars.chi);

        // Define initial trace of K and A_ij
        double one_third = 1./3.;
        FOR2(i,j) vars.h[i][j] = vars.chi * gammaLL[i][j];
        FOR4(i,j,k,l) KLL[i][j] += gammaLL[i][l] * (gammaUU_1[l][k] * KLL_1[k][j] + gammaUU_2[l][k] * KLL_2[k][j]);
        FOR2(i,j) vars.K += KLL[i][j] * gammaUU[i][j];
        FOR2(i,j) vars.A[i][j] = vars.chi * (KLL[i][j] - one_third * vars.K * gammaLL[i][j]);

        current_cell.store_vars(vars);
    }
    
    //Our new method of fixing the conformal factor (with arbitrary n power) to its central equilibrium value
    if (initial_data_choice == 2)
    {
        //This is to be filled in with plain superposed metric components evaluated at x_BS
        double superpose[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};
        
        //If one uses fixing conformal trick, we need to have the vales of the metric of star 1 at its centre
        //In the solution stored in m_1d_sol, this is at the origin, as this is the single star solution
        double r_11 = 0.;
        double p_11 = m_1d_sol.get_p_interp(r_11);
        double dp_11 = m_1d_sol.get_dp_interp(r_11);
        double omega_11 = m_1d_sol.get_lapse_interp(r_11);
        double omega_prime_11 = m_1d_sol.get_dlapse_interp(r_11);
        double psi_11 = m_1d_sol.get_psi_interp(r_11);
        double psi_prime_11 = m_1d_sol.get_dpsi_interp(r_11);
        double pc_os_11 = psi_11 * psi_11 * cosh(rapidity) * cosh(rapidity) - omega_11 * omega_11 * sinh(rapidity) * sinh(rapidity);

        // We need to calculate the influence of the BH at the BS centre.         
        double x_p2 = separation * c_;
        double z_p2 = 0.; //set /tilde{t} to zero
        double y_p2 = -impact_parameter;
        double r_p2 = sqrt(x_p2 * x_p2 + y_p2 * y_p2 + z_p2 * z_p2);
               
        double omega_p2 = (2. - M / r_p2) / (2. + M / r_p2);
        double psi_p2 = pow(1. + M/ (2. * r_p2), 2);
        double pc_os_p2 = psi_p2 * psi_p2 * c_ *c_ - omega_p2 * omega_p2 * s_ * s_;

        helferLL2[1][1] = psi_p2 * psi_p2;
        helferLL2[2][2] = psi_p2 * psi_p2;
        helferLL2[0][0] = pc_os_p2;               

        //Start with plain superposed metrics
        g_xx = g_xx_1 + g_xx_2 - 1.0;
        g_yy = g_yy_1 + g_yy_2 - 1.0;
        g_zz = g_zz_1 + g_zz_2 - 1.0;

        //metric components of \gamma_A(x_A)
        double g_zz_11 = psi_11 * psi_11;
        double g_yy_11 = psi_11 * psi_11;
        double g_xx_11 = pc_os_11;

        // This  is \gamma_{ij}(x_A) = \gamma_A(x_A) + \gamma_B(x_A) - 1
        superpose[0][0] = g_xx_11 + helferLL2[0][0] - 1.;
        superpose[1][1] = g_yy_11 + helferLL2[1][1] - 1.;
        superpose[2][2] = g_zz_11 + helferLL2[2][2] - 1.;

        double n_power = conformal_power / 12.0;

        //This is \chi(x_BS)
        double chi_1 = pow(superpose[0][0] * superpose[1][1] * superpose[2][2], n_power);

        //This is \chi^BS(x_BS)
        double chi1_1 = pow(g_xx_11 * g_yy_11 * g_zz_11, n_power);

        //This is \delta_BS
        double delta_1 = chi1_1 - chi_1;

        chi_plain = pow(g_xx * g_yy * g_zz, n_power);

        // Create the weight function
        WeightFunction weight(separation, x_star, y_star, z_star, m_params_BlackHole.weight_function_order);
        double profile1 = weight.profile_chi((coords.x - q * separation / (q+1)) * cosh(rapidity), coords.y + q * impact_parameter / (q + 1.), coords.z);
       
        // Adapted conformal factor
        chi_ = chi_plain + profile1 * delta_1;

        // Now, compute upper and lower components
        gammaLL[0][0] = g_xx;
        gammaLL[1][1] = g_yy;
        gammaLL[2][2] = g_zz;
        gammaUU[0][0] = 1. / g_xx;
        gammaUU[1][1] = 1. / g_yy;
        gammaUU[2][2] = 1. / g_zz;

        vars.chi = pow(chi_, - 4.0 / conformal_power);

        // Define initial lapse
        vars.lapse += sqrt(vars.chi);

        // Define initial trace of K and A_ij
        double one_third = 1./3.;
        FOR2(i,j) vars.h[i][j] = pow(chi_plain, - 4.0 / conformal_power) * gammaLL[i][j];
        FOR4(i,j,k,l) KLL[i][j] += gammaLL[i][l] * (gammaUU_1[l][k] * KLL_1[k][j] + gammaUU_2[l][k] * KLL_2[k][j]);
        FOR2(i,j) vars.K += KLL[i][j] * gammaUU[i][j];
        FOR2(i,j) vars.A[i][j] = pow(chi_plain, - 4.0 / conformal_power)  * (KLL[i][j] - one_third * vars.K * gammaLL[i][j]);

        current_cell.store_vars(vars);
    }
}

#endif /* BHBS_IMPL_HPP_ */