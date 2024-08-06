/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(BHBSBINARY_HPP_)
#error "This file should only be included through BHBSBinary.hpp"
#endif

#ifndef BHBSBINARY_IMPL_HPP_
#define BHBSBINARY_IMPL_HPP_

#include "BosonStarSolution.hpp" //for BosonStarSolution class
#include "DebuggingTools.hpp"
#include "Max.hpp"
#include "WeightFunction.hpp"
#ifdef USE_TWOPUNCTURES
#include "TwoPunctures.hpp" //for TwoPunctures-based ID method
#endif

inline BHBSBinary::BHBSBinary(BosonStar_params_t a_params_BosonStar,
                              BlackHole_params_t a_params_BlackHole,
                              Binary_params_t a_params_Binary,
                              Potential::params_t a_params_potential,
                              double a_G_Newton, double a_dx, int a_verbosity)
    : m_dx(a_dx), m_G_Newton(a_G_Newton),
      m_params_BosonStar(a_params_BosonStar),
      m_params_BlackHole(a_params_BlackHole),
      m_params_Binary(a_params_Binary),
      m_params_potential(a_params_potential), m_verbosity(a_verbosity)
{
}

#ifdef USE_TWOPUNCTURES
inline BHBSBinary::BHBSBinary(BosonStar_params_t a_params_BosonStar,
                              BlackHole_params_t a_params_BlackHole,
                              Binary_params_t a_params_Binary,
                              Potential::params_t a_params_potential,
                              double a_G_Newton, double a_dx, int a_verbosity,
                              TP::TwoPunctures *a_two_punctures)
    : m_dx(a_dx), m_G_Newton(a_G_Newton),
      m_params_BosonStar(a_params_BosonStar),
      m_params_BlackHole(a_params_BlackHole),
      m_params_Binary(a_params_Binary),
      m_params_potential(a_params_potential), m_verbosity(a_verbosity),
      m_two_punctures(a_two_punctures)
    {
    }
#endif

void BHBSBinary::compute_1d_BS_solution(const double max_r)
/** This function computes the 1d solution for the BS in the binary
 */
{
    try
    {
        // set initial parameters and then run the solver (didnt put it in the
        // constructor)

        pout() << "Setting initial conditions for Star 1" << endl;
        m_1d_sol.set_initialcondition_params(m_params_BosonStar,
                                             m_params_potential, max_r);
        int initial_data_construction = m_params_Binary.id_choice;
        pout() << "I am running initial data choice No: "
               << initial_data_construction << endl;
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

    /*
     *            INITIALIZATION OF THE PARAMETERS
     */

    //pout() << "Starting initial data calculation" << endl;

    // Load variables (should be set to zero if this is a single BS)
    MatterCCZ4<ComplexScalarField<>>::Vars<data_t> vars;
    current_cell.load_vars(vars);

    // Coordinates for centre of mass
    Coordinates<data_t> coords(current_cell, m_dx,
                               m_params_Binary.centre_of_mass);

    // Import BS parameters
    double BS_rapidity = m_params_BosonStar.rapidity;

    bool antiboson = m_params_BosonStar.antiboson;

    double BS_radius_width = m_params_BosonStar.radius_width;
    double R_BS = m_params_BosonStar.bump_radius;

    // pout() << "Boson star parameters:" << endl;
    // pout() << "BS mass = " << m_params_BosonStar.mass << endl;
    // pout() << "BS_rapidity = " << BS_rapidity << endl;
    // pout() << "R_BS = " << R_BS << endl;

    // Import BH parameters
    double BH_rapidity = m_params_BlackHole.rapidity;
    double M = m_params_BlackHole.mass;

    double BH_radius_width = m_params_BlackHole.radius_width;
    double R_BH = m_params_BlackHole.bump_radius;

    // pout() << "Black hole parameters:" << endl;
    // pout() << "BH mass = " << M << endl;
    // pout() << "BH_rapidity = " << BH_rapidity << endl;
    // pout() << "R_BH = " << R_BH << endl;

    // Import binary parameters. Note that the mass ratio is defined as q =
    // mBH/mBS
    int initial_data_choice = m_params_Binary.id_choice;
    double separation = m_params_Binary.separation;
    double impact_parameter = m_params_Binary.impact_parameter;
    double q = m_params_Binary.mass_ratio;

    int conformal_power = m_params_Binary.conformal_factor_power;
    double epsilon = m_params_Binary.epsilon;
    int weight_function_choice = m_params_Binary.weight_function_choice;

    // Initialise extrinsic curvature and metric with upper indices
    double KLL_1[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
    double KLL_2[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
    double KLL[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
    double gammaLL[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
    double gammaUU[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
    double gammaUU_1[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
    double gammaUU_2[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
    double K1 = 0., K2 = 0.;

    /*
     *            SET UP OF THE BOSON STAR
     */

    /* Define boosts and coordinate objects, suppose BS is on the right of
     * the centre of mass and BH is on the left of centre of mass,
     * i.e.taking the centre of mass to be the origin, then BH --------
     * (origin) -------- BS */

    // Boson star positioning. Recall tanh(rapidity) = v/c, and x' = Lambda x
    // with Lambda = ( c -s \\ -s c) 2x2 matrix. We determine (x,y,z,t) wrt to
    // the star centre
    double c_ = cosh(BS_rapidity);
    double s_ = sinh(BS_rapidity);
    double v_ = tanh(BS_rapidity);

    double t =
        (coords.x - q * separation / (q + 1.)) * s_; // set /tilde{t} to zero
    double x = (coords.x - q * separation / (q + 1.)) * c_;
    double y = coords.y + q * impact_parameter / (q + 1.);
    double z = coords.z; // set /tilde{t} to zero

    double r = sqrt(x * x + y * y + z * z);

    // Save relative coordinates to the star.
    // These are the coordinates w.r.t the new origin that is the center of the
    // star, in its own restframe.
    double x_star{x};
    double y_star{y};
    double z_star{z};
    double r_star{r};

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
    double w_ = m_1d_sol.get_w(); // frequency

    // Write in phase, shift, metric components of star 1 and initialise metric
    // components of star 2
    double phase_ = m_params_BosonStar.phase * M_PI + w_ * t;
    double beta_x = s_ * c_ * (psi_ * psi_ - omega_ * omega_) / (pc_os);
    vars.shift[0] += beta_x;
    double g_zz_1 = psi_ * psi_;
    double g_yy_1 = psi_ * psi_;
    double g_xx_1 = pc_os;

    // Add on to evolution equations
    vars.phi_Re += p_ * cos(phase_);
    vars.phi_Im += p_ * sin(phase_);
    vars.Pi_Re +=
        -(1. / lapse_1) * ((x / r) * (s_ - beta_x * c_) * dp_ * cos(phase_) -
                           w_ * (c_ - beta_x * s_) * p_ * sin(phase_));
    vars.Pi_Im +=
        -(1. / lapse_1) * ((x / r) * (s_ - beta_x * c_) * dp_ * sin(phase_) +
                           w_ * (c_ - beta_x * s_) * p_ * cos(phase_));

    // Metric upper indices
    gammaUU_1[0][0] = 1. / g_xx_1;
    gammaUU_1[1][1] = 1. / g_yy_1;
    gammaUU_1[2][2] = 1. / g_zz_1;

    // Extrinsic curvature boson star
    KLL_1[2][2] = -lapse_1 * s_ * x * psi_prime_ / (r * psi_);
    KLL_1[1][1] = KLL_1[2][2];
    KLL_1[0][1] = lapse_1 * c_ * s_ * (y / r) *
                  (psi_prime_ / psi_ - omega_prime_ / omega_);
    KLL_1[0][2] = lapse_1 * c_ * s_ * (z / r) *
                  (psi_prime_ / psi_ - omega_prime_ / omega_);
    KLL_1[1][0] = KLL_1[0][1];
    KLL_1[2][0] = KLL_1[0][2];
    KLL_1[0][0] = lapse_1 * (x / r) * s_ * c_ * c_ *
                  (psi_prime_ / psi_ - 2. * omega_prime_ / omega_ +
                   v_ * v_ * omega_ * omega_prime_ * pow(psi_, -2));
    FOR2(i, j) K1 += gammaUU_1[i][j] * KLL_1[i][j];

    /*  Here we use Thomas Helfer's trick and find the corresponding fixed
       values to be substracted in the initial guess Prepare the matrices that
       will store the corrections. Note that for equal mass helferLL = helferLL2
     */

    //! Helfer matrix to store influence of the BS on the BH
    double helferLL[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};

    //! Helfer matrix to store influence of the BH on the BS
    double helferLL2[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};

    /*! Here is the conformal factor that will be differently calculated
       depending on the choice of the initial data. Note that chi = 1 / lambda,
       with lambda what is used in Tamara's paper.*/
    double chi_;
    double chi_plain;

    // This is the effect of the BS (object 1) on the BH (object 2) and hence
    // represents the value to be substracted in the initial data from the
    // position of the BH In the second block, stuff gets calculated at the
    // position of the second star that got calculated before already.
    double t_p = (-separation) * s_; // set /tilde{t} to zero
    double x_p = (-separation) * c_;
    double z_p = 0.; // set /tilde{t} to zero
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

    // Uncomment if you want to output the metric at infinity
    // double chi_inf = pow((2. - helferLL[0][0]) * (2. - helferLL[1][1]) *
    // (2. - helferLL[2][2]), -1./3.);
    // double h00_inf = (2. - helferLL[0][0]) * chi_inf;
    // double h11_inf = (2. - helferLL[1][1]) * chi_inf;
    // double h22_inf = (2. - helferLL[2][2]) * chi_inf;
    /*if (r<3){
    std::cout << "h00 = " << h00_inf << ", h11 = " << h11_inf
                        << ", h22 = " << h22_inf << ", chi inf = " <<
                        chi_inf << std::endl;}*/

    /*
     *                  SET UP OF THE BLACK HOLE
     */

    // BH positioning
    c_ = cosh(-BH_rapidity);
    s_ = sinh(-BH_rapidity);
    v_ = tanh(-BH_rapidity);
    t = (coords.x + separation / (q + 1.)) * s_; // set /tilde{t} to zero
    x = (coords.x + separation / (q + 1.)) * c_;
    z = coords.z;
    y = coords.y - impact_parameter / (q + 1.);
    r = sqrt(x * x + y * y + z * z);

    // Save relative coordinates to the BH
    double x_BH{x};
    double y_BH{y};
    double z_BH{z};
    double r_BH{r};

    // BH contributions
    double r_tilde;
    r_tilde = sqrt(r * r + 10e-10);

    omega_ = (2. - M / r_tilde) / (2. + M / r_tilde);
    omega_prime_ = 4. * M / pow(2. * r_tilde + M, 2);
    psi_ = pow(1. + M / (2. * r_tilde), 2);
    psi_prime_ = -(M / (r_tilde * r_tilde)) * (1. + M / (2. * r_tilde));

    pc_os = psi_ * psi_ * c_ * c_ - omega_ * omega_ * s_ * s_;
    double lapse_2 = omega_ * psi_ / (sqrt(pc_os));

    if (antiboson)
    {
        w_ = -m_1d_sol.get_w();
    }
    else
    {
        w_ = m_1d_sol.get_w();
    }

    beta_x = s_ * c_ * (psi_ * psi_ - omega_ * omega_) / (pc_os);
    vars.shift[0] += beta_x;
    double g_zz_2 = psi_ * psi_;
    double g_yy_2 = psi_ * psi_;
    double g_xx_2 = pc_os;
    gammaUU_2[0][0] = 1. / g_xx_2;
    gammaUU_2[1][1] = 1. / g_yy_2;
    gammaUU_2[2][2] = 1. / g_zz_2;

    // do not need to modify scalar field or momentum if we have a black hole

    KLL_2[2][2] = -lapse_2 * s_ * x * psi_prime_ / (r * psi_);
    KLL_2[1][1] = KLL_2[2][2];
    KLL_2[0][1] = lapse_2 * c_ * s_ * (y / r) *
                  (psi_prime_ / psi_ - omega_prime_ / omega_);
    KLL_2[0][2] = lapse_2 * c_ * s_ * (z / r) *
                  (psi_prime_ / psi_ - omega_prime_ / omega_);
    KLL_2[1][0] = KLL_2[0][1];
    KLL_2[2][0] = KLL_2[0][2];
    KLL_2[2][1] = 0.;
    KLL_2[1][2] = 0.;
    KLL_2[0][0] = lapse_2 * (x / r) * s_ * c_ * c_ *
                  (psi_prime_ / psi_ - 2. * omega_prime_ / omega_ +
                   v_ * v_ * omega_ * omega_prime_ * pow(psi_, -2));
    FOR2(i, j) K2 += gammaUU_2[i][j] * KLL_2[i][j];

    // If one uses fixing conformal trick, we need to have the values of the
    // metric of the BS at its centre In the solution stored in m_1d_sol, this
    // is at the origin, as this is the single star solution
    double r_11 = 0.;
    double p_11 = m_1d_sol.get_p_interp(r_11);
    double dp_11 = m_1d_sol.get_dp_interp(r_11);
    double omega_11 = m_1d_sol.get_lapse_interp(r_11);
    double omega_prime_11 = m_1d_sol.get_dlapse_interp(r_11);
    double psi_11 = m_1d_sol.get_psi_interp(r_11);
    double psi_prime_11 = m_1d_sol.get_dpsi_interp(r_11);
    double pc_os_11 = psi_11 * psi_11 * cosh(BS_rapidity) * cosh(BS_rapidity) -
                      omega_11 * omega_11 * sinh(BS_rapidity) * sinh(BS_rapidity);

    // We need to calculate the influence of the BH at the BS centre.
    double x_p2 = separation * c_;
    double z_p2 = 0.; // set /tilde{t} to zero
    double y_p2 = -impact_parameter;
    double r_p2 = sqrt(x_p2 * x_p2 + y_p2 * y_p2 + z_p2 * z_p2 + 10e-10); // add 1e-10 similar to r_tilde above

    double omega_p2 = (2. - M / r_p2) / (2. + M / r_p2);
    double psi_p2 = pow(1. + M / (2. * r_p2), 2);
    double pc_os_p2 = psi_p2 * psi_p2 * c_ * c_ - omega_p2 * omega_p2 * s_ * s_;

    helferLL2[1][1] = psi_p2 * psi_p2;
    helferLL2[2][2] = psi_p2 * psi_p2;
    helferLL2[0][0] = pc_os_p2;

    /*
     *            SET UP OF THE BINARY INITIAL DATA
     */

    double g_xx, g_yy, g_zz;

    // Plain superposition
    if (initial_data_choice == 0)
    {
        g_xx = g_xx_1 + g_xx_2 - 1.0;
        g_yy = g_yy_1 + g_yy_2 - 1.0;
        g_zz = g_zz_1 + g_zz_2 - 1.0;

        // Define initial conformal factor
        chi_ = pow(g_xx * g_yy * g_zz, -1. / 3.);
    }

    // Thomas' trick
    if (initial_data_choice == 1)
    {
        g_xx = g_xx_1 + g_xx_2 - helferLL2[0][0];
        g_yy = g_yy_1 + g_yy_2 - helferLL2[1][1];
        g_zz = g_zz_1 + g_zz_2 - helferLL2[2][2];

        // Define initial conformal factor
        chi_ = pow(g_xx * g_yy * g_zz, -1. / 3.);
    }

    // Our new method of fixing the conformal factor (with arbitrary n power) to
    // its central equilibrium value
    if (initial_data_choice == 2)
    {
        // This is to be filled in with plain superposed metric components
        // evaluated at x^BS
        double superpose[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};

        // Start with plain superposed metrics
        double g_xx_plain = g_xx_1 + g_xx_2 - 1.0;
        double g_yy_plain = g_yy_1 + g_yy_2 - 1.0;
        double g_zz_plain = g_zz_1 + g_zz_2 - 1.0;

        // metric components of \gamma^BS(x^BS)
        double g_zz_11 = psi_11 * psi_11;
        double g_yy_11 = psi_11 * psi_11;
        double g_xx_11 = pc_os_11;

        // This  is \gamma_{ij}(x^BS) = \gamma^BS(x^BS) + \gamma^BH(x^BS) - 1
        superpose[0][0] = g_xx_11 + helferLL2[0][0] - 1.;
        superpose[1][1] = g_yy_11 + helferLL2[1][1] - 1.;
        superpose[2][2] = g_zz_11 + helferLL2[2][2] - 1.;

        double n_power = conformal_power / 12.0;

        // This is \chi(x^BS)
        double chi_1 =
            pow(superpose[0][0] * superpose[1][1] * superpose[2][2], n_power);

        // This is \chi^BS(x_BS)
        double chi1_1 = pow(g_xx_11 * g_yy_11 * g_zz_11, n_power);

        // This is \delta_BS
        double delta_1 = chi1_1 - chi_1;

        chi_plain = pow(g_xx_plain * g_yy_plain * g_zz_plain, n_power);

        // Create the weight function
        double profile1;

        if (weight_function_choice == 1)
        {
            WeightFunction weight(separation, x_star, y_star, z_star,
                                  m_params_Binary.weight_function_order);
            profile1 = weight.profile_chi();
        }
        else if (weight_function_choice == 2)
        {
            WeightFunctionAngle weight(separation, impact_parameter, x_star,
                                       y_star, z_star, epsilon, BS_radius_width);
            profile1 = weight.profile_chi();
        }

        // Adapted conformal factor
        chi_ = chi_plain + profile1 * delta_1;

        g_xx = chi_plain * g_xx_plain / chi_;
        g_yy = chi_plain * g_yy_plain / chi_;
        g_zz = chi_plain * g_zz_plain / chi_;
    }

    // weight functions on metric components
    if (initial_data_choice == 3)
    {
        // Create the weight function
        WeightFunction weight(separation, x_star, y_star, z_star,
                              m_params_Binary.weight_function_order);
        double profile_star =
            weight.profile_chi_2(x_star, y_star, z_star, BS_radius_width);
        double profile_hole =
            weight.profile_chi_2(x_BH, y_BH, z_BH, BH_radius_width);

        // Start with plain superposed metrics
        g_xx = g_xx_1 + g_xx_2 - 1.0 + profile_star * (1 - helferLL2[0][0]) +
               profile_hole * (1 - helferLL[0][0]);
        g_yy = g_yy_1 + g_yy_2 - 1.0 + profile_star * (1 - helferLL2[1][1]) +
               profile_hole * (1 - helferLL[1][1]);
        g_zz = g_zz_1 + g_zz_2 - 1.0 + profile_star * (1 - helferLL2[2][2]) +
               profile_hole * (1 - helferLL[2][2]);

        // Define initial conformal factor
        chi_ = pow(g_xx * g_yy * g_zz, -1. / 3.);
    }

    if (initial_data_choice == 4)
    {
        double g_xx_H = g_xx_1 + g_xx_2 - helferLL2[0][0];
        double g_yy_H = g_yy_1 + g_yy_2 - helferLL2[1][1];
        double g_zz_H = g_zz_1 + g_zz_2 - helferLL2[2][2];

        double chi_Helfer = pow(g_xx_H * g_yy_H * g_zz_H, -1. / 3.);

        double chi_inf = pow((2. - helferLL[0][0]) * (2. - helferLL[1][1]) *
                                 (2. - helferLL[2][2]),
                             -1. / 3.);

        double diff = 1-chi_inf;

       	// double profile1 = (tanh(pow(r_star / (separation * 5), 2.) - BS_radius_width) + tanh(BS_radius_width)) / (1+tanh(BS_radius_width));
	double profile1 = 1/BS_radius_width - 1/sqrt(BS_radius_width*BS_radius_width + r_star * r_star); 

        chi_ = chi_Helfer + profile1 * diff;

        g_xx = chi_Helfer * g_xx_H / chi_;
        g_yy = chi_Helfer * g_yy_H / chi_;
        g_zz = chi_Helfer * g_zz_H / chi_;
    }

    if (initial_data_choice != 6) 
    {	    
    // Now, compute upper and lower components
    gammaLL[0][0] = g_xx;
    gammaLL[1][1] = g_yy;
    gammaLL[2][2] = g_zz;
    gammaUU[0][0] = 1. / g_xx;
    gammaUU[1][1] = 1. / g_yy;
    gammaUU[2][2] = 1. / g_zz;

    vars.chi = chi_;

    // Define initial lapse
    vars.lapse += sqrt(vars.chi);

    // Define initial trace of K and A_ij
    double one_third = 1. / 3.;
    FOR2(i, j) vars.h[i][j] = chi_ * gammaLL[i][j];
    FOR4(i, j, k, l)
    KLL[i][j] += gammaLL[i][l] * (gammaUU_1[l][k] * KLL_1[k][j] +
                                  gammaUU_2[l][k] * KLL_2[k][j]);
    FOR2(i, j) vars.K += KLL[i][j] * gammaUU[i][j];
    FOR2(i, j)
    vars.A[i][j] = chi_ * (KLL[i][j] - one_third * vars.K * gammaLL[i][j]);

    current_cell.store_vars(vars);
    }

#ifdef USE_TWOPUNCTURES
     //TwoPunctures based approach
    if (initial_data_choice == 6)
    {

	    //start with data corrected according to standard Thomas' trick as in id_choice = 1
	    g_xx = g_xx_1 + g_xx_2 - helferLL2[0][0];
        g_yy = g_yy_1 + g_yy_2 - helferLL2[1][1];
        g_zz = g_zz_1 + g_zz_2 - helferLL2[2][2];

        //Now, compute upper and lower components
        gammaLL[0][0] = g_xx;
        gammaLL[1][1] = g_yy;
        gammaLL[2][2] = g_zz;
        gammaUU[0][0] = 1. / g_xx;
        gammaUU[1][1] = 1. / g_yy;
        gammaUU[2][2] = 1. / g_zz;

        // Define initial conformal factor
        double  chiThomas = pow(g_xx * g_yy * g_zz, -1. / 3.);

        //physical metric and ext. curvature for thomas' trick
	    double gammaThomas[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};
	    double KLLThomas[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};	
		
	    FOR2(i,j) gammaThomas[i][j] = gammaLL[i][j];
        FOR4(i,j,k,l) KLLThomas[i][j] += gammaLL[i][l] * (gammaUU_1[l][k] * KLL_1[k][j] + gammaUU_2[l][k] * KLL_2[k][j]);

	    //TwoPunctures and final versions of CCZ4 vars
	    Tensor<2, double> gamma_TP, K_TP, KLLFinal, gammaLLFinal, gammaUUFinal;
    	Tensor<1, double> shiftTP, Z3TP;
    	double lapseTP, ThetaTP;

        //radial distance to BH and weight function value
        double r_hole = sqrt(pow(x_BH,2) + pow(y_BH,2) + pow (z_BH,2));

    	double TPFactor = 1 - tanh(r_hole * r_hole /( R_BH * R_BH));

	    //ensure TPFactor dies off entirely around BS center (where TP solution diverges)
	    double r_star = sqrt(pow(x_star,2) + pow(y_star,2) + pow(z_star,2) );
	    if (r_star < R_BS)
	        TPFactor = 0.;

	    // rotated coordinates to extract data from BHBSAMR
	    double coords_array[CH_SPACEDIM];
    	coords_array[0] = coords.x;
    	coords_array[1] = coords.y;
    	coords_array[2] = coords.z;

        //amount the original configuration has been  rotated relative to x-axis parallel config. assumes impact param/ separation > 0
        double rotation_angle = asin(impact_parameter / sqrt(separation * separation + impact_parameter * impact_parameter ));

        //rotate by rotation_angle to obtain coordinates of source point in TP system
        if (m_params_Binary.do_rotation)
        {
            coords_array[0] = cos(rotation_angle) * coords.x - sin(rotation_angle) * coords.y;
            coords_array[1] = sin(rotation_angle) * coords.x + cos(rotation_angle) * coords.y;
        }

	    using namespace TP::Z4VectorShortcuts;
        double TP_state[Qlen];	

        // pout() << "masses: " << (*m_two_punctures).mp << " " << (*m_two_punctures).mp << " " << (*m_two_punctures).mm_adm << " " << (*m_two_punctures).mp_adm << endl;
	    
        //Read TwoPunctures data from the TPAMR initialized in Main
        (*m_two_punctures).Interpolate(coords_array, TP_state);	
	
        // pout() << "gamma_TP[0][0]" << gamma_TP[0][0] << ", TP_state[g11] " << TP_state[g11] << endl;

	    // TP metric
        gamma_TP[0][0] = TP_state[g11];
        gamma_TP[0][1] = gamma_TP[1][0] = TP_state[g12];
        gamma_TP[0][2] = gamma_TP[2][0] = TP_state[g13];
        gamma_TP[1][1] = TP_state[g22];
        gamma_TP[1][2] = gamma_TP[2][1] = TP_state[g23];
        gamma_TP[2][2] = TP_state[g33];

        // TP extrinsic curvature
        K_TP[0][0] = TP_state[K11];
        K_TP[0][1] = K_TP[1][0] = TP_state[K12];
        K_TP[0][2] = K_TP[2][0] = TP_state[K13];
        K_TP[1][1] = TP_state[K22];
        K_TP[1][2] = K_TP[2][1] = TP_state[K23];
        K_TP[2][2] = TP_state[K33];
	
	    lapseTP = TP_state[lapse];	

        //apply TP correction to physical metric and extrinsic curvature
        FOR2(i,j) gammaLLFinal[i][j] = gammaThomas[i][j] + TPFactor * (gamma_TP[i][j] - gammaThomas[i][j]);
        //FOR2(i,j) KLL[i][j] = KLLThomas[i][j] + TPFactor * (K_TP[i][j] - KLLThomas[i][j]);
        FOR2(i,j) KLL[i][j] = KLLThomas[i][j]; // seems to work better without K fix
        
        //vars.lapse = vars.lapse + TPFactor * (lapseTP - vars.lapse); 

        //get corrected chi and inverse metric
        vars.chi = pow(TensorAlgebra::compute_determinant_sym(gammaLLFinal), -1.0 / 3.0);
        gammaUUFinal = TensorAlgebra::compute_inverse_sym(gammaLLFinal);
            
        //initial lapse
        vars.lapse += sqrt(vars.chi);

        //finally reconstruct h and A vars
        double one_third = 1./3.;
        FOR2(i,j) vars.h[i][j] = vars.chi * gammaLLFinal[i][j];
        FOR2(i,j) vars.K += KLL[i][j] * gammaUUFinal[i][j];
            FOR2(i,j) vars.A[i][j] = vars.chi * (KLL[i][j] - one_third * vars.K * gammaLLFinal[i][j]);

            current_cell.store_vars(vars);
	}
#endif
}

#endif /* BHBS_IMPL_HPP_ */
