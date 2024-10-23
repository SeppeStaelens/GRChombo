
#ifndef BOSONSTARSOLUTION_HPP_
#define BOSONSTARSOLUTION_HPP_

#include "cmath"
#include <vector>

class BosonStarSolution
{

  private:         // private member variables/arrays
    double MM, A0; // Klein Gordon mass squared, KG scalr field central
                   // amplitude
    double PSC, OMC;
    // double PSC=2.52749283231717, OMC=0.160056080542207; // central density of
    // scalar field (0.272 for kaup)  PSC and OMC are central values of
    // conformal factor and lapse, not important as long as they are sensible
    // (i.e. order 1) PSC=2.0, OMC=0.5
    double lambda;  // phi 4 coupling in Klein gordon potential
    double sigma;   // 0.2 works with PC = 0.05 // parameter for solitonic stars
    bool solitonic; // false fro mini/lambda star. true for solitonic star
    bool BS_verbosity; // outputs more messages whilst finding the solution
    double EIGEN = 0;  // the desired eigenstate, 0 for ground
    int gridsize, adaptive_buffer;                // anywhere from 2k-200k is ok
    const int adaptive_stepsize_repetitions = 20; // 50; // 0 for no adaptive
    double L, dx;                                 // L, length of domain, dx.
    double omega_lower_ansatz, omega_ansatz, omega_true; // omega_bisection
    double OM_INF, PSI_INF; // asymptotics of lapse and cpnformal factpr
    int matching_index;     // integer where growing mode becomes relevant
    double eps = 10e-20;
    double omega_tolerance = 1e-20;
    int niter;           // number of iterations for finding the solution
    bool use_own_ansatz; // whether to set own custom upper omega

    void rk4(const double ww_);
    void rk4_asymp(const int iter, const bool adaptive, const double ww_);
    void rk4_match(const int iter, const bool adaptive, const double ww_);
    double small_A_RHS(const double x, const double A, const double DA,
                       const double PSI, const double DPSI, const double OM,
                       const double ww_);
    double A_RHS(const double x, const double A, const double DA,
                 const double PSI, const double DPSI, const double OM,
                 const double ww_);
    double DA_RHS(const double x, const double A, const double DA,
                  const double PSI, const double DPSI, const double OM,
                  const double ww_);
    double PSI_RHS(const double x, const double A, const double DA,
                   const double PSI, const double DPSI, const double OM,
                   const double ww_);
    double DPSI_RHS(const double x, const double A, const double DA,
                    const double PSI, const double DPSI, const double OM,
                    const double ww_);
    double OMEGA_RHS(const double x, const double A, const double DA,
                     const double PSI, const double DPSI, const double OM,
                     const double ww_);
    void initialise();
    double amplitude_criterion(double omega, double index);
    int zero_crossings();
    void truncate_solution();
    void force_flat(const int iter_crit);
    double bisect_omega(double omega_min, double omega_max);
    int found_zero_crossing();
    double find_upper_omega();
    int find_matching_index();
    int iofr(double rtarget);
    double V(const double A);
    double DV(const double A);
    void calculate_aspect_mass();
    void calculate_adm_mass();
    double calculate_radius();

  public:
    BosonStarSolution();

    std::vector<double> A;                  // scalar field modulus
    std::vector<double> dA;                 // scalar field modulus gradient
    std::vector<double> psi;                // conformal factor
    std::vector<double> dpsi;               // conformal factor gradient
    std::vector<double> omega;              // lapse
    std::vector<double> radius_array;       // isotropic radius
    std::vector<double> areal_radius_array; // areal radius
    std::vector<double> boson_mass;         // mass
    std::vector<double> adm_mass;           // mass
    double radius,
        compactness_value; // for storing the radius and compactness of the BS
    bool initialise_from_data_file;

    void set_initialcondition_params(BosonStar_params_t m_params_BosonStar,
                                     Potential::params_t m_params_potential,
                                     const double max_r);
    void initialise_from_file();
    double interpolate_vars(const double r, std::vector<double> vars) const;
    double get_A_interp(const double r) const;
    double get_lapse_interp(const double r) const;
    double get_psi_interp(const double r) const;
    double get_dA_interp(const double r) const;
    double get_dlapse_interp(const double r) const;
    double get_dpsi_interp(const double r) const;
    double get_BSfrequency() const;
    void output_csv();
    void main();
};

#include "BosonStarSolution.impl.hpp"

#endif /* BOSONSTARSOLUTION_HPP_ */
