#ifndef WEIGHTFUNCTION_HPP_
#define WEIGHTFUNCTION_HPP_

#include <array>

class WeightFunction
{
    public:

        double m_separation;
        double m_x_BS, m_y_BS, m_z_BS;
        int m_order;

        WeightFunction(double separation, double x_BS, double y_BS, double z_BS , int n){
            m_separation = separation;
            m_x_BS = x_BS;
            m_y_BS = y_BS;
            m_z_BS = z_BS;
            m_order = n;
        }

        double profile_chi()
        {
            double dist_BS = sqrt(m_x_BS*m_x_BS + m_y_BS*m_y_BS + m_z_BS*m_z_BS);
            double numer = m_separation * pow(m_separation - dist_BS, m_order);
            double denom = pow(m_separation, m_order + 1) + pow(dist_BS, m_order + 1);
            return numer / denom;
        }

        double profile_chi_2(double x, double y, double z, double R)
        {
            double r = sqrt(x*x + y*y + z*z);
            double ratio = r / R;
            return 1 - tanh( pow(ratio, m_order) );
        }
};

class WeightFunctionAngle
{
    public:

        double m_x_sep, m_y_sep, m_distance;
        double m_x_BS, m_y_BS, m_z_BS, m_r_BS;
        double m_epsilon, m_R;
        int m_order;

        WeightFunctionAngle(double separation, double impact_param, 
                       double x_BS, double y_BS, double z_BS , 
                       double eps, double R){
            m_x_sep = separation;
            m_y_sep = impact_param;
            m_distance = sqrt(m_x_sep*m_x_sep + m_y_sep*m_y_sep)
            m_x_BS = x_BS;
            m_y_BS = y_BS;
            m_z_BS = z_BS;
            m_r_BS = sqrt(m_x_BS*m_x_BS + m_y_BS*m_y_BS + m_z_BS*m_z_BS);
            m_epsilon = eps;
            m_R = R;
            m_order = n;
        }

        double profile_chi()
        {
            double correction = (-m_x_sep * x_BS + m_y_sep * y_BS) / m_distance;
            double arg = m_R*m_R + m_r_BS*m_r_BS + m_epsilon * correction;
            return 1 / sqrt(arg);
        }
};

#endif /* WEIGHTFUNCTION_HPP_ */
