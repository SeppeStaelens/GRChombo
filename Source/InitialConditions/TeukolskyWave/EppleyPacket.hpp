
#ifndef EPPLEYPACKET_HPP_
#define EPPLEYPACKET_HPP_

#include "cmath"
#include "EppleyPacketParams.hpp"
#include <vector>

class EppleyPacket
{
  private:          // private member variables/arrays
    double amplitude; //!< amplitude of the packet
    double sigma; //!< width of the packet

    //! F function and its derivatives, where x = r \pm t
    double get_F(double x) const;
    double get_Fd1(double x) const;
    double get_Fd2(double x) const;
    double get_Fd3(double x) const;
    double get_Fd4(double x) const;

    //! Auxiliary functions
    double get_A(double r, double t, int sign) const;
    double get_B(double r, double t, int sign) const;
    double get_C(double r, double t, int sign) const;

  public:
    double get_A_tot(double r, double t) const;
    double get_B_tot(double r, double t) const;
    double get_C_tot(double r, double t) const;

    EppleyPacket(EppleyPacket_params_t m_params);
};

class EppleyPacketM0 : public EppleyPacket
{
  public:
    EppleyPacketM0(EppleyPacket_params_t m_params) : EppleyPacket(m_params) {}

    //! Metric functions
    double get_gxx(double x, double y, double z, double r, double t) const;
    double get_gxy(double x, double y, double z, double r, double t) const;
    double get_gxz(double x, double y, double z, double r, double t) const;
    double get_gyy(double x, double y, double z, double r, double t) const;
    double get_gyz(double x, double y, double z, double r, double t) const;
    double get_gzz(double x, double y, double z, double r, double t) const;
};


#include "EppleyPacket.impl.hpp"

#endif /* EPPLEYPACKET_HPP_ */
