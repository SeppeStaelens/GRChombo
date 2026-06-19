
#ifndef EPPLEYPACKET_HPP_
#define EPPLEYPACKET_HPP_

#include "cmath"
#include "EppleyPacketParams.hpp"
#include <vector>

//! Base EppleyPacket class
class EppleyPacket
{
  private:          // private member variables/arrays
    double amplitude; //!< amplitude of the packet
    double sigma; //!< width of the packet

  public:
    //! F function and its derivatives, where x = r \pm t
    double get_F(double x) const;
    double get_Fd1(double x) const;
    double get_Fd2(double x) const;
    double get_Fd3(double x) const;
    double get_Fd4(double x) const;

    EppleyPacket(EppleyPacket_params_t m_params):
        amplitude(m_params.amplitude),
        sigma(m_params.sigma) {};
};

//! Base class for the even parity Eppley packets
class EvenEppleyPacket : EppleyPacket
{
  private:          // private member variables/arrays
    double amplitude; //!< amplitude of the packet
    double sigma; //!< width of the packet

    //! Auxiliary functions
    double get_A(double r, double t, int sign) const;
    double get_B(double r, double t, int sign) const;
    double get_C(double r, double t, int sign) const;

  public:
    double get_A_tot(double r, double t) const;
    double get_B_tot(double r, double t) const;
    double get_C_tot(double r, double t) const;

    EvenEppleyPacket(EppleyPacket_params_t m_params) : EppleyPacket(m_params) {};
};

//! Base class for the odd parity Eppley packets
class OddEppleyPacket : EppleyPacket
{
  private:          // private member variables/arrays
    double amplitude; //!< amplitude of the packet
    double sigma; //!< width of the packet

    //! Auxiliary functions
    double get_K(double r, double t, int sign) const;
    double get_L(double r, double t, int sign) const;

  public:
    double get_K_tot(double r, double t) const;
    double get_L_tot(double r, double t) const;

    OddEppleyPacket(EppleyPacket_params_t m_params) : EppleyPacket(m_params) {};
};

//! Specific Eppley Packet classes
//! m = 0 and m = 2 are the only ones implemented so far for the even parity case
class EppleyPacketM0 : public EvenEppleyPacket
{
  public:
    EppleyPacketM0(EppleyPacket_params_t m_params) : EvenEppleyPacket(m_params) {}

    //! Metric functions
    double get_gxx(double x, double y, double z, double r, double t) const;
    double get_gxy(double x, double y, double z, double r, double t) const;
    double get_gxz(double x, double y, double z, double r, double t) const;
    double get_gyy(double x, double y, double z, double r, double t) const;
    double get_gyz(double x, double y, double z, double r, double t) const;
    double get_gzz(double x, double y, double z, double r, double t) const;
};

class EppleyPacketM2 : public EvenEppleyPacket
{
  public:
    EppleyPacketM2(EppleyPacket_params_t m_params) : EvenEppleyPacket(m_params) {}

    //! Metric functions
    double get_gxx(double x, double y, double z, double r, double t) const;
    double get_gxy(double x, double y, double z, double r, double t) const;
    double get_gxz(double x, double y, double z, double r, double t) const;
    double get_gyy(double x, double y, double z, double r, double t) const;
    double get_gyz(double x, double y, double z, double r, double t) const;
    double get_gzz(double x, double y, double z, double r, double t) const;
};

//! m = 2 class for the odd parity case
class OddEppleyPacketM2 : public OddEppleyPacket
{
  public:
    OddEppleyPacketM2(EppleyPacket_params_t m_params) : OddEppleyPacket(m_params) {}

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
