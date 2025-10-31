/** \file    potential_analytic.h
    \brief   Several common analytic potential models
    \author  Eugene Vasiliev
    \date    2009-2015
*/
#pragma once
#include "potential_base.h"

namespace potential{

/** Spherical Plummer potential:
    \f$  \Phi(r) = - M / \sqrt{r^2 + b^2}  \f$. */
class EXP Plummer: public BasePotentialSphericallySymmetric{
public:
    Plummer(double _mass, double _scaleRadius) :
        mass(_mass), scaleRadius(_scaleRadius) {}
    virtual const char* name() const { return myName(); }
    static const char* myName() { static const char* text = "Plummer"; return text; }
    virtual double enclosedMass(const double radius) const;
    virtual double totalMass() const { return mass; }
private:
    const double mass;         ///< total mass  (M)
    const double scaleRadius;  ///< scale radius of the Plummer model  (b)

    /** Evaluate potential and up to two its derivatives by spherical radius. */
    virtual void evalDeriv(double r,
        double* potential, double* deriv, double* deriv2) const;
    double Phi(double r) const;
    double dPhidr(double r) const;
    double Vc2(double r) const{
	    return r*dPhidr(r);
    }
    /** explicitly define the density function, instead of relying on the potential derivatives
        (they suffer from cancellation errors already at r>1e5) */
    virtual double densitySph(const coord::PosSph &pos) const;
};

/** Spherical Isochrone potential:
    \f$  \Phi(r) = - M / (b + \sqrt{r^2 + b^2})  \f$. */
class EXP Isochrone: public BasePotentialSphericallySymmetric{
public:
    Isochrone(double _mass, double _scaleRadius) :
        mass(_mass), scaleRadius(_scaleRadius) {}
    virtual const char* name() const { return myName(); }
    static const char* myName() { static const char* text = "Isochrone"; return text; }
    virtual double totalMass() const { return mass; }
    virtual void evalDeriv(double r,
			   double* potential, double* deriv, double* deriv2) const;
    double density(double r) const;
    double Phi(double r) const;
    double dPhidr(double r) const;
    double Vc2(double r) const{
	    return r*dPhidr(r);
    }
    double Mass(double r) const;
private:
    const double mass;         ///< total mass  (M)
    const double scaleRadius;  ///< scale radius of the Isochrone model  (b)
};
/** Spherical Hernquist potential
 **/
class EXP Hernquist: public BasePotentialSphericallySymmetric{
	private:
		const double mass;
		const double scaleRadius;
		double rho0;
	public:
		Hernquist(double _M, double _a):
		    mass(_M), scaleRadius(_a),
		    rho0(mass/(2*M_PI*pow(scaleRadius,3))){}
		virtual const char* name() const { return myName(); }
		static const char* myName() { static const char* text = "Hernquist"; return text; }
		virtual double totalMass() const { return mass; }
		virtual void evalDeriv(double r,
				       double* potential, double* deriv, double* deriv2) const;
		double density(double r) const;
		double Phi(double r) const;
		double dPhidr(double r) const;
		double Vc2(double r) const{
			return r*dPhidr(r);
		}
		double Mass(double r) const;
		double a(void) const{
			return scaleRadius;
		}
		virtual double densitySph(const coord::PosSph &pos) const {
			return density(pos.r);
		}
};
/** Spherical Navarro-Frenk-White potential:
    \f$  \Phi(r) = - M \ln{1 + (r/r_s)} / r  \f$  (note that total mass is infinite and not M). */
class EXP NFW: public BasePotentialSphericallySymmetric{
	public:
		NFW(double _mass, double _scaleRadius) :
		    mass(_mass), scaleRadius(_scaleRadius) {}
		virtual const char* name() const { return myName(); }
		static const char* myName() { static const char* text = "NFW"; return text; }
		virtual double totalMass() const { return INFINITY; }
	private:
		const double mass;         ///< normalization factor  (M);  equals to mass enclosed within ~5.3r_s
		const double scaleRadius;  ///< scale radius of the NFW model  (r_s)
		virtual void evalDeriv(double r,
				       double* potential, double* deriv, double* deriv2) const;
};

/* Class for a spherical Phi(r)=r^a
 */
class EXP purePower: public BasePotentialSphericallySymmetric{
	public:
		purePower(double _index) :
		    index(_index) {}
		virtual const char* name() const { return myName(); }
		static const char* myName() { static const char* text = "purePower"; return text; }
		virtual double totalMass() const { return INFINITY; }
	private:
		const double index;         ///< powe-law index
		virtual void evalDeriv(double r,
				       double* potential, double* deriv, double* deriv2) const;
};

/** Axisymmetric Miyamoto-Nagai potential:
    \f$  \Phi(r) = - M / \sqrt{ R^2 + (A + \sqrt{z^2+b^2})^2 }  \f$. */
class EXP MiyamotoNagai: public BasePotentialCyl{
public:
    MiyamotoNagai(double _mass, double _scaleRadiusA, double _scaleRadiusB) :
        mass(_mass), scaleRadiusA(_scaleRadiusA), scaleRadiusB(_scaleRadiusB) {};
    virtual coord::SymmetryType symmetry() const { return coord::ST_AXISYMMETRIC; }
    virtual const char* name() const { return myName(); }
    static const char* myName() { static const char* text = "MiyamotoNagai"; return text; }
    virtual double totalMass() const { return mass; }
private:
    const double mass;         ///< total mass  (M)
    const double scaleRadiusA; ///< first scale radius  (A),  determines the extent in the disk plane
    const double scaleRadiusB; ///< second scale radius (B),  determines the vertical extent

    virtual void evalCyl(const coord::PosCyl &pos,
			 double* potential, coord::GradCyl* deriv, coord::HessCyl* deriv2) const;
};

/** Triaxial logarithmic potential:
    \f$  \Phi(r) = (1/2) \sigma^2 \ln[ r_c^2 + x^2 + (y/p)^2 + (z/q)^2 ]  \f$. */
class EXP Logarithmic: public BasePotentialCar{
public:
	Logarithmic(double sigma, double coreRadius=0, double axisRatioYtoX=1,
		    double axisRatioZtoX=1, double _Lm=10) :
	    sigma2(pow_2(sigma)), Rc2(pow_2(coreRadius)),
	    p2(pow_2(axisRatioYtoX)), q2(pow_2(axisRatioZtoX)), Lm(_Lm){}
	virtual coord::SymmetryType symmetry() const {
        return p2==1 ? (q2==1 ? coord::ST_SPHERICAL : coord::ST_AXISYMMETRIC) : coord::ST_TRIAXIAL; }
    virtual const char* name() const { return myName(); }
    static const char* myName() { static const char* text = "Logarithmic"; return text; }
    virtual double totalMass() const { return INFINITY; }
private:
    const double sigma2;       ///< squared asymptotic circular velocity (sigma)
    const double Rc2;          ///< squared core radius (r_c)
    const double p2;           ///< squared y/x axis ratio (p)
    const double q2;           ///< squared z/x axis ratio (q)
    const double Lm;		///< asymptotic value of 2*Phi/sigma^2

    virtual void evalCar(const coord::PosCar &pos,
        double* potential, coord::GradCar* deriv, coord::HessCar* deriv2) const;
};

/** Triaxial harmonic potential:
    \f$  \Phi(r) = (1/2) \Omega^2 [ x^2 + (y/p)^2 + (z/q)^2 ]  \f$. */
class EXP Harmonic: public BasePotentialCar{
	public:
		Harmonic(double Omega, double axisRatioYtoX=1, double axisRatioZtoX=1) :
		    Omega2(pow_2(Omega)), p2(pow_2(axisRatioYtoX)), q2(pow_2(axisRatioZtoX)) {}
		virtual coord::SymmetryType symmetry() const {
			return p2==1 ? (q2==1 ? coord::ST_SPHERICAL : coord::ST_AXISYMMETRIC) : coord::ST_TRIAXIAL; }
		virtual const char* name() const { return myName(); }
		static const char* myName() { static const char* text = "Harmonic"; return text; }
		virtual double totalMass() const { return INFINITY; }
	private:
		const double Omega2;       ///< squared oscillation frequency (Omega)
		const double p2;           ///< squared y/x axis ratio (p)
		const double q2;           ///< squared z/x axis ratio (q)

		virtual void evalCar(const coord::PosCar &pos,
				     double* potential, coord::GradCar* deriv, coord::HessCar* deriv2) const;
};

/** Henon & Heiles (1964) potential turned through 90 deg and z added:
    \f$  \Phi(r) = (1/2)[ x^2 + y^2 + (z/q)^2] + y^2x - x^3/3  \f$. */
class EXP HenonHeiles: public BasePotentialCyl{
	public:
		HenonHeiles(void) {}
		virtual coord::SymmetryType symmetry() const {
			return coord::ST_ZROTATION; }
		virtual const char* name() const { return myName(); }
		static const char* myName() { static const char* text = "HenonHeiles"; return text; }
		virtual double totalMass() const { return INFINITY; }
	private:
		virtual void evalCyl(const coord::PosCyl &pos,
				     double* potential, coord::GradCyl* deriv, coord::HessCyl* deriv2) const;
};

/** Potential from Binney (1982) illustrating irregular orbits
 **/
class EXP logRe: public BasePotentialCyl{
	private:
		const double sigma2;
		const double coreRadius2;
		const double q2;
		const double Re;
		const bool ifRe;
		virtual void evalCyl(const coord::PosCyl &pos, double* potential,
				     coord::GradCyl* deriv, coord::HessCyl* deriv2) const;
	public:
		logRe(double _sigma, double _Rc, double _q, double _Re=NAN) :
		    sigma2(_sigma*_sigma), coreRadius2(_Rc*_Rc), q2(_q*_q), 
		    Re(_Re), ifRe(!std::isnan(_Re)) {}
		virtual coord::SymmetryType symmetry() const {
			return coord::ST_ZROTATION; }
		virtual const char* name() const { return myName(); }
		static const char* myName() { static const char* text = "HenonHeiles"; return text; }
		virtual double totalMass() const { return INFINITY; }
};
/** Approximation to Sormani et al (2015 bar from Binney 2020
 **/
class EXP Sormani : public BasePotentialCyl {
	public:
		Sormani(const double rq, const double V0, const double Rb, const double q) :
		K(3*exp(2)/20.*pow_2(V0)*pow_3(rq)*0.4), q2(q*q), Rb2(Rb*Rb) {}
		virtual coord::SymmetryType symmetry() const {return coord::ST_TRIAXIAL;}
		virtual const char* name() const { return myName(); }
		static const char* myName() { static const char* text = "Sormani"; return text; }
		virtual double totalMass() const { return 0; }
	private:
		const double K,q2,Rb2;
		virtual void evalCyl(const coord::PosCyl &pos, double* potential,
				     coord::GradCyl* deriv, coord::HessCyl* deriv2) const;
//		A=3*exp(2)/20*pow_2(V0/Vc)*pow_3(rq/Rb)*0.4;
//		K=A*pos_2(Vc)*pow_3(Rb);
};
}//namespace