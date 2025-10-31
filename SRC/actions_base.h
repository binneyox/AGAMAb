/** \file    actions_base.h
    \brief   Base classes for actions, angles, and action/angle finders
    \author  Eugene Vasiliev
    \date    2015
*/
#pragma once
#include "coord.h"

/** Classes and routines for transformations between position/velocity and action/angle phase spaces */
namespace actions {

/** Actions in arbitrary potential */
struct Actions {
    double Jr;       ///< radial action or its analog, [0..infinity)
    double Jz;       ///< vertical action or its analog, [0..infinity)
    double Jphi;     ///< azimuthal action (equal to the z-component of angular momentum in
                     ///< axisymmetric case, can have any value)
    Actions() {};
    Actions(double _Jr, double _Jz, double _Jphi) : Jr(_Jr), Jz(_Jz), Jphi(_Jphi) {};
    Actions& operator *= (const double a){
	    Jr*=a; Jz*=a; Jphi*=a;
	    return *this;
    }
    const Actions operator * (const double a) const{
	    return Actions(Jr*a, Jz*a, Jphi*a);
    }
    const Actions operator / (const double a) const{
	    return Actions(Jr/a, Jz/a, Jphi/a);
    }
    Actions& operator += (const Actions Jp){
	    Jr += Jp.Jr; Jz += Jp.Jz; Jphi += Jp.Jphi;
	    return *this;
    }
    const Actions operator + (const Actions Jp) const{
	    return Actions(Jr+Jp.Jr, Jz+Jp.Jz, Jphi+Jp.Jphi);
    }
};

/** Frequencies of motion (Omega = dH/dJ) */
struct Frequencies {
	double Omegar;    ///< frequency of radial motion, dH/dJr
	double Omegaz;    ///< frequency of vertical motion, dH/dJz
	double Omegaphi;  ///< frequency of azimuthal motion, dH/dJphi

	Frequencies() {};
	Frequencies(double omr, double omz, double omphi) : Omegar(omr), Omegaz(omz), Omegaphi(omphi) {};
	Frequencies& operator *= (const double a){
		Omegar*=a; Omegaz*=a; Omegaphi*=a;
		return *this;
	}
	const Frequencies operator * (const double a) const{
		return Frequencies(Omegar*a, Omegaz*a, Omegaphi*a);
	}
	const Frequencies operator / (const double a) const{
		return Frequencies(Omegar/a, Omegaz/a, Omegaphi/a);
	}
	Frequencies& operator += (const Frequencies fr){
		Omegar += fr.Omegar; Omegaz += fr.Omegaz; Omegaphi += fr.Omegaphi;
		return *this;
	}
	const Frequencies operator + (const Frequencies fr) const{
		return Frequencies(Omegar+fr.Omegar, Omegaz+fr.Omegaz,
				   Omegaphi+fr.Omegaphi);
	}
};

/** Angles in arbitrary potential */
struct Angles {
    double thetar;   ///< phase angle of radial motion
    double thetaz;   ///< phase angle of vertical motion
    double thetaphi; ///< phase angle of azimuthal motion

    Angles() {};
    Angles(double tr, double tz, double tphi) : thetar(tr), thetaz(tz), thetaphi(tphi) {};
    Angles(Frequencies Om) : thetar(Om.Omegar), thetaz(Om.Omegaz), thetaphi(Om.Omegaphi) {};
    Angles& operator *= (const double a){
	    thetar*=a; thetaz*=a; thetaphi*=a;
	    return *this;
    }
    const Angles operator * (const double a) const{
	    return Angles(thetar*a, thetaz*a, thetaphi*a);
    }
    Angles& operator += (const Angles thetap){
	    thetar += thetap.thetar; thetaz += thetap.thetaz; thetaphi += thetap.thetaphi;
	    return *this;
    }
    const Angles operator + (const Angles thetap) const{
	    return Angles(thetar+thetap.thetar, thetaz+thetap.thetaz, thetaphi+thetap.thetaphi);
    }
};

/** A combination of both actions and angles */
struct ActionAngles: Actions, Angles {
    ActionAngles() {};
    ActionAngles(const Actions& acts, const Angles& angs) : Actions(acts), Angles(angs) {};
    ActionAngles(const double Jr,const double Jz,const double Jphi,
		 const double thetar,const double thetaz,const double thetaphi):
	    Actions(Jr,Jz,Jphi), Angles(thetar,thetaz,thetaphi){}
};


/** Derivatives of coordinate/momentum variables w.r.t actions:
    each of three member fields stores the derivative of 6 pos/vel elements by the given action,
    in an inverted notation:  e.g.,  d(v_phi)/d(J_z) = dbyJz.vphi */
template <typename CoordT> struct DerivAct {
    coord::PosMomT<CoordT> dbyJr, dbyJz, dbyJphi;
};
typedef struct EXP DerivAct<coord::Cyl> DerivActCyl;

/** Derivatives of coordinate/momentum variables w.r.t angles:
    each of three member fields stores the derivative of 6 pos/vel elements by the given angle */
template <typename CoordT> struct DerivAng {
    coord::PosMomT<CoordT> dbythetar, dbythetaz, dbythetaphi;
};
typedef struct EXP DerivAng<coord::Cyl> DerivAngCyl;

/** Base class for action finders, which convert position/velocity pair to action/angle pair */
class EXP BaseActionFinder{
public:
    BaseActionFinder() {};
    virtual ~BaseActionFinder() {};

    /** Evaluate actions for a given position/velocity point in cylindrical coordinates */
    virtual Actions actions(const coord::PosVelCyl& point) const = 0;

    /** Evaluate actions and angles for a given position/velocity point in cylindrical coordinates;
        if the output argument freq!=NULL, also store the frequencies */
    virtual ActionAngles actionAngles(const coord::PosVelCyl& point, Frequencies* freq=NULL) const = 0;

private:
    /// disable copy constructor and assignment operator
    BaseActionFinder(const BaseActionFinder&);
    BaseActionFinder& operator= (const BaseActionFinder&);
};

/** Base class for action/angle mappers, which convert action/angle variables to position/velocity point */
class EXP BaseActionMapper{
public:
    BaseActionMapper() {};
    virtual ~BaseActionMapper() {};

    /** Map a point in action/angle space to a position/velocity in physical space;
        if the output argument freq!=NULL, also store the frequencies */
    virtual coord::PosVelCyl map(const ActionAngles& actAng, Frequencies* freq=NULL) const = 0;
private:
    /// disable copy constructor and assignment operator
    BaseActionMapper(const BaseActionMapper&);
    BaseActionMapper& operator= (const BaseActionMapper&);
};

/** Base class for canonical maps in action/angle space, which transform from one set of a/a
    variables to another one */
class EXP BaseCanonicalMap{
public:
    BaseCanonicalMap() {};
    virtual ~BaseCanonicalMap() {};

    virtual unsigned int numParams() const = 0;

    /** Map a point in action/angle space to a point in another action/angle space */
    virtual ActionAngles map(const ActionAngles& actAng) const = 0;
private:
    /// disable copy constructor and assignment operator
    BaseCanonicalMap(const BaseCanonicalMap&);
    BaseCanonicalMap& operator= (const BaseCanonicalMap&);
};

/** Base class for toy maps used in torus machinery, which provide conversion from action/angle
    to coordinate/momentum variables, and also provide the derivatives of this transformation */
/*class BaseToyMap{
public:
    virtual ~BaseToyMap() {};
    virtual coord::PosMomT<CoordT> map(
        const ActionAngles& actAng,
        Frequencies* freq=NULL,
        DerivAct<CoordT>* derivAct=NULL,
        DerivAng<CoordT>* derivAng=NULL,
        coord::PosMomT<CoordT>* derivParam=NULL) const = 0;
};
//*/
/** Base class for point transformations that map canonically conjugate coordinate/momentum
    in some intrinsic coord system into position/momentum in cylindrical coordinates */

class EXP BasePointTransform{
public:
    virtual ~BasePointTransform() {};
    /** maps in both given coordinate system (either cylindrical or spherical) and other coordinate system. 
    In other coordinate system there is the usual transformation from spherical or cylindrical to given coordinate system
    before mapping*/
    virtual coord::PosMomCyl map(const coord::PosMomCyl &point) const=0;
    virtual coord::PosMomCyl map(const coord::PosMomSph &point) const=0;
    /*map as well as derivates of transformated phase space point w.r.t original phase space point. 
    Note derivative is always in given coordinate system*/
    virtual coord::PosMomCyl map(const coord::PosMomCyl &point, math::Matrix<double>& dRzdrt) const=0;
    virtual coord::PosMomCyl map(const coord::PosMomSph &point, math::Matrix<double>& dRzdrt) const=0;
    //reverse map.
    virtual coord::PosMomCyl revmap(const coord::PosMomCyl &point) const=0;
    virtual coord::PosMomSph revmapSph(const coord::PosMomCyl &point) const=0;
    //Gives parameters and/or fourier series.
    virtual void getParams(double* params=NULL, double* Fr=NULL, double* Fz=NULL) const=0;
    //Name of PT:either Harmonic oscillator or Isochrone
    virtual const char* name() const=0;
    //gives number of parameters for Point Transformation not counting fourier series.
    virtual int numParams() const=0;
    //number of coefficients in fourier series expansion inr r.
    virtual int FourierSizer() const=0;
    //number of coefficients in fourier series expansion in v.
    virtual int FourierSizez() const=0;
};
typedef std::shared_ptr<const BasePointTransform> PtrPointTransform;
enum ToyPotType {None, Is, HO};
class EXP BaseToyMap{
    public:
    virtual ~BaseToyMap() {};
    virtual Actions pq2J(const coord::PosMomCyl Rzp) const=0;
    virtual ActionAngles pq2aa(const coord::PosMomCyl& Rz) const=0;
    virtual coord::PosMomCyl from_aaT(const ActionAngles& aaT) const=0;
    virtual coord::PosMomCyl from_aaT(const Actions& J, const Angles& theta) const=0;
    virtual coord::PosMomCyl from_aaT(const ActionAngles& aaT, DerivAct<coord::Cyl>& dRzdJ) const=0;
    virtual coord::PosMomCyl from_aaT(const ActionAngles& aaT, DerivAng<coord::Cyl>& dRzdT) const=0;
    virtual coord::PosMomCyl from_aaT(const Actions& J, const Angles& theta, DerivAct<coord::Cyl>& dRzdJ) const=0;
    virtual coord::PosMomCyl from_aaT(const Actions& J, const Angles& theta, DerivAng<coord::Cyl>& dRzdt) const=0;
	//returns number of parameters in point transformation
    virtual PtrPointTransform getPointTrans() const=0;
    virtual ToyPotType getToyMapType() const=0;
	virtual int PTParamSize() const =0;
	virtual int potParamSize() const=0;
	virtual const char* name() const=0;
    // Gets parameters of Toy potential. For isochrone x=Js, y=b. For HO x=omegar, y=omegaz; 
    virtual void getParams(double* x=NULL,double* y=NULL) const=0;
};
typedef std::shared_ptr<const BaseToyMap> PtrToyMap;
}  // namespace action