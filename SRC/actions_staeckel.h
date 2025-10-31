/** \file    actions_staeckel.h
    \brief   Action-angle finders using Staeckel potential approximation
    \author  Eugene Vasiliev
    \date    2015

Computation of actions and angles for axisymmetric systems, using two related methods:
1) exact expressions for Staeckel potentials  (e.g. for the Perfect Ellipsoid),
2) "Staeckel Fudge" approximation applicable for any axisymmetric potential.
To avoid confusion, the latter is called simply "Fudge" in what follows.

Most sub-steps are shared between the two methods; 
only the computation of integrals of motion, and the auxiliary function that enters 
the expression for canonical momentum, are specific to each method.

The implementation is inspired by the code written by Jason Sanders,
but virtually nothing of the original code remains.
*/
#pragma once
#include "actions_base.h"
#include "actions_focal_distance_finder.h"
#include "smart.h"

namespace actions {

/** aggregate class that contains the point in prolate spheroidal coordinates,
    integrals of motion, and reference to the potential,
    shared between Axisymmetric Staeckel  and Axisymmetric Fudge action finders;
    only the coordinates and the two classical integrals are in common between them.
    It also implements the IFunction interface, providing the "auxiliary function" 
    that is used in finding the integration limits and computing actions and angles;
    this function F(tau) is related to the canonical momentum p(tau)  as 
    \f$  p(tau)^2 = F(tau) / (2*(tau+alpha)^2*(tau+gamma))  \f$,
    and the actual implementation of this function is specific to each descendant class.
*/
class EXP AxisymFunctionBase: public math::IFunction {
	public:
		coord::PosVelProlSph point; ///< position/derivative in prolate spheroidal coordinates
		const double E;             ///< total energy
		const double Lz;            ///< z-component of angular momentum
		const double I3;            ///< third integral
		AxisymFunctionBase(const coord::PosVelProlSph& _point, double _E, double _Lz, double _I3) :
		    point(_point), E(_E), Lz(_Lz), I3(_I3) {};
};

// ------ SPECIALIZED functions for Staeckel action finder -------

/** parameters of potential, integrals of motion, and prolate spheroidal coordinates 
    SPECIALIZED for the Axisymmetric Staeckel action finder */
class EXP AxisymFunctionStaeckel: public AxisymFunctionBase {
	public:
		const math::IFunction& fncG;  ///< single-variable function of a Staeckel potential
		AxisymFunctionStaeckel(const coord::PosVelProlSph& _point, double _E, double _Lz, double _I3,
				       const math::IFunction& _fncG) :
		    AxisymFunctionBase(_point, _E, _Lz, _I3), fncG(_fncG) {};

    /** auxiliary function that enters the definition of canonical momentum for 
        for the Staeckel potential: it is the numerator of eq.50 in de Zeeuw(1985);
        the argument tau is replaced by tau+gamma >= 0. */
		virtual void evalDeriv(const double tauplusgamma, 
				       double* value=0, double* deriv=0, double* deriv2=0) const;
		virtual unsigned int numDerivs() const { return 2; }
};

/** compute integrals of motion in the Staeckel potential of an oblate perfect ellipsoid, 
    together with the coordinates in its prolate spheroidal coordinate system 
*/
EXP AxisymFunctionStaeckel findIntegralsOfMotionOblatePerfectEllipsoid(
	const potential::OblatePerfectEllipsoid& potential, 
	const coord::PosVelCyl& point);
/// \name  ------- Stand-alone driver routines that compute actions for a single point -------
///@{

/** Find exact actions in the Staeckel potential of oblate Perfect Ellipsoid.
    \param[in]  potential is the input Staeckel potential.
    \param[in]  point     is the position/velocity point.
    \return     actions for the given point, or Jr=Jz=NAN if the energy is positive.
*/
EXP Actions actionsAxisymStaeckel(
	const potential::OblatePerfectEllipsoid& potential, 
	      const coord::PosVelCyl& point);

/** Find exact actions and angles in the Staeckel potential of oblate Perfect Ellipsoid.
    \param[in]  potential is the input Staeckel potential.
    \param[in]  point     is the position/velocity point.
    \param[out] freq      if not NULL, store the frequencies of motion in this variable.
    \return     actions and angles for the given point, or Jr=Jz=NAN if the energy is positive.
*/
EXP ActionAngles actionAnglesAxisymStaeckel(
    const potential::OblatePerfectEllipsoid& potential, 
    const coord::PosVelCyl& point,
    Frequencies* freq=NULL);

/** Find approximate actions in a given axisymmetric potential, using the Staeckel Fudge method.
    \param[in]  potential is the arbitrary axisymmetric potential.
    \param[in]  point     is the position/velocity point.
    \param[in]  focalDistance is the geometric parameter of best-fit coordinate system:
    the accuracy of the method depends on this parameter, which should be estimated by one of 
    the methods from actions_focal_distance_finder.h.
    \return     actions for the given point, or Jr=Jz=NAN if the energy is positive.
    \throw      std::invalid_argument exception if the potential is not axisymmetric.
*/
EXP Actions actionsAxisymFudge(
    const potential::BasePotential& potential, 
    const coord::PosVelCyl& point,
    double focalDistance);
//    double lambda_dip = NAN);

/** Find approximate actions and angles in a given axisymmetric potential, 
    using the Staeckel Fudge method.
    \param[in]  potential is the arbitrary axisymmetric potential.
    \param[in]  point     is the position/velocity point.
    \param[in]  focalDistance is the geometric parameter of best-fit coordinate system.
    \param[out] freq      if not NULL, store the frequencies of motion in this variable.
    \return     actions and angles for the given point, or Jr=Jz=NAN if the energy is positive.
    \throw      std::invalid_argument exception if the potential is not axisymmetric.
*/
EXP ActionAngles actionAnglesAxisymFudge(
    const potential::BasePotential& potential, 
    const coord::PosVelCyl& point, 
    double focalDistance,
//    double lambda_dip=NAN,
    Frequencies* freq=NULL);

///@}
/// \name  ------- Class interface to action/angle finders  -------
///@{

/** Action/angle finder for an Oblate Perfect Ellipsoid potential. */
class EXP ActionFinderAxisymStaeckel: public BaseActionFinder {
public:
    explicit ActionFinderAxisymStaeckel(const potential::PtrOblatePerfectEllipsoid& potential) :
        pot(potential) {};

    virtual Actions actions(const coord::PosVelCyl& point) const {
        return actionsAxisymStaeckel(*pot, point); }

    virtual ActionAngles actionAngles(const coord::PosVelCyl& point, Frequencies* freq=NULL) const {
        return actionAnglesAxisymStaeckel(*pot, point, freq); }

private:
    const potential::PtrOblatePerfectEllipsoid pot;  ///< the potential in which actions are computed
};

/** Action/angle finder for a generic axisymmetric potential, based on Staeckel Fudge approximation.
    It is more suitable for massive computation in a fixed potential than just using 
    the standalone routines, because it estimates the focal distance using a pre-computed 
    interpolation grid, rather than doing it individually for each point. This results in 
    a considerable speedup in action computation, for a minor overhead during initialization.
    Additionally, it may set up an interpolation table for actions as functions of three integrals
    of motion (one of them being approximate), which speeds up the evaluation by another order of
    magnitude, for a moderate decrease in accuracy. Interpolated actions have small but non-negligible
    systematic bias, which depends on the potential as well as the phase-space location,
    hence they cannot be used for comparing the likelihood of a DF in different potentials.
*/
class EXP ActionFinderAxisymFudge: public BaseActionFinder {
public:
    /** set up the action finder: an interpolator for the focal distance, and optionally
        interpolators for actions (if interpolate==true).
        \throw std::invalid_argument exception if the potential is not axisymmetric
        or std::runtime_error in case of other problems in initialization.
    */
    ActionFinderAxisymFudge(const potential::PtrPotential& potential, bool interpolate = false);

    virtual Actions actions(const coord::PosVelCyl& point) const;

    virtual ActionAngles actionAngles(const coord::PosVelCyl& point, Frequencies* freq=NULL) const;

    /** return the best-suitable focal distance for the given point, obtained by interpolation */
    double focalDistance(const coord::PosVelCyl& point) const;

private:
    const double invPhi0;                   ///< 1 / Phi(r=0)
    const potential::PtrPotential pot;      ///< the potential Phi in which actions are computed
    const potential::Interpolator interp;   ///< 1d interpolator for Lcirc(E)
    math::LinearInterpolator2d interpD;     ///< 2d interpolator for the focal distance Delta(E,Lz)
//   math::LinearInterpolator3d interpD;     ///< 3d interpolator for the focal distance Delta(E,Lz)
//    math::LinearInterpolator3d interpDip;   ///< 3d interpolator for min lam_dip in f(lam) contributed by I3
    math::CubicSpline2d interpR;            ///< 2d interpolator for Rshell(E,Lz) / Rcirc(E)
    math::CubicSpline3d intJr, intJz;       ///< 3d interpolators for Jr and Jz as functions of (E,Lz,I3)
};
// Class to find the largest FD that yields a centrifugal barrier even
// with Lz=0. Rsh is the radius of the shell orbit and vR is vR is
// what's required to create the orbit of interest shooting from
// (Rsh,0). Delta is an estimate of FD. The method bestFD() returns the
// FD at second argument; the first is the value of u at which p_u^2=0
//
class EXP FDfinder{
	private:
		double E,Rsh,vR,P0,Delta0,EmP0;
		const potential::PtrPotential pot;
	public:
		FDfinder(const double _E,const double _Rsh,const double _vR,
			 const double _Delta,
			 const potential::PtrPotential _pot) :
		    E(_E), Rsh(_Rsh), vR(_vR), Delta0(_Delta), pot(_pot){
			P0=pot->value(coord::PosCyl(Rsh,0,0));
			EmP0=E-P0;
		}
		math::Matrix<double> derivs(double u,double Delta,double& d2p2du2,
					    double* p2=NULL,double* p2prime=NULL);
		double bestFD(double&, double&);
};

///@}
}  // namespace actions
