/** \file    potential_interpolators
    \brief   Parameters of shell orbits for use with the Staeckel Fudge
    and Torus generation
    \author  Eugene Vasiliev, James binney & Tom Wright
    \date    2015 - 2025

   Routines in this file determine the radii Rsh at which shell orbits cut
   the meridional plane and the focal distances of the ellipses that
   fit them. These numbers are needed for the Staeckel
   Fudge actionFinder and the TorusGenerator. Values obtained at points
   on a grid in E and Xip are packaged into
   instances of lineartInterpolators2d in  different ways for the
   actionFinder and the TorusGenerator: the former uses (scaledE,scalesXip)
   the latter uses (scaledL,scaledXip) where L=Jz+|Jphi|. For the Fudge
   we define Xip=|Jphi|/Lcirc(E), while for the TorusGenerator we
   define Xip=|Jphi|/(Jz+|Jphi|) (because we don't know E)
   Interpolation is done in scaled variables. Finally we store Jz(Jf)
   along the box/loop transition in a PolarInterpolator. 
*/
#pragma once
#include <iostream>
#include "math_core.h"
#include "math_fit.h"
#include "math_spline.h"
#include "math_ode.h"
#include "utils.h"
#include "orbit.h"
#include "potential_base.h"
#include "potential_utils.h"
#include "actions_base.h"

namespace potential{

/// return scaledE as a function of E and invPhi0 = 1/Phi(0)
//inline double scaleE(const double E, const double invPhi0) { return log(invPhi0 - 1/E); }

		//scaling of E for finding planar orbit energy
inline double scaleE(const double E, const double invPhi0, /*output*/ double* dEdscaledE=NULL)
{
	double expX = invPhi0 - 1/E;
	if(dEdscaledE)
		*dEdscaledE = E * E * expX;
	return log(expX);
}
		//return E and optionally dE/d(scaledE) as a function of scaledE and invPhi0
inline double unscaleE(const double scaledE, const double invPhi0, /*output*/ double* dEdscaledE=NULL)
{
	double expX = exp(scaledE);
	double E = 1 / (invPhi0 - expX);
	if(dEdscaledE)
		*dEdscaledE = E * E * expX;
	return E;
}

/// function to be used in root-finder for locating a shell orbit in R-z plane
class EXP FindClosedOrbitRZplane: public math::IFunction {
	public:
		FindClosedOrbitRZplane(const potential::BasePotential& p, 
				       double _E, double _Jphi, double _Rmin, double _Rmax,
				       double& _timeCross, std::vector<std::pair<coord::PosVelCyl,double> >& _traj,
				       double& _Jz)
				:
		    poten(p), E(_E), Jphi(_Jphi), Rmin(_Rmin), Rmax(_Rmax), 
		    Jz(_Jz), timeCross(_timeCross), traj(_traj)
		{}
    /// report the difference in R between starting point (R0, z=0, vz>0)
    /// and return point (Rcross, z=0, vz<0)
		    virtual void evalDeriv(const double R0, double* val, double* der, double*) const;
		    virtual unsigned int numDerivs() const { return 1; }
	private:
		const potential::BasePotential& poten;
		const double E, Jphi;               ///< parameters of motion in the R-z plane
		const double Rmin, Rmax;          ///< boundaries of interval in R (to skip the first two calls)
		double& Jz;
		double& timeCross;                ///< keep track of time required to complete orbit
		std::vector<std::pair<coord::PosVelCyl,double> >& traj; ///< store the trajectory
};

/// function to be used in root-finder for locating a shell orbit in R-z plane
class EXP FindRzClosedOrbitV: public math::IFunction {
	public:
		FindRzClosedOrbitV(const potential::BasePotential& p, double _R0,
				   double _Jphi, double& _timeCross,
				   std::vector<std::pair<coord::PosVelCyl,double> >& _traj,
				   double& _Jz)
				:
		    poten(p), R0(_R0), Jphi(_Jphi), 
		    timeCross(_timeCross), traj(_traj), Jz(_Jz)
		{}
    /// report the difference in R between starting point (R0, z=0, vz>0)
    /// and return point (Rcross, z=0, vz<0)
		    virtual void evalDeriv(const double V0, double* val, double* der, double*) const;
		    virtual unsigned int numDerivs() const { return 1; }
	private:
		const potential::BasePotential& poten;
		const double R0, Jphi;               ///< parameters of motion in the R-z plane
		double& timeCross;                ///< keep track of time required to complete orbit
		std::vector<std::pair<coord::PosVelCyl,double> >& traj; ///< store the trajectory
		double& Jz;
};

void findCrossingPointV(
			const potential::BasePotential& poten, double R0, double Jphi, double V0,
			double& timeCross, std::vector<std::pair<coord::PosVelCyl,double> >& traj,
			double& Rcross, double& dRcrossdV0, double& Jz);

/* Returns focal distance of the shell orbit with given E & Jphi. Returns Rsh and
 * optionally Jz and the entire orbit
 */
EXP double estimateFocalDistanceShellOrbit(
					   const potential::BasePotential& poten, double E, double Jphi, 
					   double* R=0, double* Jz=NULL, std::vector<coord::PosVelCyl>* shell=NULL);

/* To construct a ShellInterpolator we integrate shell orbits on grid
 * in E,Xi=Jphi/Jc(E) (energy and inclination)
 * For each of D,Rsh we produce 2 types of LinearInterpolator2d:
 * For actionFinder x,y are scaledE and scaledXi
 * For torusGenerator x,y are L,Xi
 * Note: Xi has different meaning in the 2 cases:
    for actionFinder Xi=Jphi/Lc(E)
    for TorusGenerator Xi=Jphi/(Jz+Jphi)
*/

class EXP ShellInterpolator{
	private:
		math::LinearInterpolator2d interpDE, interpRE;
		math::LinearInterpolator2d interpDL, interpRL;
	public:
		ShellInterpolator(){}
		ShellInterpolator(const BasePotential&);
		ShellInterpolator(const math::LinearInterpolator2d& _intDE,
					   const math::LinearInterpolator2d& _intRE,
					   const math::LinearInterpolator2d& _intDL,
					   const math::LinearInterpolator2d& _intRL):
		    interpDE(_intDE), interpRE(_intRE),  interpDL(_intDL), interpRL(_intRL) {
		}
		
		void getRshDelta(const double L,const double Xip,
				 double& Rsh, double& Delta) const {
			const math::ScalingSemiInf scalingL;
			const double scaledL = math::clip(scale(scalingL,L),
				interpDL.xmin(), interpDL.xmax());
			const math::ScalingCub scalingXi(0,1);
			const double scaledXip = math::clip(scale(scalingXi, Xip), 0., 1.);
			Rsh = interpRL.value(scaledL, scaledXip);
			Delta = interpDL.value(scaledL, scaledXip);
		}
		void getRshDelta(const double E, const double Xip,
				 const double invPhi0, double& Rsh, double& Delta) const {
			const double scaledE = math::clip(scaleE(E, invPhi0),
				interpDE.xmin(), interpDE.xmax());
			const math::ScalingCub scalingXi(0,1);
			const double scaledXip = math::clip(scale(scalingXi, Xip), 0., 1.);
			Rsh = interpRE.value(scaledE, scaledXip);
			Delta = interpDE.value(scaledE, scaledXip);
		}
		double getRsh(const double L, const double Xip) const{
			double Rsh, Delta;
			getRshDelta(L, Xip, Rsh, Delta);
			return Rsh;
		}
		double getRsh(const double E, const double Xip,
				const double invPhi0) const{
			double Rsh, Delta;
			getRshDelta(E, Xip, invPhi0, Rsh, Delta);
			return Rsh;
		}
		double getDelta(const double L, const double Xip) const{
			double Rsh, Delta;
			getRshDelta(L, Xip, Rsh, Delta);
			return Delta;
		}
		double getDelta(const double E, const double Xip,
				const double invPhi0) const{
			double Rsh, Delta;
			getRshDelta(E, Xip, invPhi0, Rsh, Delta);
			return Delta;
		}
		std::pair<double,double> getXiChi(const double E,const double Xip, const double invPhi0) const{
			const double scaledE = math::clip(scaleE(E, invPhi0),
				interpDE.xmin(), interpDE.xmax());
			const math::ScalingCub scalingXi(0,1);
			const double scaledXip = math::clip(scale(scalingXi, Xip), 0., 1.);
			return std::make_pair(scaledE,scaledXip);
		}
};

/* The main job of PolarInterpolator is to hold the curve Jz(Jf) along
 * which the box/loop transition lies. In adition it holds the values
 * of Delta(E) that cause the I3 centrifugal barrier to vanish on this
 * curve and te associaed I3(E), where I3 is computed from the velocity
 * of the transition orbit at (Rsh,0) 
*/
class EXP  PolarInterpolator{
	private:
		math::LinearInterpolator interpI3, interpFD;
		std::vector<double> coeffsJz;
		math::ScalingSemiInf Sc;
	public:
		PolarInterpolator(){}
		PolarInterpolator(const BasePotential&, const PtrShellInterpolator);
		PolarInterpolator(const std::vector<double>& gridEscaled, const std::vector<double>& gridI3,
				  const std::vector<double>& gridFD, const std::vector<double>& gridJfScaled,
				  const std::vector<double>& gridJz) :
		    interpI3(gridEscaled,gridI3), interpFD(gridEscaled,gridFD),
		    coeffsJz(math::fitPoly(15,gridJfScaled,gridJz)) {}
		void getDeltaI3(const double E,const double invPhi0,
			       double& Delta, double& I3) const{
			double scaledE = math::clip(scaleE(E, invPhi0),
				interpFD.xmin(), interpFD.xmax());
			Delta = interpFD.value(scaledE);
			scaledE = math::clip(scaleE(E, invPhi0),
				interpI3.xmin(), interpI3.xmax());
			I3 = interpI3.value(scaledE);
		}			
		double getI3(const double E,const double invPhi0) const{
			const double scaledE = math::clip(scaleE(E, invPhi0),
				interpI3.xmin(), interpI3.xmax());
			return interpI3.value(scaledE);
		}
		double getDelta(const double E,const double invPhi0) const{
			const double scaledE = math::clip(scaleE(E, invPhi0),
				interpFD.xmin(), interpFD.xmax());
			return interpFD.value(scaledE);
		}
		double getJz(const double Jf) const{
			return math::evalPoly(coeffsJz, scale(Sc,Jf));
		}
};

}//namespace