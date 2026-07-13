/* DFs of form f(E)
 */
#pragma once
#include "potential_base.h"
#include "potential_analytic.h"
#include "math_spline.h"
#include "utils.h"
#include "actions_base.h"
#include "units.h"
#include <fstream>
#include <iostream>
#include <vector>

namespace df{
class EXP ergodicDF: public math::IFunctionNoDeriv{
	public:
		ergodicDF(void) {}
		virtual double value(const double E) const=0;
		//fraction of probability contributed as  v<Vmax
		virtual double fraction(const double Phi,const double Vmax) const;
};
class EXP BaseBWDF: public ergodicDF{
	public:
		BaseBWDF(){}
		virtual double value(const double E) const=0;
};
//Normalised to produce rho0 at r0
class EXP BahcallWolfDF_E: public BaseBWDF {
	private:
		const double MBH, rho0, r0, fac;
	public:
		BahcallWolfDF_E(const double _MBH, const double _rho0, const double _r0):
		    MBH(_MBH), rho0(_rho0), r0(_r0), fac(pow(.5*r0/MBH,1.5)) {}
		virtual double value(const double E) const{
			return E<0? rho0/M_PI*fac*pow(-E*r0/MBH,.25) : 0;
		}
};

//Normalised to produce rho0 at r0
class EXP BahcallWolfDF_Lc: public BaseBWDF {
	private:
		const double MBH, rho0, r0, fac, rcut;
		const math::LinearInterpolator& LcE;
	public:
		BahcallWolfDF_Lc(const double _MBH, const double _rho0,
				 const double _r0, const double _rcut,
				 math::LinearInterpolator& _LcE):
		    MBH(_MBH), rho0(_rho0), r0(_r0), rcut(_rcut),
		    LcE(_LcE), fac(pow(.5*r0/MBH,1.5)) {}
		virtual double value(const double E) const{
			//if(E>PhiCut) return 0;
			const double Lc=exp(LcE.value(E));//Lc in actual Phi
			double E0s=-.5*MBH*r0/pow_2(Lc);//corresponding dimless Kepler E
			double f=rho0/M_PI*fac*pow(-E0s,.25);
			if(rcut<1e3) f *= exp(-1e0*pow_2(Lc)/(MBH*rcut));
			//if(std::isinf(f) || f<=0) printf("(%g %g %g %g)",Lc,E,E0s,f);
			return f;
		}
};

/*
 *DFs are normalised to int d^3J f = 1
*/
class EXP HernquistDF: public ergodicDF {
	private:
		const double mass, scaleRadius, factor;
	public:
		HernquistDF(double _mass, double _scaleRadius):
		    mass(_mass), scaleRadius(_scaleRadius),
		    factor(sqrt(2)*pow_3(2*M_PI*sqrt(mass*scaleRadius))) {}
		virtual double value(const double E) const;
		double M() const{ return mass;}
};
class EXP BaseIsochroneDF: public ergodicDF {
	public:
		BaseIsochroneDF(){}
		virtual double value(const double E) const=0;
		virtual double M() const = 0;
};
class EXP IsochroneDF_E: public BaseIsochroneDF {
	private:
		const double mass, scaleRadius, factor;
	public:
		IsochroneDF_E(double _mass, double _scaleRadius):
		    mass(_mass), scaleRadius(_scaleRadius),
		    factor(sqrt(2)*pow_3(2*M_PI*sqrt(mass*scaleRadius))) {}
		virtual double value(const double E) const;
		virtual double M() const{ return mass;}
};
class EXP IsochroneDF_Lc: public BaseIsochroneDF {
	private:
		const double mass, scaleRadius, factor;
		const math::LinearInterpolator& LcE;
	public:
		IsochroneDF_Lc(double _mass, double _scaleRadius,
			    const math::LinearInterpolator& _LcE):
		    mass(_mass), scaleRadius(_scaleRadius), LcE(_LcE),
		    factor(sqrt(2)*pow_3(2*M_PI*sqrt(mass*scaleRadius))) {}
		virtual double value(const double E) const;
		virtual double M() const{ return mass;}
};

}//namespace

