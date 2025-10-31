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
		virtual double fraction(const double Phi,const double Vmax) const;
};
class EXP HernquistDF: public ergodicDF {
	private:
		const double mass, scaleRadius, factor;
	public:
		HernquistDF(double _mass, double _scaleRadius):
		    mass(_mass), scaleRadius(_scaleRadius),
		    factor(sqrt(2)*pow_3(2*M_PI*sqrt(mass*scaleRadius))) {}
		virtual double value(const double E) const;
};
class EXP IsochroneDF: public ergodicDF {
	private:
		const double mass, scaleRadius, factor;
	public:
		IsochroneDF(double _mass, double _scaleRadius):
		    mass(_mass), scaleRadius(_scaleRadius),
		    factor(sqrt(2)*pow_3(2*M_PI*sqrt(mass*scaleRadius))) {}
		virtual double value(const double E) const;
};

}//namespace

