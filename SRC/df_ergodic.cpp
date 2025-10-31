#include "df_ergodic.h"

namespace{
/* Helper class for evaluation of density and frictional drag
 */
class rho_integrand: public math::IFunctionNoDeriv {
	private:
		const df::ergodicDF* DF;
		const double Phi;
	public:
		rho_integrand(const df::ergodicDF* _DF,double _Phi): DF(_DF), Phi(_Phi){}
		virtual double value(const double v) const{
			double E = fmin(0,.5*v*v + Phi);//shun positive E
			return v*v*(*DF)(E);
		}
};

}//anon namespace

namespace df{

//fraction of density contributed as  v<Vmax
double ergodicDF::fraction(const double Phi,const double Vmax) const{
	rho_integrand rhoI(this,Phi);
	return 4*M_PI*math::integrate(rhoI,0,Vmax,1e-5);
}

double HernquistDF::value(const double E) const{
	double cE = -E*scaleRadius/mass;
	if(cE<=0) return 0;
	double sqE = sqrt(cE), cEp = 1-cE;
	double f = sqE/pow_2(cEp)*((2*cE-1)*(8*cE*cEp+3)
				   +3*asin(sqE)/sqrt(cE*cEp));
	return f/factor;
}

double IsochroneDF::value(const double E) const{
	double cE = -E*scaleRadius/mass;
	double sqE = sqrt(cE), cEp = 1-cE;
	double f = sqE/pow_4(2*cEp)*(27-(66-(320-(240-64*cE)*cE)*cE)*cE
				   +3*(cE*(cE*16+28)-9)*asin(sqE)/sqrt(cE*cEp));
	return f/factor;
}

}//namespace df
