#include "potential_bars.h"
#include "math_specfunc.h"

namespace potential{


void bar_potS::evalCyl(const coord::PosCyl& Rz, double* Phi, coord::GradCyl* deriv, coord::HessCyl* Hs) const{
	const double r=sqrt(Rz.R*Rz.R+Rz.z*Rz.z), x=Rz.R/r, y=Rz.z/r;
	if(Phi) *Phi = -F(r/rq)*A*pow(V0e*Rz.R/r,2)*cos(2*Rz.phi);
	if(deriv){
		deriv->dR = x*(dPhi2(r)*x*x+2*Phi2(r)/r*(1-x*x))*cos(2*Rz.phi);
		deriv->dz = y*x*x*(dPhi2(r)-2*Phi2(r)/r)*cos(2*Rz.phi);
		deriv->dphi= 2*F(r/rq)*A*pow(V0e*Rz.R/r,2)*sin(2*Rz.phi);
	}
}
double bar_potS::rho(const coord::PosCyl& Rz) const{
	const double r=sqrt(Rz.R*Rz.R+Rz.z*Rz.z);
	return A*pow(V0e/rq*Rz.R/r,2)*exp(-2*r/rq)*cos(2*Rz.phi);
}
double bar_potS::F(const double& x) const{
	const double e2x=exp(2*x),TINY=.05;
	if(x>TINY)
		return 0.05*(3-(3+x*(6+x*(6+x*(4+x*2))))/e2x)/pow(x,3)+.2*x*x*math::E1(2*x);
	else
		return x*x*(.8/(20*e2x)+.2*math::E1(2*x));
}
double bar_potS::dF(const double& x) const {
	const double e2x=exp(2*x),TINY=.05;
	if(x>TINY)
		return .1*(3+x*(6+x*(6+x*(4+2*x))))/e2x/pow(x,3)-.2*x/e2x
				-(9-(9+x*(12+x*(6-2*x*x)))/e2x)/(20*x*x*x*x)
				+.4*x*math::E1(2*x);
	else
		return 2*x*(.8/(20*e2x)+.2*math::E1(2*x))
				-x*x*(1.6/(20*e2x)+.2/(x*e2x));
}
double bar_potS::Phi2(const double& r) const {
	return -F(r/rq)*A*pow(V0e,2);
}
double bar_potS::dPhi2(const double& r) const{
	return -dF(r/rq)/rq*A*pow(V0e,2);
}	

}//namespace