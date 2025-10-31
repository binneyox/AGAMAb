/* AngleAction coords of the isochrone sphere. Provides H, PE Omegar and OmegaL
 * in addition to p22J, pq2aa and aa2pq.
 * aa2pq optionally computes derivatives wrt (i) J, (ii) theta, (iii) M and b
*/
#pragma once
#include "math_core.h"
#include "math_specfunc.h"
#include "actions_base.h"
#include "coord.h"

namespace actions {

class EXP Isochrone {
	public:
		double Js, b;
		Isochrone(double _Js=1, double _b=1) :  Js(_Js), b(_b) {
			if(Js<=0 || b<=0){
				printf("Illegal Isochrone parameter: %g %g\n",Js,b);
			}
		}
		void reset(double _Js,double _b){
			Js=_Js; b=_b;
			if(Js<=0 || b<=0){
				printf("Illegal Isochrone parameter: %g %g\n",Js,b);
			}
		}
		double H(coord::PosMomSph& rp) const;
		double H(Actions J) const{
			const double L=J.Jz+fabs(J.Jphi);
			return -.5*pow_2(pow_2(Js)/b)/pow_2(J.Jr+.5*(L+sqrt(L*L+4*Js*Js)));
		}
		double Omegar(const double Jr,const double L) const{
			return pow(Js,4)/(pow_2(b)*pow(Jr+.5*(L+sqrt(L+L+4*Js*Js)),3));
		}
		double OmegaL(const double Jr,const double L) const{
			return .5*(1+L/sqrt(L*L+4*Js*Js))*Omegar(Jr,L);
		}
		double PE(coord::PosMomSph& rp) const;
		coord::PosMomSph aa2pq(const ActionAngles& aa, Frequencies* freqs=NULL,
				       DerivAct<coord::Sph>* dJ=NULL, DerivAng<coord::Sph>* dA=NULL) const;
		coord::PosMomSph aa2pq(const Actions& J, const Angles& theta, Frequencies* freqs=NULL,
				       DerivAct<coord::Sph>* dJ=NULL, DerivAng<coord::Sph>* dA=NULL) const{
			return aa2pq(ActionAngles(J,theta), freqs, dJ, dA);
		}
		coord::PosMomSph aa2pq(const ActionAngles& aa,
				       coord::PosMomSph& drdM, coord::PosMomSph& drdb) const;
		coord::PosMomSph aa2pq(const Actions& J, const Angles& theta,
				       coord::PosMomSph& drdM, coord::PosMomSph& drdb) const{
			return  aa2pq(ActionAngles(J,theta), drdM, drdb);
		}
		ActionAngles pq2aa(const coord::PosMomSph& rtheta, Frequencies* freqs=NULL) const;
		Actions pq2J(const coord::PosMomSph& rtheta) const;
		Isochrone& operator *= (const double a) {
			Js *= a; b *= a;
			return *this;
		}
		Isochrone& operator += (const Isochrone& I){
			Js += I.Js; b+= I.b;
			return *this;
		}
		const Isochrone operator * (const double a) const {
			Isochrone I2(Js*a,b*a);
			return I2;
		}
		const Isochrone operator + (const Isochrone I) const{
			Isochrone I2(Js+I.Js,b+I.b);
			return I2;
		}
};
EXP Isochrone interpIsochrone(const double, const Isochrone&, const Isochrone&);

}//namespace