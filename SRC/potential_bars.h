#include "potential_base.h"

namespace potential{
class EXP bar_potS: public potential::BasePotential{
	private:
		const double A,V0e,rq;
		double Phi2(const double&) const;
		double dPhi2(const double&) const;
		double F(const double&) const;
		double dF(const double&) const;
	public:
		bar_potS(const double _A, const double _V0e, const double _rq) :
		    A(_A), V0e(_V0e), rq(_rq) {}
		virtual void evalCyl(const coord::PosCyl &pos, double* potential=NULL,
				     coord::GradCyl* deriv=NULL, coord::HessCyl* deriv2=NULL) const;
		virtual void evalCar(const coord::PosCar& p, double* Phi, coord::GradCar* deriv, coord::HessCar* deriv2) const {
			coord::GradCyl derivR;
			coord::PosCyl Rz(coord::toPosCyl(p));
			evalCyl(Rz,Phi,&derivR,NULL);
			deriv->dx=derivR.dR*cos(Rz.phi)-derivR.dphi*sin(Rz.phi);
			deriv->dy=derivR.dR*sin(Rz.phi)+derivR.dphi*cos(Rz.phi);
			deriv->dz=derivR.dz;
		}
		virtual void evalSph(const coord::PosSph &p, double* Phi, coord::GradSph* deriv, coord::HessSph* deriv2) const {
			coord::GradCyl derivR;
			coord::PosCyl Rz(coord::toPosCyl(p));
			evalCyl(Rz,Phi,&derivR,NULL);
			deriv->dr=derivR.dR*sin(p.theta)+derivR.dz*cos(p.theta);
			deriv->dtheta=derivR.dR*cos(p.theta)-derivR.dz*sin(p.theta);
			deriv->dphi=derivR.dphi;
		}
		double rho(const coord::PosCyl&) const;
		virtual coord::SymmetryType symmetry() const { return coord::ST_REFLECTION; }
		virtual const char* name() const { return "SormaniBar"; }
};
}//namespace