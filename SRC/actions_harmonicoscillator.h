#pragma once
#include "potential_base.h"
#include "potential_utils.h"
#include "potential_factory.h"
#include "actions_base.h"
namespace actions {
	class EXP HarmonicOscilator {

	public:
		double omegaR, omegaz;
		HarmonicOscilator(double _omegaR = 1, double _omegaz = 1) : omegaR(fabs(_omegaR)), omegaz(fabs(_omegaz)) {};
		double H(coord::PosMomCyl& rp) const;
		double H(Actions J) const {
			return ((2 * J.Jr + fabs(J.Jphi)) * omegaR + J.Jz*omegaz);
		}
		double Omegar() const {
			return omegaR;
		}
		double Omegaz() const {
			return omegaz;
		}
		double PE(coord::PosMomCyl& rp) const { return .5 * (pow_2(omegaR * rp.R) + pow_2(omegaz * rp.z)); }
		ActionAngles pq2aa(const coord::PosMomCyl& Rz, Frequencies* freqs = NULL) const;
		Actions pq2J(const coord::PosMomCyl& rtheta) const;
		coord::PosMomCyl aa2pq(const ActionAngles& aa,
			coord::PosMomCyl& drdomegar, coord::PosMomCyl& drdomegaz) const;
		coord::PosMomCyl aa2pq(const ActionAngles& aa, Frequencies* freqs=NULL,
			DerivAct<coord::Cyl>* dJ=NULL, DerivAng<coord::Cyl>* dA=NULL) const;
		coord::PosMomCar aa2pqCar(const ActionAngles& aa, Frequencies* freqs = NULL,
			DerivAct<coord::Car>* dJ = NULL, DerivAng<coord::Car>* dA = NULL) const;
		coord::PosMomCar aa2pqCar(const ActionAngles& aa,
			coord::PosMomCar& drdomegar, coord::PosMomCar& drdomegaz) const;
		coord::PosMomCar aa2pq(const ActionAngles& aa,
			coord::PosMomCar& drdomegar, coord::PosMomCar& drdomegaz) const;
		HarmonicOscilator& operator *= (const double a) {
			omegaR *= a; omegaz *= a;
			return *this;
		}
		HarmonicOscilator& operator += (const HarmonicOscilator& os) {
			omegaR += os.omegaR; omegaz += os.omegaz;
			return *this;
		}
		const HarmonicOscilator operator * (const double a) const {
			HarmonicOscilator I2(omegaR * a, omegaz * a);
			return I2;
		}
		const HarmonicOscilator operator + (const HarmonicOscilator ho) const {
			HarmonicOscilator ho2(omegaR + ho.omegaR, ho.omegaz + omegaz);
			return ho2;
		}
	};
	EXP HarmonicOscilator interpHarmonicOscilator(const double, const HarmonicOscilator&, const HarmonicOscilator&);
}
