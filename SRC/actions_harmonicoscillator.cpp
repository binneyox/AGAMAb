#include "actions_harmonicoscillator.h"

namespace actions {
	double HarmonicOscilattor::H(coord::PosMomCyl& rp) const {
		double H = .5 * (pow_2(rp.pR) + pow_2(rp.pz))
			+ .5 * (pow_2(omegaR * rp.R) + pow_2(omegaz * rp.z));
		return rp.pphi == 0 ? H : H + .5 * pow_2(rp.pphi / (rp.R));
	}
	Actions HarmonicOscilattor::pq2J(const coord::PosMomCyl& pc) const {
		double omegaR2 = omegaR * omegaR, omegaz2 = omegaz * omegaz;
		double Jz = .5 * (pow_2(pc.pz) + pow_2(omegaz * pc.z)) / omegaz;
		double HR = .5 * (pow_2(pc.pR))
			+ .5 * (pow_2(omegaR * pc.R));
		if (pc.pphi != 0) HR += .5 * pow_2(pc.pphi / (pc.R));
		double Jr = .5 * (HR / omegaR - fabs(pc.pphi));
		return Actions(Jr, Jz, pc.pphi);
	}
	coord::PosMomCyl HarmonicOscilattor::aa2pq(const ActionAngles& aa, Frequencies* freqs,
		DerivAct<coord::Cyl>* dJ, DerivAng<coord::Cyl>* dA) const {
		//freqz=omegaz, freqr=2*omega
		double z = sqrt(2 * aa.Jz / omegaz) * sin(aa.thetaz);
		double pz = sqrt(2 * aa.Jz * omegaz) * cos(aa.thetaz);
		double J1 = 2 * aa.Jr + fabs(aa.Jphi);
		double e = 2. * sqrt(aa.Jr * (aa.Jr + fabs(aa.Jphi))) / J1;
		double R = aa.Jphi == 0 ? sqrt(2 * J1 / omegaR) * (sin(aa.thetar / 2)) : sqrt(J1 / omegaR * (1 - e * cos(aa.thetar)));
		double pR = aa.Jphi == 0 ? sqrt(2 * J1 * omegaR) * cos(aa.thetar / 2) : J1 * e * sin(aa.thetar) / R;
		double K = (aa.Jphi == 0) ? 0 : sqrt((1 + e) / (1 - e));
		double phi = (aa.Jphi == 0) ? 0 : aa.thetaphi - .5 * aa.thetar +
			math::sign(aa.Jphi) * atan(K * tan(aa.thetar / 2));
		if (aa.Jphi == 0 && aa.thetar > M_PI)phi += M_PI;
		coord::PosMomCyl xpcyl(R, z, phi, pR, pz, aa.Jphi);
		if (freqs)*freqs = Frequencies(2. * omegaR, omegaz, omegaR);
		if (dJ || dA) {
			double dedJr = pow_2(aa.Jphi) / (sqrt(aa.Jr * (aa.Jr + fabs(aa.Jphi))) * J1 * J1);
			double dedJphi = -aa.Jphi * sqrt(aa.Jr) / (sqrt(aa.Jr + fabs(aa.Jphi)) * pow_2(J1));
			double dzdJz = .5 * z / aa.Jz;double dpzdJz = .5 * pz / aa.Jz;
			double dRdJr = aa.Jphi == 0 ? R / J1 : R / J1 - .5 * J1 / omegaR * dedJr * cos(aa.thetar) / R;
			double dpRdJr = aa.Jphi == 0 ? pR / J1 : 2 * pR / J1 + dedJr * J1 * sin(aa.thetar) / R - pR / R * dRdJr;
			double dRdJphi = aa.Jphi == 0 ? .5 * R / J1 : .5 * R / J1 - .5 * J1 / omegaR * dedJphi * cos(aa.thetar) / R;
			double dpRdJphi = aa.Jphi == 0 ? .5 * pR / J1 :
				pR / J1 + dedJphi * J1 * sin(aa.thetar) / R - pR / R * dRdJphi;
			double dzdthetaz = pz / omegaz;
			double dRdthetar = pR / (2 * omegaR);
			double K1 = -pow_2(omegaR) * R;
			K1 += aa.Jphi == 0 ? 0 : pow_2(aa.Jphi / R) / R;
			double dpRdthetar = K1 / (2 * omegaR);
			double dpzdthetaz = -omegaz * z;

			double dphidJr = (aa.Jphi == 0) ? 0 : math::sign(aa.Jphi) * tan(aa.thetar / 2) / (1 + pow_2(K * tan(aa.thetar / 2))) * dedJr / (sqrt(1 - e * e) * (1 - e));
			double dphidJphi = (aa.Jphi == 0) ? 0 : math::sign(aa.Jphi) * tan(aa.thetar / 2) / (1 + pow_2(K * tan(aa.thetar / 2))) * dedJphi / (sqrt(1 - e * e) * (1 - e));
			double dphidthetar = (aa.Jphi == 0) ? 0 : -.5 + .5 * math::sign(aa.Jphi) / (1 + pow_2(K * tan(aa.thetar / 2))) * K / pow_2(cos(aa.thetar / 2));
			if (dJ) {
				dJ->dbyJz.z = dzdJz;
				dJ->dbyJz.R = 0.0;
				dJ->dbyJz.phi = 0.0;
				dJ->dbyJz.pR = 0.0;
				dJ->dbyJz.pz = dpzdJz;
				dJ->dbyJz.pphi = 0.0;

				dJ->dbyJr.R = dRdJr;
				dJ->dbyJr.z = 0.0;
				dJ->dbyJr.phi = dphidJr;
				dJ->dbyJr.pR = dpRdJr;
				dJ->dbyJr.pz = 0.0;
				dJ->dbyJr.pphi = 0.0;

				dJ->dbyJphi.R = dRdJphi;
				dJ->dbyJphi.z = 0.0;
				dJ->dbyJphi.phi = dphidJphi;
				dJ->dbyJphi.pR = dpRdJphi;
				dJ->dbyJphi.pz = 0.0;
				dJ->dbyJphi.pphi = 1.0;
			}
			if (dA) {
				//need to sort out derivatives w.r.t phi for both actions and angles
				dA->dbythetaz.z = dzdthetaz;
				dA->dbythetaz.R = 0.0;
				dA->dbythetaz.phi = 0.0;
				dA->dbythetaz.pR = 0.0;
				dA->dbythetaz.pz = dpzdthetaz;
				dA->dbythetaz.pphi = 0.0;

				dA->dbythetar.R = dRdthetar;
				dA->dbythetar.z = 0.0;
				dA->dbythetar.phi = dphidthetar;//not correct
				dA->dbythetar.pR = dpRdthetar;
				dA->dbythetar.pz = 0.0;
				dA->dbythetar.pphi = 0.0;

				dA->dbythetaphi.R = 0.0;
				dA->dbythetaphi.z = 0.0;
				dA->dbythetaphi.phi = 1.0;
				dA->dbythetaphi.pR = 0.0;
				dA->dbythetaphi.pz = 0.0;
				dA->dbythetaphi.pphi = 0.0;
			}
		}
		return xpcyl;
	}
	coord::PosMomCar HarmonicOscilattor::aa2pqCar(const ActionAngles& aa, Frequencies* freqs,
		DerivAct<coord::Car>* dJ, DerivAng<coord::Car>* dA) const {
		//freqz=omegaz, freqr=2*omega
		double z = sqrt(2 * aa.Jz / omegaz) * sin(aa.thetaz);
		double pz = sqrt(2 * aa.Jz * omegaz) * cos(aa.thetaz);
		double x = sqrt(4 * aa.Jr / omegaR) * sin(aa.thetar / 2);
		double px = sqrt(4 * aa.Jr * omegaR) * cos(aa.thetar / 2);
		coord::PosMomCar xpcar(x, 0.0, z, px, 0.0, pz);
		if (freqs)*freqs = Frequencies(2. * omegaR, omegaz, omegaR);
		if (dJ || dA) {
			int sgn = aa.Jphi != 0 ? math::sign(aa.Jphi) : 0;
			double dzdJz = .5 * z / aa.Jz;double dpzdJz = .5 * pz / aa.Jz;
			double dxdJr = .5 * x / aa.Jr;
			double dpxdJr = .5 * px / aa.Jr;
			double dzdthetaz = pz / omegaz;
			double dxdthetar = px / (2 * omegaR);
			double K1 = -pow_2(omegaR) * x;
			double dpxdthetar = -.5 * omegaR * x;
			double dpzdthetaz = -omegaz * z;
			if (dJ) {
				dJ->dbyJz.z = dzdJz;
				dJ->dbyJz.x = 0.0;
				dJ->dbyJz.y = 0.0;
				dJ->dbyJz.px = 0.0;
				dJ->dbyJz.pz = dpzdJz;
				dJ->dbyJz.py = 0.0;

				dJ->dbyJr.x = dxdJr;
				dJ->dbyJr.z = 0.0;
				dJ->dbyJr.y = 0.0;//phi angles not done correctly
				dJ->dbyJr.px = dpxdJr;
				dJ->dbyJr.pz = 0.0;
				dJ->dbyJr.py = 0.0;

				dJ->dbyJphi.x = 0.0;
				dJ->dbyJphi.z = 0.0;
				dJ->dbyJphi.y = 0.0;
				dJ->dbyJphi.px = 0.0;
				dJ->dbyJphi.pz = 0.0;
				dJ->dbyJphi.py = 0.0;
			}
			if (dA) {
				//need to sort out derivatives w.r.t phi for both actions and angles
				dA->dbythetaz.z = dzdthetaz;
				dA->dbythetaz.x = 0.0;
				dA->dbythetaz.y = 0.0;
				dA->dbythetaz.px = 0.0;
				dA->dbythetaz.pz = dpzdthetaz;
				dA->dbythetaz.py = 0.0;

				dA->dbythetar.x = dxdthetar;
				dA->dbythetar.z = 0.0;
				dA->dbythetar.y = 0.0;//not correct
				dA->dbythetar.px = dpxdthetar;
				dA->dbythetar.pz = 0.0;
				dA->dbythetar.py = 0.0;

				dA->dbythetaphi.x = 0.0;
				dA->dbythetaphi.z = 0.0;
				dA->dbythetaphi.y = 0.0;
				dA->dbythetaphi.px = 0.0;
				dA->dbythetaphi.pz = 0.0;
				dA->dbythetaphi.pz = 0.0;
			}
		}
		return xpcar;
	}
	ActionAngles HarmonicOscilattor::pq2aa(const coord::PosMomCyl& pc, Frequencies* freqs) const {
		double omegaR2 = omegaR * omegaR, omegaz2 = omegaz * omegaz;
		double Jz = .5 *  (pow_2(pc.pz) + pow_2(omegaz * pc.z))/omegaz;
		double HR = .5 * (pow_2(pc.pR))
			+ .5 * (pow_2(omegaR * pc.R));
		if (pc.pphi != 0) HR += .5 * pow_2(pc.pphi / (pc.R));
		double Jr = .5 * (HR / omegaR - fabs(pc.pphi));
		double J1 = (2 * Jr + fabs(pc.pphi));
		double e = 2. * sqrt(Jr * (Jr + fabs(pc.pphi))) / J1;
		double thetaz = atan2(omegaz * pc.z, pc.pz);
		double thetar = atan2(pc.pR * pc.R, (J1 - pow_2(pc.R) * omegaR));
		double K = (pc.pphi == 0) ? 0 : sqrt((1 + e) / (1 - e));
		double thetaphi = 0;
		if (pc.pphi != 0) {
			if (thetar == .5 * M_PI||thetar==-.5*M_PI)thetaphi= pc.phi + .5 * thetar -
				math::sign(pc.pphi) * math::sign(thetar)*.5*M_PI;
			else thetaphi = pc.phi + .5 * thetar -math::sign(pc.pphi) * atan(K * tan(thetar / 2));
		}
		ActionAngles aa(Actions(Jr, Jz, pc.pphi), Angles(thetar, thetaz, thetaphi));
		if (freqs)*freqs = Frequencies(2. * omegaR, omegaz, omegaR * math::sign(pc.pphi));
		return aa;
	}
	coord::PosMomCyl HarmonicOscilattor::aa2pq(const ActionAngles& aa,
		coord::PosMomCyl& drdomegar, coord::PosMomCyl& drdomegaz) const {
		double z = sqrt(2 * aa.Jz / omegaz) * sin(aa.thetaz);
		double pz = sqrt(2 * aa.Jz * omegaz) * cos(aa.thetaz);
		double J1 = 2 * aa.Jr + fabs(aa.Jphi);
		double e = 2. * sqrt(aa.Jr * (aa.Jr + fabs(aa.Jphi))) / J1;
		double R = aa.Jphi == 0 ? sqrt(2 * J1 / omegaR) * (sin(aa.thetar / 2)) : sqrt(J1 / omegaR * (1 - e * cos(aa.thetar)));
		double pR = aa.Jphi == 0 ? sqrt(2 * J1 * omegaR) * cos(aa.thetar / 2) : J1 * e * sin(aa.thetar) / R;
		double K = sqrt((1 + e) / (1 - e));
		double phi = (aa.Jphi == 0) ? 0 : aa.thetaphi - .5 * aa.thetar +
			math::sign(aa.Jphi) * atan(K * tan(aa.thetar / 2));
		if (aa.Jphi == 0 && aa.thetar > M_PI)phi += M_PI;
		coord::PosMomCyl xpcyl(R, z, phi, pR, pz, aa.Jphi);
		double dzdomegaz = -.5 * z / omegaz;
		double dpzdomegaz = .5 * pz / omegaz;
		double dRdomegaR = -.5 * R / omegaR;
		double dpRdomegaR = .5 * pR / omegaR;
		double dphidomegaR = 0.0;
		drdomegar = coord::PosMomCyl(dRdomegaR, 0, 0.0, dpRdomegaR, 0.0, 0.0);
		drdomegaz = coord::PosMomCyl(0, dzdomegaz, 0.0, 0.0, dpzdomegaz, 0.0);
		return xpcyl;
	}
	coord::PosMomCar HarmonicOscilattor::aa2pq(const ActionAngles& aa,
		coord::PosMomCar& drdomegar, coord::PosMomCar& drdomegaz) const {

		double z = sqrt(2 * aa.Jz / omegaz) * sin(aa.thetaz);
		double pz = sqrt(2 * aa.Jz * omegaz) * cos(aa.thetaz);
		double x = sqrt(4 * aa.Jr / omegaR) * sin(aa.thetar / 2);
		double px = sqrt(4 * aa.Jr * omegaR) * cos(aa.thetar / 2);
		coord::PosMomCar xv(x, 0, z, px, 0, pz);
		double dzdomegaz = -.5 * z / omegaz;
		double dpzdomegaz = .5 * pz / omegaz;
		double dxdomegaR = -.5 * x / omegaR;
		double dpxdomegaR = .5 * px / omegaR;
		drdomegar = coord::PosMomCar(dxdomegaR, 0, 0.0, dpxdomegaR, 0.0, 0.0);
		drdomegaz = coord::PosMomCar(0, 0, dzdomegaz, 0.0, 0.0, dpzdomegaz);
		return xv;
	}
	//for Jphi=0
	coord::PosMomCar HarmonicOscilattor::aa2pqCar(const ActionAngles& aa,
		coord::PosMomCar& drdomegar, coord::PosMomCar& drdomegaz) const {

		double z = sqrt(2 * aa.Jz / omegaz) * sin(aa.thetaz);
		double pz = sqrt(2 * aa.Jz * omegaz) * cos(aa.thetaz);
		double J1 = 2 * aa.Jr + fabs(aa.Jphi);
		double e = 2. * sqrt(aa.Jr * (aa.Jr + fabs(aa.Jphi))) / J1;
		double x = sqrt(4 * aa.Jr / omegaR) * sin(aa.thetar / 2);
		double px = sqrt(4 * aa.Jr * omegaR) * cos(aa.thetar / 2);
		coord::PosMomCar xpcar(x, z, 0.0, px, pz, 9.9);
		double dzdomegaz = -.5 * z / omegaz;
		double dpzdomegaz = .5 * pz / omegaz;
		double dRdomegaR = -.5 * x / omegaR;
		double dpRdomegaR = .5 * px / omegaR;
		double dphidomegaR = 0.0;
		drdomegar = coord::PosMomCar(dRdomegaR, 0, 0.0, dpRdomegaR, 0.0, 0.0);
		drdomegaz = coord::PosMomCar(0, dzdomegaz, 0.0, 0.0, dpzdomegaz, 0.0);
		return  xpcar;

	}
	HarmonicOscilattor interpHarmonicOscillator(const double x, const HarmonicOscilattor& os0, const HarmonicOscilattor& os1) {
		const double xp = 1 - x;
		return HarmonicOscilattor(os0.omegaR * x + os1.omegaR * xp, os0.omegaz * x + os1.omegaR * xp);
	}
	//*/

}//namespace actions