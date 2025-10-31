#include "actions_newgenfnc.h"
#include "math_fit.h"
#include "math_core.h"
#include <cmath>

#define bmax .95

namespace actions {

	namespace { // internal

		/// return the absolute value of an element in a map, or zero if it doesn't exist
		static inline double absvalue(const std::map< std::pair<int, int>, double >& indPairs, int ir, int iz)
		{
			if (indPairs.find(std::make_pair(ir, iz)) != indPairs.end())
				return fabs(indPairs.find(std::make_pair(ir, iz))->second);
			else
				return 0;
		}

		/// return the value of an element in a map, or zero if it doesn't exist
		static inline double value(const std::map< std::pair<int, int>, double >& indPairs, int ir, int iz)
		{
			if (indPairs.find(std::make_pair(ir, iz)) != indPairs.end())
				return (indPairs.find(std::make_pair(ir, iz))->second);
			else
				return 0;
		}

		/// Does this series already have a fraction?
		static bool haveFrac(const int iz, GenFncFitFracs& fracs) {// true if we already have iz frac
			bool need = true;
			for (unsigned int j = 0; j < fracs.size(); j++)
				need &= (fracs[j].mz != iz);
			return !need;
		}

		/// Should we add a fraction to this series?
		static bool addFrac(std::map< std::pair<int, int>, double>& indPairs,
			const int maxmr, const int iz, const int nrStart, const int nrBack,
			double& B, double& bmean) {
			unsigned int mr_last = 0;
			for (int ir = maxmr; ir > 2; ir--) {//find highest term in this series
				if (indPairs.find(std::make_pair(ir, iz)) != indPairs.end()) {
					mr_last = ir; break;
				}
			}
			if (mr_last < nrStart) return false; // too soon to introduce a frac
			std::vector<double> b;
			for (int j = 1; j < nrBack + 1; j++) {
				b.push_back(value(indPairs, mr_last - j, iz) / value(indPairs, mr_last - j - 1, iz));
				if (std::isnan(b.back()) || fabs(b.back()) > 0.9) return false; // still rising?
			}
			double bvar = 0; bmean = 0;
			for (int i = 0;i < nrBack;i++) {
				bvar += pow_2(b[i]); bmean += b[i];
			}
			bmean /= (double)nrBack;
			if (bmean > 0.95) return false;
			bvar = sqrt(bvar / (double)nrBack - pow_2(bmean));
			if (bvar > 0.15) return false;
			B = value(indPairs, mr_last - nrBack, iz) / pow(bmean, mr_last - nrBack - 1);
			return true;
		}

		/// create an array of angles uniformly covering the range [0:pi]  (NB: why not 2pi?)
		static std::vector<Angles> makeGridAngles(unsigned int nr, unsigned int nz, unsigned int nphi = 1)
		{
			std::vector<Angles> vec(nr * nz * nphi);
			for (unsigned int ir = 0; ir < nr; ir++) {
				double thetar = ir * M_PI / nr;
				for (unsigned int iz = 0; iz < nz; iz++) {
					double thetaz = iz * M_PI / nz;
					for (unsigned int iphi = 0; iphi < nphi; iphi++)
						vec[(ir * nz + iz) * nphi + iphi] =
						Angles(thetar, thetaz, iphi * M_PI / fmax(nphi - 1, 1));
				}
			}
			return vec;
		}

		/// create grid in angles with size determined by the maximal Fourier harmonic in the indices array
		static std::vector<Angles> makeGridAngles(const GenFncIndices& indices,
			const int timesr, const int timesz)
		{
			//	int maxmr=4, maxmz=4, maxmphi=0;
			int maxmr = 0, maxmz = 0, maxmphi = 0;
			for (unsigned int i = 0; i < indices.size(); i++) {
				maxmr = std::max<int>(maxmr, math::abs(indices[i].mr));
				maxmz = std::max<int>(maxmz, math::abs(indices[i].mz));
				maxmphi = std::max<int>(maxmphi, math::abs(indices[i].mphi));
			}
			maxmz /= 2;//because GF ~ sin(2*n*theta_z)
			//	return makeGridAngles(8*(maxmr/4+1), 6*(maxmz/4+1), maxmphi>0 ? 6*(maxmphi/4+1) : 1);
			//	return makeGridAngles(maxmr>0? 6*(maxmr/4+1) : 2, 6*(maxmz/3+1), maxmphi>0 ? 6*(maxmphi+1) : 1);
			return makeGridAngles(maxmr > 0 ? timesr * (maxmr / 4 + 1) : 2, timesz * (maxmz / 4 + 1), maxmphi > 0 ? 6 * (maxmphi + 1) : 1);
		}

		/** Helper class to be used in the iterative solution of a nonlinear system of equations
			that implicitly define the toy angles as functions of real angles and the derivatives
			of the generating function.
		*/
		class AngleFinder : public math::IFunctionNdimDeriv {
		public:
			AngleFinder(const GenFnc* _GF, const Angles& _ang) : GF(_GF), ang(_ang) {}
			virtual unsigned int numVars() const { return 3; }
			virtual unsigned int numValues() const { return 3; }

			virtual void evalDeriv(const double vars[], double values[], double* derivs = 0) const
			{
				if (values) {
					values[0] = vars[0] - ang.thetar;
					values[1] = vars[1] - ang.thetaz;
					values[2] = vars[2] - ang.thetaphi;
				}
				if (derivs) {
					derivs[0] = derivs[4] = derivs[8] = 1.;  // diagonal
					derivs[1] = derivs[2] = derivs[3] = derivs[5] = derivs[6] = derivs[7] = 0;  // off-diag
				}
				for (unsigned int i = 0; i < GF->indices.size(); i++) {
					double arg = GF->indices[i].mr * vars[0] + GF->indices[i].mz * vars[1] +
						GF->indices[i].mphi * vars[2];    // argument of trig functions
					if (values) {
						double s = sin(arg);
						values[0] += s * GF->derivs[i].Jr;
						values[1] += s * GF->derivs[i].Jz;
						values[2] += s * GF->derivs[i].Jphi;
					}
					if (derivs) {
						double c = cos(arg);
						derivs[0] += c * GF->derivs[i].Jr * GF->indices[i].mr;
						derivs[1] += c * GF->derivs[i].Jr * GF->indices[i].mz;
						derivs[2] += c * GF->derivs[i].Jr * GF->indices[i].mphi;
						derivs[3] += c * GF->derivs[i].Jz * GF->indices[i].mr;
						derivs[4] += c * GF->derivs[i].Jz * GF->indices[i].mz;
						derivs[5] += c * GF->derivs[i].Jz * GF->indices[i].mphi;
						derivs[6] += c * GF->derivs[i].Jphi * GF->indices[i].mr;
						derivs[7] += c * GF->derivs[i].Jphi * GF->indices[i].mz;
						derivs[8] += c * GF->derivs[i].Jphi * GF->indices[i].mphi;
					}
				}
				for (unsigned int j = 0; j < GF->fracs.size(); j++) {
					Angles thetaT(vars[0], vars[1], vars[2]);
					if (values) {
						double F, BdFdpsi;  GF->fracs[j].dBFdP(thetaT, F, BdFdpsi);
						values[0] += GF->fracs[j].dBdJ.Jr * F + GF->fracs[j].dpsidJ.Jr * BdFdpsi;
						values[1] += GF->fracs[j].dBdJ.Jz * F + GF->fracs[j].dpsidJ.Jz * BdFdpsi;
						values[2] += GF->fracs[j].dBdJ.Jphi * F + GF->fracs[j].dpsidJ.Jphi * BdFdpsi;
					}
					if (derivs) {//Note that b is replaced by psi in dSby and dJTb
						Actions dJTdB, dJTdpsi; GF->fracs[j].derivJ(thetaT, dJTdB, dJTdpsi);
						derivs[0] += GF->fracs[j].dBdJ.Jr * dJTdB.Jr + GF->fracs[j].dpsidJ.Jr * dJTdpsi.Jr;
						derivs[1] += GF->fracs[j].dBdJ.Jr * dJTdB.Jz + GF->fracs[j].dpsidJ.Jr * dJTdpsi.Jz;
						derivs[2] += GF->fracs[j].dBdJ.Jr * dJTdB.Jphi + GF->fracs[j].dpsidJ.Jr * dJTdpsi.Jphi;
						derivs[3] += GF->fracs[j].dBdJ.Jz * dJTdB.Jr + GF->fracs[j].dpsidJ.Jr * dJTdpsi.Jr;
						derivs[4] += GF->fracs[j].dBdJ.Jz * dJTdB.Jz + GF->fracs[j].dpsidJ.Jr * dJTdpsi.Jz;
						derivs[5] += GF->fracs[j].dBdJ.Jz * dJTdB.Jphi + GF->fracs[j].dpsidJ.Jr * dJTdpsi.Jphi;
						derivs[6] += GF->fracs[j].dBdJ.Jphi * dJTdB.Jr + GF->fracs[j].dpsidJ.Jr * dJTdpsi.Jr;
						derivs[7] += GF->fracs[j].dBdJ.Jphi * dJTdB.Jz + GF->fracs[j].dpsidJ.Jr * dJTdpsi.Jz;
						derivs[8] += GF->fracs[j].dBdJ.Jphi * dJTdB.Jphi + GF->fracs[j].dpsidJ.Jr * dJTdpsi.Jphi;
					}
				}
			}
		private:
			const GenFnc* GF;
			const Angles ang;             ///< true angles
		};

	} // internal namespace

	//contribution to JT = B dF/dtheta
	Actions GenFncFitFrac::dJ(const double* p, const Angles& theta) const {
		double B = p[0], b = bmax * tanh(p[1]);
		std::complex<double> tz(exp(std::complex<double>(0, mz * theta.thetaz)));
		std::complex<double> tr(exp(std::complex<double>(0, theta.thetar)));
		std::complex<double> onembtr(std::complex<double>(1, 0) - b * tr);
		return Actions(B * real(tz * tr / pow_2(onembtr)), B * mz * real(tr * tz / onembtr), 0);
	}

	// dJT/dB=dF/dthetaT or dJT/dpsi = dJT/db*(1-b^2) and note dJT/db=Bd^2F/dbdthetaT
	Actions GenFncFitFrac::derivJ(const double* p, bool wrtB, const Angles& theta) const {
		const double B = p[0], b = bmax * tanh(p[1]), sechsq = 1 - pow_2(b / bmax);
		std::complex<double> tz(exp(std::complex<double>(0, mz * theta.thetaz)));
		std::complex<double> tr(exp(std::complex<double>(0, theta.thetar)));
		std::complex<double> onembtr(std::complex<double>(1, 0) - b * tr);
		Actions dJ(0, 0, 0);
		if (wrtB) {		//dif JT wrt B
			dJ.Jz = mz * real(tr * tz / onembtr);
			dJ.Jr = real(tr * tz / pow_2(onembtr));
		}
		else {		//dif JT wrt b
			std::complex<double> Z(tr * tr * tz / pow_2(onembtr));
			dJ.Jz = mz * B * bmax * sechsq * real(Z);
			dJ.Jr = 2 * B * bmax * sechsq * real(Z / onembtr);
		}
		return dJ;
	}

	GenFncFrac::GenFncFrac(GenFncFitFrac& frac, double* params, Actions* derivs) :
		mz(frac.mz), mphi(frac.mphi), B(params[0]), b(bmax* tanh(params[1])),
		sechsq(1 - pow_2(b / bmax)), dBdJ(derivs[0]), dpsidJ(derivs[1]) {
	}
	GenFncFrac::GenFncFrac(double _B, double _b, Actions _dBdJ, Actions _dpsidJ,
		int _mz, int _mphi) :
		B(_B), b(_b), dBdJ(_dBdJ), dpsidJ(_dpsidJ), sechsq(1 - pow_2(b / bmax)), mz(_mz), mphi(_mphi) {
	}

	// contribution to JT
	Actions GenFncFrac::dJ(const Angles& theta) const {//contribution to JT
		std::complex<double> tz(exp(std::complex<double>(0, mz * theta.thetaz)));
		std::complex<double> tr(exp(std::complex<double>(0, theta.thetar)));
		std::complex<double> onembtr(std::complex<double>(1, 0) - b * tr);
		return Actions(B * real(tz * tr / pow_2(onembtr)), B * mz * real(tr * tz / onembtr), 0);
	}

	// dJT/dB=dF/dthetaT or dJT/dpsi = dJT/db*(1-b^2) and note dJT/db=Bd^2F/dbdthetaT
	void GenFncFrac::derivJ(const Angles& theta, Actions& dJTdB, Actions& dJTdpsi) const {
		std::complex<double> tz(exp(std::complex<double>(0, mz * theta.thetaz)));
		std::complex<double> tr(exp(std::complex<double>(0, theta.thetar)));
		std::complex<double> onembtr(std::complex<double>(1, 0) - b * tr);
		std::complex<double> Z(tr * tr * tz / pow_2(onembtr));
		dJTdB.Jz = mz * real(tr * tz / onembtr);
		dJTdB.Jr = real(tr * tz / pow_2(onembtr));
		dJTdB.Jphi = 0;
		dJTdpsi.Jz = mz * B * bmax * sechsq * real(Z);
		dJTdpsi.Jr = 2 * B * bmax * sechsq * real(Z / onembtr);
		dJTdpsi.Jphi = 0;
	}

	//coeffs of dB/dJ, db/dJ in findAngle
	void GenFncFrac::dBFdP(const Angles& theta, double& rF, double& rBdFdpsi) const {
		std::complex<double> tz(exp(std::complex<double>(0, mz * theta.thetaz)));
		std::complex<double> tr(exp(std::complex<double>(0, theta.thetar)));
		std::complex<double> onembtr(std::complex<double>(1, 0) - b * tr);
		std::complex<double> F = tr * tz / onembtr;
		rF = imag(F);
		rBdFdpsi = B * bmax * sechsq * imag(F * tr / onembtr);
	}
	void GenFncFrac::d2dtheta2(const Angles& theta, double& drr, double& dzz, double& dphiphi,
		double& drz, double& drphi, double& dzphi) const {
		std::complex<double> tz(exp(std::complex<double>(0, mz * theta.thetaz)));
		std::complex<double> tr(exp(std::complex<double>(0, theta.thetar)));
		std::complex<double> onembtr(std::complex<double>(1, 0) - b * tr);
		std::complex<double> BF = B * tr * tz / onembtr;
		drr = -imag(BF * (std::complex<double>(1, 0) + b * tr) / pow_2(onembtr));
		dzz = -mz * mz * imag(BF);
		dphiphi = -mphi * mphi * imag(BF);
		drz = -mz * imag(BF / onembtr);
		drphi = -mphi * imag(BF / onembtr);
		dzphi = -mz * mphi * imag(BF);
	}
	int GenFnc::FracNo(const int kz) const {
		for (int i = 0; i < fracs.size(); i++)
			if (fracs[i].mz == kz) return i;
		return -1;
	}
	double GenFnc::strength(void) const {
		std::pair<int, int> maxIs(maxIndices());
		double strngth = 0;
		for (int kz = -maxIs.second; kz <= maxIs.second; kz++) {
			int fracNo = FracNo(kz);
			int N = 0;
			for (int kr = kz > 0 ? 1 : 0; kr <= maxIs.first; kr++) {
				GenFncIndex ind(kr, kz, 0);
				double val(giveValue(ind));
				if (std::isnan(val)) continue;
				N = kr;
				if (fracNo >= 0 && kr > 0)
					val += fracs[fracNo].B * pow(fracs[fracNo].b, kr - 1);
				strngth += val * val;
			}
			if (fracNo >= 0) strngth += pow_2(fracs[fracNo].B *
				pow(fracs[fracNo].b, N)) / (1 - pow_2(fracs[fracNo].b));
		}
		return strngth;
	}

	void GenFnc::write(FILE* ofile) const {
		fprintf(ofile, "%zd %zd\n", indices.size(), fracs.size());
		for (int i = 0; i < indices.size(); i++)
			fprintf(ofile, "%d %d %d %g %g %g %g\n",
				indices[i].mr, indices[i].mz, indices[i].mphi,
				values[i], derivs[i].Jr, derivs[i].Jz, derivs[i].Jphi);
		for (int j = 0; j < fracs.size(); j++) {
			fprintf(ofile, "%d %g %g\n", fracs[j].mz, fracs[j].B, fracs[j].b);
			fprintf(ofile, "%g %g %g %g %g %g\n",
				fracs[j].dBdJ.Jr, fracs[j].dBdJ.Jz, fracs[j].dBdJ.Jphi,
				fracs[j].dpsidJ.Jr, fracs[j].dpsidJ.Jz, fracs[j].dpsidJ.Jphi);
		}
	}
	void GenFnc::read(FILE* ifile) {
		int nt, nf;
		fscanf_s(ifile, "%d %d", &nt, &nf);
		indices.resize(nt); values.resize(nt); derivs.resize(nt + 2 * nf); fracs.resize(nf);
		for (int i = 0; i < indices.size(); i++)
			fscanf_s(ifile, "%d %d %d %lg %lg %lg %lg\n",
				&indices[i].mr, &indices[i].mz, &indices[i].mphi,
				&values[i], &derivs[i].Jr, &derivs[i].Jz, &derivs[i].Jphi);
		for (int j = 0; j < nf; j++) {
			fscanf_s(ifile, "%d %lg %lg", &fracs[j].mz, &fracs[j].B, &fracs[j].b);
			fscanf_s(ifile, "%lg %lg %lg %lg %lg %lg\n",
				&fracs[j].dBdJ.Jr, &fracs[j].dBdJ.Jz, &fracs[j].dBdJ.Jphi,
				&fracs[j].dpsidJ.Jr, &fracs[j].dpsidJ.Jz, &fracs[j].dpsidJ.Jphi);
		}
	}
	void GenFnc::print(void) const {
		for (int i = 0; i < values.size(); i++)
			printf("(%d %d %d) %g  ", indices[i].mr, indices[i].mz, indices[i].mphi, values[i]);
		printf("\n");
	}

	Actions GenFnc::toyJ(const Actions& J, const Angles& thetaT) const {//returns toyJ(thetaT)
		Actions JT(J);
		for (unsigned int i = 0; i < indices.size(); i++) {// compute toy actions
			double val = values[i] * cos(indices[i].mr * thetaT.thetar +
				indices[i].mz * thetaT.thetaz +
				indices[i].mphi * thetaT.thetaphi);
			JT.Jr += val * indices[i].mr;
			JT.Jz += val * indices[i].mz;
			JT.Jphi += val * indices[i].mphi;
		}
		for (unsigned int i = 0; i < fracs.size(); i++) {
			Actions dJ(fracs[i].dJ(thetaT));
			JT.Jr += dJ.Jr;
			JT.Jz += dJ.Jz;
			JT.Jphi += dJ.Jphi;
		}
		if (JT.Jr < 0 || JT.Jz < 0) {// prevent non-physical negative values
			JT.Jr = fmax(JT.Jr, 0); JT.Jz = fmax(JT.Jz, 0);
		}
		return JT;
	}
	Angles GenFnc::trueA(const Angles& thetaT) const {// toy Angs -> real Angs
		Angles theta(thetaT);
		for (unsigned int i = 0; i < indices.size(); i++) {
			double sinT = sin(indices[i].mr * thetaT.thetar + indices[i].mz * thetaT.thetaz
				+ indices[i].mphi * thetaT.thetaphi);
			theta.thetar += derivs[i].Jr * sinT;
			theta.thetaz += derivs[i].Jz * sinT;
			theta.thetaphi += derivs[i].Jphi * sinT;
		}
		for (unsigned int i = 0; i < fracs.size(); i++) {
			double F, BdFdpsi; fracs[i].dBFdP(thetaT, F, BdFdpsi);
			theta.thetar += fracs[i].dBdJ.Jr * F + fracs[i].dpsidJ.Jr * BdFdpsi;
			theta.thetaz += fracs[i].dBdJ.Jz * F + fracs[i].dpsidJ.Jz * BdFdpsi;
			theta.thetaphi += fracs[i].dBdJ.Jphi * F + fracs[i].dpsidJ.Jphi * BdFdpsi;
		}
		return Angles(math::wrapAngle(theta.thetar), math::wrapAngle(theta.thetaz),
			math::wrapAngle(theta.thetaphi));
	}

	Angles GenFnc::toyA(const Angles& theta) const {//real Angs -> toy Angs
		double thetaT[3], Theta[3] = { theta.thetar, theta.thetaz, theta.thetaphi };
		AngleFinder AF(this, theta);
		const int maxIter = 15;
		int numIter = math::findRootNdimDeriv(AF, Theta, 1e-6, maxIter, thetaT);
		//if(numIter>=maxIter) printf("GenFnc toyA: max iterations in findRootNdimDeriv\n");
		return Angles(math::wrapAngle(thetaT[0]),
			math::wrapAngle(thetaT[1]), math::wrapAngle(thetaT[2]));
	}

	ActionAngles GenFnc::true2toy(const ActionAngles& aa) const {
		Angles thetaT(toyA(Angles(aa)));
		Actions JT(toyJ(Actions(aa), thetaT));
		return ActionAngles(JT, thetaT);
	}

	//Derivs of actions wrt thetaT
	DerivAng<coord::Cyl> GenFnc::dJdt(const Angles& thetaT) const {
		DerivAng<coord::Cyl> D;
		for (unsigned int i = 0; i < indices.size(); i++) {
			D.dbythetar.R = 0; D.dbythetaz.R = 0; D.dbythetaphi.R = 0;
			D.dbythetar.z = 0; D.dbythetaz.z = 0; D.dbythetaphi.z = 0;
			D.dbythetar.phi = 0; D.dbythetaz.phi = 0; D.dbythetaphi.phi = 0;
		}
		for (unsigned int i = 0; i < indices.size(); i++) {
			double sinT = values[i] * sin(indices[i].mr * thetaT.thetar
				+ indices[i].mz * thetaT.thetaz
				+ indices[i].mphi * thetaT.thetaphi);
			D.dbythetar.R -= indices[i].mr * indices[i].mr * sinT;//derivs of Jr
			D.dbythetaz.R -= indices[i].mr * indices[i].mz * sinT;
			D.dbythetaphi.R -= indices[i].mr * indices[i].mphi * sinT;
			D.dbythetar.z -= indices[i].mz * indices[i].mr * sinT;//derivs of Jz
			D.dbythetaz.z -= indices[i].mz * indices[i].mz * sinT;
			D.dbythetaphi.z -= indices[i].mz * indices[i].mphi * sinT;
			D.dbythetar.phi -= indices[i].mphi * indices[i].mr * sinT;//derivs of Jphi
			D.dbythetaz.phi -= indices[i].mphi * indices[i].mz * sinT;
			D.dbythetaphi.phi -= indices[i].mphi * indices[i].mphi * sinT;
		}
		for (unsigned int j = 0; j < fracs.size(); j++) {
			double drr, dzz, dphiphi, drz, drphi, dzphi;
			fracs[j].d2dtheta2(thetaT, drr, dzz, dphiphi, drz, drphi, dzphi);
			D.dbythetar.R += drr;
			D.dbythetaz.R += drz;
			D.dbythetaphi.R += drphi;
			D.dbythetar.z += drz;
			D.dbythetaz.z += dzz;
			D.dbythetaphi.z += dzphi;
			D.dbythetar.phi += drphi;
			D.dbythetaz.phi += dzphi;
			D.dbythetaphi.phi += dphiphi;
		}
		return D;
	}

	// dtheta_i/dthetaT_j  * NEES UPDATING FOR FRACs *
	double GenFnc::dtbydtT_Jacobian(const Angles& thetaT, math::Matrix<double>& M) const {
		M(0, 0) = M(1, 1) = M(2, 2) = 1; M(0, 1) = M(0, 2) = M(1, 2) = M(1, 0) = M(2, 0) = M(2, 1) = 0;
		for (unsigned int i = 0; i < indices.size(); i++) {
			double cosT = cos(indices[i].mr * thetaT.thetar
				+ indices[i].mz * thetaT.thetaz
				+ indices[i].mphi * thetaT.thetaphi);
			M(0, 0) += derivs[i].Jr * indices[i].mr * cosT;
			M(0, 1) += derivs[i].Jr * indices[i].mz * cosT;
			M(0, 2) += derivs[i].Jr * indices[i].mphi * cosT;
			M(1, 0) += derivs[i].Jz * indices[i].mr * cosT;
			M(1, 1) += derivs[i].Jz * indices[i].mz * cosT;
			M(1, 2) += derivs[i].Jz * indices[i].mphi * cosT;
			M(2, 0) += derivs[i].Jphi * indices[i].mr * cosT;
			M(2, 1) += derivs[i].Jphi * indices[i].mz * cosT;
			M(2, 2) += derivs[i].Jphi * indices[i].mphi * cosT;
		}
		for (unsigned int j = 0; j < fracs.size(); j++) {
			Actions dBdJ(fracs[j].dBdJ), dpsidJ(fracs[j].dpsidJ);
			Actions dJTdB, dJTdpsi; fracs[j].derivJ(thetaT, dJTdB, dJTdpsi);
			M(0, 0) += dBdJ.Jr * dJTdB.Jr + dpsidJ.Jr * dJTdpsi.Jr;
			M(0, 1) += dBdJ.Jr * dJTdB.Jz + dpsidJ.Jr * dJTdpsi.Jz;
			M(0, 2) += dBdJ.Jr * dJTdB.Jphi + dpsidJ.Jr * dJTdpsi.Jphi;
			M(1, 0) += dBdJ.Jz * dJTdB.Jr + dpsidJ.Jz * dJTdpsi.Jr;
			M(1, 1) += dBdJ.Jz * dJTdB.Jz + dpsidJ.Jz * dJTdpsi.Jz;
			M(1, 2) += dBdJ.Jz * dJTdB.Jphi + dpsidJ.Jz * dJTdpsi.Jphi;
			M(2, 0) += dBdJ.Jphi * dJTdB.Jr + dpsidJ.Jphi * dJTdpsi.Jr;
			M(2, 1) += dBdJ.Jphi * dJTdB.Jz + dpsidJ.Jphi * dJTdpsi.Jz;
			M(2, 2) += dBdJ.Jphi * dJTdB.Jphi + dpsidJ.Jphi * dJTdpsi.Jphi;
		}
		double det = M(0, 0) * (M(1, 1) * M(2, 2) - M(2, 1) * M(1, 2))
			- M(1, 0) * (M(0, 1) * M(2, 2) - M(2, 1) * M(0, 2))
			+ M(2, 0) * (M(0, 1) * M(1, 2) - M(1, 1) * M(0, 2));
		return fabs(det);
	}

	double GenFnc::giveValue(const GenFncIndex& ind) const {
		for (int i = 0; i < indices.size(); i++) {
			if (ind.mr != indices[i].mr) continue;
			if (ind.mz != indices[i].mz) continue;
			if (ind.mphi == indices[i].mphi) return values[i];
		}
		return NAN;
	}

	std::pair<int, int> GenFnc::maxIndices() const {
		int mr = 0, mz = 0;
		for (int i = 0; i < indices.size(); i++) {
			mr = std::max<int>(mr, std::abs(indices[i].mr));
			mz = std::max<int>(mz, std::abs(indices[i].mz));
		}
		return std::make_pair(mr, mz);
	}

	GenFnc interpGenFnc(const double x, const GenFnc& GF1, const GenFnc& GF2) {
		if (x == 1) return GenFnc(GF1);
		if (x == 0) return GenFnc(GF2);
		const double eps = .05, xp = 1 - x;
		GenFnc G2(GF2);
		GenFncIndices indices = GF1.indices;
		std::vector<double> values = GF1.values;
		GenFncDerivs derivs = GF1.derivs;
		for (int i = 0; i < values.size(); i++) {
			values[i] *= x; derivs[i] *= x;
		}
		for (int i = 0; i < indices.size(); i++) {
			int mr = indices[i].mr, mz = indices[i].mz, mphi = indices[i].mphi;
			std::vector<double>::iterator jt = G2.values.begin();
			std::vector<Actions>::iterator kt = G2.derivs.begin();//find terms matching present ones
			for (GenFncIndices::iterator it = G2.indices.begin(); it != G2.indices.end();) {
				if (mr == (*it).mr && mz == (*it).mz && mphi == (*it).mphi) {
					values[i] += xp * (*jt); derivs[i] += (*kt) * xp;
					it = G2.indices.erase(it);
					jt = G2.values.erase(jt);
					kt = G2.derivs.erase(kt);
					break;
				}
				else {
					it++; jt++; kt++;
				}
			}
		}
		if (G2.indices.size() > 0) {
			for (int j = 0; j < G2.indices.size(); j++) {
				indices.push_back(G2.indices[j]);
				values.push_back(G2.values[j]);
				derivs.push_back(G2.derivs[j]);
			}
		}
		GenFncFracs fracs = GF1.fracs;
		for (int j = 0; j < fracs.size(); j++) {
			int mz = fracs[j].mz, mphi = fracs[j].mphi; double b = fracs[j].b;
			fracs[j].B *= x; fracs[j].dBdJ *= x;
			for (GenFncFracs::iterator it = G2.fracs.begin(); it != G2.fracs.end();) {
				if (mz == (*it).mz && mphi == (*it).mphi) {
					if (fabs(b - (*it).b) < eps) {//essentally the same fraction;
						fracs[j].B += xp * (*it).B;
						fracs[j].dBdJ += (*it).dBdJ * xp;
						fracs[j].dpsidJ += (*it).dpsidJ * xp;
						fracs[j].b = x * fracs[j].b + xp * (*it).b;
						fracs[j].sechsq = 1 - pow_2(fracs[j].b / bmax);
						fracs[j].dpsidJ = fracs[j].dpsidJ * x + (*it).dpsidJ * xp;
						it = G2.fracs.erase(it);
					}
					break;
				}
				else
					it++;
			}
		}
		if (G2.fracs.size() > 0) {
			for (int j = 0; j < GF2.fracs.size(); j++)
				fracs.push_back(GF2.fracs[j]);
		}
		return GenFnc(indices, values, derivs, fracs);
	}

	GenFncFitSeries::GenFncFitSeries(const GenFncIndices& _indices, const int timesr,
		const int timesz, const GenFncFitFracs& _fracs,
		const Actions& _acts) :
		indices(_indices), fracs(_fracs), acts(_acts)
	{
		angs = makeGridAngles(indices, timesr, timesz);
		coefs = math::Matrix<double>(angs.size(), indices.size());

		for (unsigned int indexAngle = 0; indexAngle < angs.size(); indexAngle++)
			for (unsigned int indexCoef = 0; indexCoef < indices.size(); indexCoef++)
				coefs(indexAngle, indexCoef) = cos(
					indices[indexCoef].mr * angs[indexAngle].thetar +
					indices[indexCoef].mz * angs[indexAngle].thetaz +
					indices[indexCoef].mphi * angs[indexAngle].thetaphi);
	}

	std::pair<int, int> GenFncFitSeries::maxIndices() const {
		int mr = 0, mz = 0;
		for (int i = 0; i < indices.size(); i++) {
			mr = std::max<int>(mr, std::abs(indices[i].mr));
			mz = std::max<int>(mz, std::abs(indices[i].mz));
		}
		return std::make_pair(mr, mz);
	}


	/* To handle long series in theta_r we add
	 * exp(i*(theta_r + m*theta_z))/(1-b*exp(i*theta_r)) with b=bmax*tanh(psi) to keep |b|<1
	*/
	ActionAngles GenFncFitSeries::toyActionAngles(unsigned int indexAngle, const double values[]) const
	{
		ActionAngles aa(acts, angs[indexAngle]);
		for (unsigned int indexCoef = 0; indexCoef < indices.size(); indexCoef++) {
			double val = values[indexCoef] * coefs(indexAngle, indexCoef);
			aa.Jr += val * indices[indexCoef].mr;
			aa.Jz += val * indices[indexCoef].mz;
			aa.Jphi += val * indices[indexCoef].mphi;
		}
		unsigned int fracNum = 0;
		for (unsigned int indexCoef = indices.size(); indexCoef < numParams(); indexCoef += 2) {
			Actions dJ(fracs[fracNum].dJ(&values[indexCoef], angs[indexAngle]));
			aa.Jr += dJ.Jr;
			aa.Jz += dJ.Jz;
			fracNum++;
		}
		return aa;
	}
	Actions GenFncFitSeries::deriv(const unsigned int indexAngle, const unsigned int indexCoef,
		const double values[]) const {
		Actions dJ;
		if (indexCoef < indices.size()) {//Dif wrt S_k
			double val = coefs(indexAngle, indexCoef);  // no range check performed!
			return Actions(
				val * indices[indexCoef].mr,
				val * indices[indexCoef].mz,
				val * indices[indexCoef].mphi);
		}
		else {//dif wrt parameter B or psi of a fraction
			int fracNum = (indexCoef - indices.size()) / 2, k = (indexCoef - indices.size()) % 2;
			return fracs[fracNum].derivJ
			(&values[indices.size() + 2 * fracNum], k == 0, angs[indexAngle]);
		}
	}
	void GenFncFitSeries::print(const std::vector<double>& params) const {
		printf("Fourier oeffs\n");
		for (int i = 0; i < indices.size(); i++)
			printf("(%d %d %d) %g  ", indices[i].mr, indices[i].mz, indices[i].mphi, params[i]);
		printf("\n");
		printf("Fractions\n");
		for (int i = 0; i < fracs.size(); i++)
			printf("m: %d, B: %g, b: %g\n",
				fracs[i].mz, params[indices.size() + 2 * i], params[indices.size() + 2 * i + 1]);
	}
	void GenFncFitSeries::write(FILE* ofile, const std::vector<double>& params) const {
		fprintf(ofile, "%zd\n", indices.size());
		for (int i = 0; i < indices.size(); i++)
			fprintf(ofile, "%d %d %d %g\n", indices[i].mr, indices[i].mz, indices[i].mphi, params[i]);
		fprintf(ofile, "%zd\n", fracs.size());
		for (int i = 0; i < fracs.size(); i++)
			fprintf(ofile, "%d %g %g\n",//mz, B, b
				fracs[i].mz, params[indices.size() + 2 * i], params[indices.size() + 2 * i + 1]);
	}
	void GenFncFitSeries::read(FILE* ifile, std::vector<double>& params) {
		indices.clear(); params.clear();
		int mr, mz, mphi; double p;
		int sI; fscanf_s(ifile, "%d", &sI);
		for (int i = 0; i < sI; i++) {
			fscanf_s(ifile, "%d %d %d %lg", &mr, &mz, &mphi, &p);
			indices.push_back(GenFncIndex(mr, mz, mphi));
			params.push_back(p);
		}
		printf("Read %zd terms\n", params.size());
		int nf; double B, b;
		fscanf_s(ifile, "%d", &nf);
		for (int i = 0; i < nf; i++) {
			fscanf_s(ifile, "%d %lg %lg", &mz, &B, &b);
			fracs.push_back(GenFncFitFrac(mz));
			params.push_back(B); params.push_back(b);
		}
		printf("Read %zd fracs\n", fracs.size());
	}

	GenFncIndices GenFncFitSeries::expand(std::vector<double>& params, GenFncFitFracs& fracs,bool fracts)
	{   /// NOTE: here we specialize for the case of axisymmetric systems!
		assert(params.size() == numParams());
		std::map< std::pair<int, int>, double > indPairs;
		int numaddTerm = 0, numaddFrac = 0;
		GenFncIndices newIndices(indices);

		// 1. Store amplitudes & determine the extent of existing grid in (mr,mz)
		int maxmr = 0, maxmz = 0;
		double maxP = 0;
		for (unsigned int i = 0; i < indices.size(); i++) {
			indPairs[std::make_pair(indices[i].mr, indices[i].mz)] = params[i];
			maxmr = std::max<int>(maxmr, math::abs(indices[i].mr));
			maxmz = std::max<int>(maxmz, math::abs(indices[i].mz));
			maxP = fmax(maxP, fabs(params[i]));
		}

		/*if(maxmz==0) {  // dealing with the case Jz==0 -- add only two elements in m_r
			newIndices.push_back(GenFncIndex(maxmr+1, 0, 0));
			newIndices.push_back(GenFncIndex(maxmr+2, 0, 0));
			return newIndices;
		}*/

		// 2. determine the largest amplitude of coefs that are at the boundary of existing values
		double maxval = 0;
		//note terms that are not themselves off end, but have a neighbour that is
		for (int ir = 0; ir <= maxmr; ir++) {
			for (int iz = -maxmz; iz <= maxmz; iz += 2) {
				bool topfrac = haveFrac(iz + 2, fracs);
				bool botfrac = haveFrac(iz - 2, fracs);
				if (indPairs.find(std::make_pair(ir, iz)) != indPairs.end() &&
					((iz <= 0 && !botfrac && indPairs.find(std::make_pair(ir, iz - 2)) == indPairs.end()) ||
						(iz >= 0 && !topfrac && indPairs.find(std::make_pair(ir, iz + 2)) == indPairs.end()) ||
						indPairs.find(std::make_pair(ir + 1, iz)) == indPairs.end()))
					maxval = fmax(fabs(indPairs[std::make_pair(ir, iz)]), maxval);
			}
		}
		// If nothing on the frontier is significant, don't expand
		const double negligible = 1e-4;
		if (maxval < negligible * maxP) return newIndices;

		// 3. Add fractions if appropriate
		std::vector<double> fracPs;
		for (unsigned int j = 0; j < fracs.size(); j++) {//store params associated with existing fracs
			int k = indices.size() + 2 * j;
			fracPs.push_back(params[k]); params[k] = 0;
			fracPs.push_back(params[k + 1]); params[k + 1] = 0;
		}
		//first attempt fit when reach S_nrStart and use terms from S_{nrStart-nrBack}
		if (fracts) {
			int nrStart = 15, nrBack = 5;
			for (int iz = -2; iz <= 8; iz += 2) {// run over ir-series
				if (haveFrac(iz, fracs)) continue;//continue if we already have a frac
				double B, b;
				if (!addFrac(indPairs, maxmr, iz, nrStart, nrBack, B, b)) continue;
				for (int i = 0; i < indices.size(); i++) {// Subtract power captured by frac
					if (indices[i].mz == iz)
						params[i] -= B * pow(b, indices[i].mr - 1);
				}
				numaddFrac++;
				fracs.push_back(GenFncFitFrac(iz));
				fracPs.push_back(B); fracPs.push_back(.5 * log((1 + b) / (1 - b)));//b = tanh(psi) 
			}
		}
		// 4. add more terms adjacent to the existing ones at the boundary, if they are large enough
		double thresh = maxval * 0.2;
		int mstop = maxmr == 0 ? 0 : maxmr + 3;
		for (int iz = -maxmz - 2; iz <= maxmz + 2; iz += 2) {
			if (haveFrac(iz, fracs)) continue;//skip if we already have a frac
			for (int ir = 0; ir <= mstop; ir++) {
				if (indPairs.find(std::make_pair(ir, iz)) != indPairs.end()
					|| (ir == 0 && iz >= 0)) continue;  // already exists or not required
				if (absvalue(indPairs, ir - 3, iz) >= thresh ||
					absvalue(indPairs, ir - 2, iz) >= thresh ||
					absvalue(indPairs, ir - 1, iz) >= thresh ||
					absvalue(indPairs, ir, iz - 2) >= thresh ||
					absvalue(indPairs, ir, iz + 2) >= thresh ||
					absvalue(indPairs, ir + 1, iz - 2) >= thresh ||
					absvalue(indPairs, ir + 1, iz + 2) >= thresh ||
					absvalue(indPairs, ir + 1, iz) >= thresh
					)
				{   // add a term if any of its neighbours are large enough
					newIndices.push_back(GenFncIndex(ir, iz, 0));
					numaddTerm++;
				}
			}
		}

		// 4. finish up
		//printf("Added %d terms, %d fracs\n",numaddTerm,numaddFrac);
		params.resize(newIndices.size() + 2 * fracs.size(), 0);
		for (int i = 0; i < fracPs.size(); i++)//copy prameters of fracs to end of params
			params[newIndices.size() + i] = fracPs[i];
		return newIndices;
	}
}  // namespace actions
