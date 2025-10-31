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
			}
		private:
			const GenFnc* GF;
			const Angles ang;             ///< true angles
		};

	} // internal namespace

	double GenFnc::strength(void) const {
		std::pair<int, int> maxIs(maxIndices());
		double strngth = 0;
		for (int kz = -maxIs.second; kz <= maxIs.second; kz++) {
			int N = 0;
			for (int kr = kz > 0 ? 1 : 0; kr <= maxIs.first; kr++) {
				GenFncIndex ind(kr, kz, 0);
				double val(giveValue(ind));
				if (std::isnan(val)) continue;
				N = kr;
				strngth += val * val;
			}
		}
		return strngth;
	}

	void GenFnc::write(FILE* ofile) const {
		fprintf(ofile, "%zd\n", indices.size()/*, fracs.size()*/);
		for (int i = 0; i < indices.size(); i++)
			fprintf(ofile, "%d %d %d %g %g %g %g\n",
				indices[i].mr, indices[i].mz, indices[i].mphi,
				values[i], derivs[i].Jr, derivs[i].Jz, derivs[i].Jphi);
/*		for (int j = 0; j < fracs.size(); j++) {
			fprintf(ofile, "%d %g %g\n", fracs[j].mz, fracs[j].B, fracs[j].b);
			fprintf(ofile, "%g %g %g %g %g %g\n",
				fracs[j].dBdJ.Jr, fracs[j].dBdJ.Jz, fracs[j].dBdJ.Jphi,
				fracs[j].dpsidJ.Jr, fracs[j].dpsidJ.Jz, fracs[j].dpsidJ.Jphi);
		}*/
	}
	void GenFnc::read(FILE* ifile) {
		int nt, nf;
		fscanf_s(ifile, "%d", &nt/*, &nf*/);
		indices.resize(nt); values.resize(nt); derivs.resize(nt/* + 2 * nf*/);
		for (int i = 0; i < indices.size(); i++)
			fscanf_s(ifile, "%d %d %d %lg %lg %lg %lg\n",
				&indices[i].mr, &indices[i].mz, &indices[i].mphi,
				&values[i], &derivs[i].Jr, &derivs[i].Jz, &derivs[i].Jphi);
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
		return D;
	}

	// dtheta_i/dthetaT_j
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
		return GenFnc(indices, values, derivs/*, fracs*/);
	}

	GenFncFit::GenFncFit(const GenFncIndices& _indices, const int timesr,
		const int timesz,
		const Actions& _acts) :
		indices(_indices), acts(_acts)
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

	std::pair<int, int> GenFncFit::maxIndices() const {
		int mr = 0, mz = 0;
		for (int i = 0; i < indices.size(); i++) {
			mr = std::max<int>(mr, std::abs(indices[i].mr));
			mz = std::max<int>(mz, std::abs(indices[i].mz));
		}
		return std::make_pair(mr, mz);
	}

	ActionAngles GenFncFit::toyActionAngles(unsigned int indexAngle, const double values[]) const
	{
		ActionAngles aa(acts, angs[indexAngle]);
		for (unsigned int indexCoef = 0; indexCoef < indices.size(); indexCoef++) {
			double val = values[indexCoef] * coefs(indexAngle, indexCoef);
			aa.Jr += val * indices[indexCoef].mr;
			aa.Jz += val * indices[indexCoef].mz;
			aa.Jphi += val * indices[indexCoef].mphi;
		}
		return aa;
	}
	Actions GenFncFit::deriv(const unsigned int indexAngle, const unsigned int indexCoef,
		const double values[]) const {
		Actions dJ;
		//Dif wrt S_k
		double val = coefs(indexAngle, indexCoef);  // no range check performed!
		return Actions(
			       val * indices[indexCoef].mr,
			       val * indices[indexCoef].mz,
			       val * indices[indexCoef].mphi);
	}
	void GenFncFit::print(const std::vector<double>& params) const {
		printf("Fourier oeffs\n");
		for (int i = 0; i < indices.size(); i++)
			printf("(%d %d %d) %g  ", indices[i].mr, indices[i].mz, indices[i].mphi, params[i]);
		printf("\n");
	}
	void GenFncFit::write(FILE* ofile, const std::vector<double>& params) const {
		fprintf(ofile, "%zd\n", indices.size());
		for (int i = 0; i < indices.size(); i++)
			fprintf(ofile, "%d %d %d %g\n", indices[i].mr, indices[i].mz, indices[i].mphi, params[i]);
	}
	void GenFncFit::read(FILE* ifile, std::vector<double>& params) {
		indices.clear(); params.clear();
		int mr, mz, mphi; double p;
		int sI; fscanf_s(ifile, "%d", &sI);
		for (int i = 0; i < sI; i++) {
			fscanf_s(ifile, "%d %d %d %lg", &mr, &mz, &mphi, &p);
			indices.push_back(GenFncIndex(mr, mz, mphi));
			params.push_back(p);
		}
		printf("Read %zd terms\n", params.size());
	}

	GenFncIndices GenFncFit::expand(std::vector<double>& params/*, GenFncFitFracs& fracs,bool fracts*/)
	{   /// NOTE: here we specialize for the case of axisymmetric systems!
		assert(params.size() == numParams());
		std::map< std::pair<int, int>, double > indPairs;
		int numaddTerm = 0/*, numaddFrac = 0*/;
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

		// 2. determine the largest amplitude of coefs that are at the boundary of existing values
		double maxval = 0;
		//note terms that are not themselves off end, but have a neighbour that is
		for (int ir = 0; ir <= maxmr; ir++) {
			for (int iz = -maxmz; iz <= maxmz; iz += 2) {
				bool topfrac = false; //haveFrac(iz + 2, fracs);
				bool botfrac = false; //haveFrac(iz - 2, fracs);
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

		// 3. add more terms adjacent to the existing ones at the boundary, if they are large enough
		double thresh = maxval * 0.2;
		int mstop = maxmr == 0 ? 0 : maxmr + 3;
		for (int iz = -maxmz - 2; iz <= maxmz + 2; iz += 2) {
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
		params.resize(newIndices.size(), 0);
		return newIndices;
	}
}  // namespace actions
