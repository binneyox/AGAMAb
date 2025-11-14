#include "actions_newtorus.h"
#include "orbit.h"
#ifdef _OPENMP
#include <omp.h>
#endif

/// accuracy parameter determining the spacing of the interpolation grid along the energy axis
static const double ACCURACY_INTERP2 = 1e-5;

namespace actions {

	namespace {//internal

		double NANfrac, NANbar;

		//Two functions used in containsPoint
		double angleDiff(Angles A1, Angles A2) {

			double dr = A1.thetar - A2.thetar;
			while (dr > M_PI) dr -= 2 * M_PI;
			while (dr < -M_PI) dr += 2 * M_PI;
			double dz = A1.thetaz - A2.thetaz;
			while (dz > M_PI) dz -= 2 * M_PI;
			while (dz < -M_PI) dz += 2 * M_PI;
			double dp = A1.thetaphi - A2.thetaphi;
			while (dp > M_PI) dp -= 2 * M_PI;
			while (dp < -M_PI) dp += 2 * M_PI;
			return sqrt(dr * dr + dz * dz + dp * dp);
		}

		bool is_new(Angles A1, std::vector<Angles> As) {
			bool ok = true; const double tiny = 1e-5;
			for (int i = 0; i < As.size(); i++) {
				double diff = angleDiff(A1, As[i]);
				ok = ok && diff > tiny;
			}
			return ok;
		}
		class locFinder : public math::IFunctionNdimDeriv {
		private:
			const Torus T;
			const coord::PosCyl P0;
		public:
			locFinder(const Torus _T, const coord::PosCyl& _P0) : T(_T), P0(_P0) {}
			virtual void evalDeriv(const double params[], double* values, double* derivs = NULL) const {
				Angles A1(params[0], params[1], 0);
				coord::PosCyl P1;  DerivAng<coord::Cyl> dA;
				if (derivs) {
					P1 = (T.new_PosDerivs(A1, dA));
					derivs[0] = dA.dbythetar.R; derivs[1] = dA.dbythetaz.R;
					derivs[2] = dA.dbythetar.z; derivs[3] = dA.dbythetaz.z;
				}
				if (values) {
					if(!derivs) P1 = T.from_toy(A1);
					values[0] = P1.R - P0.R; values[1] = P1.z - P0.z;
				}
			}
			virtual unsigned int numVars() const {
				return 2;
			}
			virtual unsigned int numValues() const {
				return 2;
			}
		};

		double insertLine(int& ntop, const int tmax, double s, GenFncIndex I,
			std::vector<double>& Hmods, GenFncIndices& Hindices) {
			if (ntop == 0) {
				Hmods.push_back(s);
				Hindices.push_back(I);
				ntop++; return Hmods[0];
			}
			int l = 0;
			while (l < ntop && s <= Hmods[l]) l++;
			if (l == ntop) {//line isn't stronger than any previous line
				if (ntop >= tmax) {
					return Hmods[ntop - 1];//no room for this line
				}
				else {//add line to end of list
					Hmods.push_back(s);
					Hindices.push_back(I);
					ntop++;
					return Hmods.back();
				}
			}
			else {//we should insert current line
				if (ntop < tmax) {//move existing terms down
					Hmods.push_back(Hmods[ntop - 1]);
					Hindices.push_back(Hindices[ntop - 1]);
				}
				for (int m = ntop - 1; m > l; m--) {
					Hmods[m] = Hmods[m - 1];
					Hindices[m] = Hindices[m - 1];
				}//we've created a space at l, so fill it
				Hmods[l] = s;
				Hindices[l] = I;
				if (ntop < tmax) ntop++;//we've added rather than replaced a line
				return Hmods.back();//return weakest retained line
			}
		}
		//scaling of E for finding planar orbit energy
		inline double scaleE(const double E, const double invPhi0, /*output*/ double* dEdscaledE=NULL)
		{
			double expX = invPhi0 - 1/E;
			if(dEdscaledE)
				*dEdscaledE = E * E * expX;
			return log(expX);
		}
		//return E and optionally dE/d(scaledE) as a function of scaledE and invPhi0
		inline double unscaleE(const double scaledE, const double invPhi0, /*output*/ double* dEdscaledE=NULL)
		{
			double expX = exp(scaledE);
			double E = 1 / (invPhi0 - expX);
			if(dEdscaledE)
				*dEdscaledE = E * E * expX;
			return E;
		}

		// H & its deriv of H wrt R,z,p_R,p_z,p_phi
		double H_dHdX(const potential::BasePotential& pot, const coord::PosMomCyl Rzphi,
			coord::PosMomCyl& dHdX) {
			double Phi; coord::GradCyl grad;
			pot.eval(Rzphi, &Phi, &grad);
			dHdX.R = grad.dR;
			if(Rzphi.pphi!=0)dHdX.R+= - pow_2(Rzphi.pphi / Rzphi.R) / Rzphi.R;
			dHdX.z = grad.dz; dHdX.phi = grad.dphi;
			dHdX.pR = Rzphi.pR; dHdX.pz = Rzphi.pz; 
			dHdX.pphi =(Rzphi.pphi!=0)? Rzphi.pphi / pow_2(Rzphi.R):0;
			double H = .5 * (pow_2(Rzphi.pR) + pow_2(Rzphi.pz)) + Phi;
			if (Rzphi.pphi != 0)H += .5 * pow_2(Rzphi.pphi / Rzphi.R);
			return H;
		}
		//Class used in choice of isochrone
		class Iso {
		private:
			double L, Jr;
		public:
			Iso(double _L, double _Jr) : L(_L), Jr(_Jr) {}
			double b2cE(double Js) const {//general E
				return .5 * pow_2(Js * Js) / pow_2(Jr + .5 * (L + sqrt(L * L + 4 * Js * Js)));
			}
			double b2cEc(double Js) const {//circular E
				return .5 * pow_2(Js * Js) / pow_2(.5 * (L + sqrt(L * L + 4 * Js * Js)));
			}
			double g(double Js) const {//Rsh=b*g(Js)
				double g2 = (pow((L + sqrt(L * L + 4 * Js * Js)) / (2 * Js), 4) - 1);
				if (g2 < 0) printf("g2<0: %f\n", g2);
				return sqrt(g2);
			}
			double cob(double Js) const {//c/b
				return .5 * pow_2(Js) / b2cE(Js) - 1;;
			}
			double ecc(double Js) const {
				double boc = 1 / cob(Js);
				double e2 = 1 - pow_2(L / Js) * boc * (1 + boc);
				return e2 < 1 ? sqrt(e2) : 0;
			}
			double f(double b, double Js, double& e) const {//ratio of forces aopo/peri
				double cb = cob(Js), boc = 1 / cb;
				e = ecc(Js);
				double up = 1 + e, um = 1 - e;
				double rp = b * cb * sqrt(up * (up + 2 * boc)), ap = sqrt(b * b + rp * rp);
				double rm = b * cb * sqrt(um * (um + 2 * boc)), am = sqrt(b * b + rm * rm);
				return pow_2((b + am) / (b + ap)) * am / ap * rp / rm;
			}
		};

		/*
		 Class to find Js (and thus b) of isochrone for which the isochrone's
		 force ratio matches that in real Phi
		*/
		class JsFinder : public math::IFunction {
			const potential::BasePotential& pot;
			const Iso Is;
			const double Rsh;
		public:
			JsFinder(const potential::BasePotential& _pot, const Iso& _Is, const double _Rsh) :
				pot(_pot), Is(_Is), Rsh(_Rsh) {
			}
			virtual void evalDeriv(double Js, double* val, double* deriv = 0, double* deriv2 = 0) const {
				double b = Rsh / Is.g(Js);
				double e, f_apo_peri = Is.f(b, Js, e);
				double c = Is.cob(Js) * b;
				double F[2];
				for (int k = -1; k < 2; k += 2) {
					double u = 1 + k * e;
					double Phi, r = c * sqrt(u * (u + 2 * b / c));
					coord::PosCyl Rz(r, 0, 0); coord::GradCyl grad;
					pot.eval(Rz, &Phi, &grad);
					F[(k + 1) / 2] = grad.dR;
				}
				if (val) *val = f_apo_peri - F[1] / F[0];
			}
			virtual unsigned int numDerivs(void) const {
				return 0;
			}

		};

		/// create the array of indices of the generating function with all terms up to the given maximum order
		static GenFncIndices makeGridIndices(int irmax, int izmax)
		{   /// NOTE: here we specialize for the case of axisymmetric systems!
			GenFncIndices indices;
			for (int ir = 0; ir <= irmax; ir++)
				for (int iz = -izmax; iz <= (ir == 0 ? -2 : izmax); iz += 2)
					indices.push_back(GenFncIndex(ir, iz, 0));
			return indices;
		}

		/// return the absolute value of an element in a map, or zero if it doesn't exist
		static inline double absvalue(const std::map< std::pair<int, int>, double >& indPairs, int ir, int iz)
		{
			if (indPairs.find(std::make_pair(ir, iz)) != indPairs.end())
				return fabs(indPairs.find(std::make_pair(ir, iz))->second);
			else
				return 0;
		}
		/** Compute the derivative of Hamiltonian by toy actions:
				dH/dJ = dH/d{x,v} d{x,v}/dJ, where the lhs is a covector of length 3,
				the first term on rhs is a covector of length 6 (the gradient dPhi/dx and the velocity),
				and the second term is a 6x3 matrix of partial derivs provided by the toy map.
		*/
		static inline Actions dHbydJ(const coord::PosMomCyl& dHby,
			const DerivAct<coord::Cyl>& dXdJ)
		{
			return Actions(
				dXdJ.dbyJr.R * dHby.R + dXdJ.dbyJr.z * dHby.z + dXdJ.dbyJr.phi * dHby.phi +
				dXdJ.dbyJr.pR * dHby.pR + dXdJ.dbyJr.pz * dHby.pz + dXdJ.dbyJr.pphi * dHby.pphi,
				dXdJ.dbyJz.R * dHby.R + dXdJ.dbyJz.z * dHby.z + dXdJ.dbyJz.phi * dHby.phi +
				dXdJ.dbyJz.pR * dHby.pR + dXdJ.dbyJz.pz * dHby.pz + dXdJ.dbyJz.pphi * dHby.pphi,
				dXdJ.dbyJphi.R * dHby.R + dXdJ.dbyJphi.z * dHby.z + dXdJ.dbyJphi.phi * dHby.phi +
				dXdJ.dbyJphi.pR * dHby.pR + dXdJ.dbyJphi.pz * dHby.pz + dXdJ.dbyJphi.pphi * dHby.pphi);
		}
		//class to help with zSoS computation
		class CrossingFinder : public math::IFunction {
		private:
			const Torus* T;
			const double thetaT_r, z0;
		public:
			CrossingFinder(const Torus* _T, const double& _thetaT, const double _z0 = 0) :
				T(_T), thetaT_r(_thetaT), z0(_z0) {
			}
			virtual void evalDeriv(const double thetaT_z, double* value, double* deriv = NULL, double* deriv2 = NULL) const {
				Angles thetaT(thetaT_r, thetaT_z, 0);
				coord::PosMomCyl Rz(T->from_toy(thetaT));
				*value = Rz.z - z0;
			}
			virtual unsigned int numDerivs(void) const {
				return 0;
			}
		};
		//class to help with RSoS computation
		class r_crossingFinder : public math::IFunction {
		private:
			const Torus* T;
			const double rbar;
			const double thetaT_z;
			const coord::ProlSph coordsys;
		public:
			r_crossingFinder(const Torus* _T, const double _rbar, const double& _thetaT, const double _delta) :
				T(_T), rbar(_rbar), thetaT_z(_thetaT), coordsys(_delta) {
			}
			virtual void evalDeriv(const double thetaT_r, double* value, double* deriv = NULL, double* deriv2 = NULL) const {
				Angles thetaT(thetaT_r, thetaT_z, 0);
				coord::PosMomCyl Rv = T->from_toy(thetaT);
				coord::PosVelProlSph pps = coord::toPosVel<coord::Cyl, coord::ProlSph>(coord::toPosVelCyl(Rv), coordsys);
				double r = sqrt(pps.lambda - coordsys.Delta2);
				*value = r - rbar;
			}
			virtual unsigned int numDerivs(void) const {
				return 0;
			}
		};
		/*
		 * Function to comput dx,dp/d true angles. It
		 * takes in (dx,dp)/dthetaT and dtheta/dthetaT, inverts the
		 * latter and multiplies on the former
		 */
		void new_assemble(DerivAng<coord::Cyl>& dA, math::Matrix<double>& M) {
			math::Matrix<double> dxdt(3, 3);
			dxdt(0, 0) = dA.dbythetar.R; dxdt(0, 1) = dA.dbythetaz.R;
			dxdt(0, 2) = dA.dbythetaphi.R;
			dxdt(1, 0) = dA.dbythetar.z; dxdt(1, 1) = dA.dbythetaz.z;
			dxdt(1, 2) = dA.dbythetaphi.z;
			dxdt(2, 0) = dA.dbythetar.phi; dxdt(2, 1) = dA.dbythetaz.phi;
			dxdt(2, 2) = dA.dbythetaphi.phi;
			math::Matrix<double> dpdt(3, 3);
			dpdt(0, 0) = dA.dbythetar.pR; dpdt(0, 1) = dA.dbythetaz.pR;
			dpdt(0, 2) = dA.dbythetaphi.pR;
			dpdt(1, 0) = dA.dbythetar.pz; dpdt(1, 1) = dA.dbythetaz.pz;
			dpdt(1, 2) = dA.dbythetaphi.pz;
			dpdt(2, 0) = dA.dbythetar.pphi; dpdt(2, 1) = dA.dbythetaz.pphi;
			dpdt(2, 2) = dA.dbythetaphi.pphi;
			math::LUDecomp L2(M);
			math::Matrix<double> M3 = (L2.inverse(3));
			math::Matrix<double> M2(3, 3);
			math::blas_dgemm(math::CblasNoTrans, math::CblasNoTrans, 1.0,
				dxdt, M3, 0.0, M2);
			dA.dbythetar.R = M2(0, 0);
			dA.dbythetaz.R = M2(0, 1);
			dA.dbythetaphi.R = M2(0, 2);
			dA.dbythetar.z = M2(1, 0);
			dA.dbythetaz.z = M2(1, 1);
			dA.dbythetaphi.z = M2(1, 2);
			dA.dbythetar.phi = M2(2, 0);
			dA.dbythetaz.phi = M2(2, 1);
			dA.dbythetaphi.phi = M2(2, 2);
			math::blas_dgemm(math::CblasNoTrans, math::CblasNoTrans, 1.0,
				dpdt, M3, 0.0, M2);
			dA.dbythetar.pR = M2(0, 0);
			dA.dbythetaz.pR = M2(0, 1);
			dA.dbythetaphi.pR = M2(0, 2);
			dA.dbythetar.pz = M2(1, 0);
			dA.dbythetaz.pz = M2(1, 1);
			dA.dbythetaphi.pz = M2(1, 2);
			dA.dbythetar.pphi = M2(2, 0);
			dA.dbythetaz.pphi = M2(2, 1);
			dA.dbythetaphi.pphi = M2(2, 2);
		}

		double avgH;
		void setH(double Hbar) { avgH = Hbar; }

		/* A helper class used by TorusGenerator. It evaluates H and its
		 * derivatives on a toy-angle grid over a proposed torus using a
		 * proposed GF. Once the GF has been
		 * optimised, it runs over the grid a final time to determine d S_k/d J
		 * so we can recover the true angles
		 */
		class torusFitter : public math::IFunctionNdimDeriv {
		private:
			const Actions& J;
			const potential::BasePotential& pot;
			const double freqScale;
			GenFncFit& GFFS;
			const PtrToyMap ptrTM;
			int N1;
		public:
			torusFitter(const Actions& _J,
				const potential::BasePotential& _pot,
				const double _freqScale,
				const PtrToyMap _ptrTM, GenFncFit& _GFFS) :
			    J(_J), pot(_pot), freqScale(1000 * _freqScale), GFFS(_GFFS), ptrTM(_ptrTM)
			{
				N1 = GFFS.numParams();
			}
			virtual unsigned int numVars() const {
				return N1;
			}
			virtual unsigned int numValues() const {
				return GFFS.numPoints();
			}
			double new_computeHamiltonianAtPoint(const double params[],
				const unsigned int indPoint, Actions& JT, Actions* dHdJ = NULL,
				double* dHdParams = NULL) const
			{
				// Generating function computes the toy actions from the real actions
				// at the given point in the grid of toy angles grid
				ActionAngles toyAA = GFFS.toyActionAngles(indPoint, params);

				// do not allow to stray into forbidden region of negative actions
				if (toyAA.Jr < 0 || toyAA.Jz < 0) {
					JT.Jr = toyAA.Jr; JT.Jz = toyAA.Jz;
					double H = toyAA.Jr < 0 ? -freqScale * toyAA.Jr : 0;
					H -= toyAA.Jz < 0 ? freqScale * toyAA.Jz : 0;
					return H;
				}

				DerivAct<coord::Cyl> dXdJ;
				coord::PosMomCyl Rzphi(ptrTM->from_aaT(toyAA, dXdJ));
		// obtain the value of the real Hamiltonian at the given point and its
		// derivatives w.r.t. coordinates/momenta
				coord::PosMomCyl dHdX;
				double H = H_dHdX(pot, Rzphi, dHdX);
				if (dHdParams) {
					// derivative of Hamiltonian by toy actions
					Actions dHby = dHbydJ(dHdX, dXdJ);
					if (dHdJ) *dHdJ = dHby;
					for (unsigned int p = 0; p < GFFS.numParams(); p++) {
						// derivs of toy actions by gen.fnc.params
						Actions dbyP = GFFS.deriv(indPoint, p, &params[0]);
						dHdParams[p] = dHby.Jr * dbyP.Jr + dHby.Jz * dbyP.Jz
							+ dHby.Jphi * dbyP.Jphi;
					}
				}
				return H;
			}
			double computeHamiltonianDisp(const std::vector<double>& params, double& Hbar)
			{
				const unsigned int max_threads = 16;
				double Hsm[max_threads] = { 0 }, Hsq[max_threads] = { 0 };
				int N[max_threads] = { 0 };
				int Nnan[max_threads] = { 0 };
			#ifdef _OPENMP
				#pragma omp parallel for schedule(dynamic)
				for (int indPoint = 0; indPoint < GFFS.numPoints(); indPoint++) {
					int nth = omp_get_thread_num();
					Actions JT(1, 1, 1);
					double H = new_computeHamiltonianAtPoint(&params[0], indPoint, JT);
					if (JT.Jr < 0 || JT.Jz < 0) Nnan[nth]++;
					{
						Hsm[nth] += H; Hsq[nth] += H * H; N[nth]++;
					}
				}
				for (int i = 1; i < max_threads; i++) {//concatenate sums
					Nnan[0] += Nnan[i]; Hsm[0] += Hsm[i]; Hsq[0] += Hsq[i]; N[0] += N[i];
				}
			#else
				for (int indPoint = 0; indPoint < GFFS.numPoints(); indPoint++) {
					Actions JT(1, 1, 1);
					double H = new_computeHamiltonianAtPoint(&params[0], indPoint, JT);
					if (JT.Jr < 0 || JT.Jz < 0) Nnan[0]++;
					Hsm[0] += H; Hsq[0] += H * H; N[0]++;
				}
			#endif
				Hbar = Hsm[0] / N[0];
				Hsq[0] = Hsq[0] / N[0] - Hbar * Hbar;
				NANfrac = (double)Nnan[0] / (double)GFFS.numPoints();
				NANbar += NANfrac;
				return Hsq[0] > 0 ? sqrt(Hsq[0]) : 0;
			}
			void evalDeriv(const double params[],
				double* deltaHvalues, double* dHdParams) const
			{
				const unsigned int numPoints = GFFS.numPoints();
				const unsigned int numParams = N1;
				const unsigned int max_threads = 16;

				// we need to store the values of Hamiltonian at grid points even if this is not requested,
				// because they are used to correct the entries of the Jacobian matrix
				// to account for the fact that the mean <H> also depends on the parameters
				std::vector<double> Hvalues(numPoints);
				double Havg[max_threads] = { 0 };  // accumulator for the average Hamiltonian
				int nNAN[max_threads] = { 0 };
				// loop over grid of toy angles
			#ifdef _OPENMP
				#pragma omp parallel for schedule(dynamic)
				for (int indPoint = 0; indPoint < numPoints; indPoint++) {
					int nth = omp_get_thread_num();
					Actions JT(1, 1, 1);
					double H = new_computeHamiltonianAtPoint(params, indPoint, JT, NULL,
						dHdParams ? dHdParams + indPoint * numParams : NULL);
					if (JT.Jr < 0 || JT.Jz < 0)nNAN[nth]++;
					else {
						// accumulate the average value and store the output
						Havg[nth] += H;
						Hvalues[indPoint] = H;
					}
				}
				for (int i = 1; i < max_threads; i++) {
					nNAN[0] += nNAN[i]; Havg[0] += Havg[i];
				}
			#else
				for (int indPoint = 0; indPoint < numPoints; indPoint++) {
					Actions JT(1, 1, 1);
					double H = new_computeHamiltonianAtPoint(params, indPoint, JT, NULL,
						dHdParams ? dHdParams + indPoint * numParams : NULL);
					if (JT.Jr < 0 || JT.Jz < 0)nNAN[0]++;
					else {
						// accumulate the average value and store the output
						Havg[0] += H;
						Hvalues[indPoint] = H;
					}
				}
			#endif
				NANfrac = (double)nNAN[0] / numPoints;
				NANbar += NANfrac;
				// convert from  H_k  to  deltaH_k = H_k - <H>
				Havg[0] /= numPoints; setH(Havg[0]);
				if (deltaHvalues) {
					//double disp = 0;
					for (unsigned int indPoint = 0; indPoint < numPoints; indPoint++) {
						deltaHvalues[indPoint] = (Hvalues[indPoint] - Havg[0]);
						//disp += pow_2(deltaHvalues[indPoint]);
					}
				}

				// convert derivatives:  d(deltaH_k) / dP_p = dH_k / dP_p - d<H> / dP_p
				if (dHdParams) {
					std::vector<double> dHavgdP(numPoints);
					for (unsigned int p = 0; p < numParams; p++) dHavgdP[p] = 0;
					for (unsigned int pp = 0; pp < numPoints * numParams; pp++) {
						if(!std::isnan(dHdParams[pp])) dHavgdP[pp % numParams] += dHdParams[pp] / numPoints;
						else printf("nan@ %d %d\n", pp % numParams, pp / numParams);
					}
					for (unsigned int pp = 0; pp < numPoints * numParams; pp++) {
						unsigned int indParam = pp % numParams;
						dHdParams[pp] = dHdParams[pp] - dHavgdP[indParam];
					}
				}
			}
			void testit(std::vector<double> params, int N) {
				int nPts = GFFS.numPoints(), nPars = GFFS.numParams();
				double* HmHbar = new double[nPts];
				double* HmX = new double[nPts];
				double* dHdP = new double[nPts * nPars];
				double dA = 1e-5;
				//explore changing Nth parameter 
				params[N] += .5 * dA;
				evalDeriv(&params[0], HmHbar, NULL);
				params[N] -= dA;
				evalDeriv(&params[0], HmX, NULL);
				for (int i = 0;i < nPts;i++)
					HmX[i] = (HmHbar[i] - HmX[i]) / dA;
				params[N] += .5 * dA;
				evalDeriv(&params[0], HmHbar, dHdP);
				for (int i = 0; i < nPts; i += 10)
					printf("(%g %g) ", HmX[i], dHdP[i * nPars + N]);
				printf("\n");
				delete[] HmHbar; delete[] HmX; delete[] dHdP;
			}
			/** Compute the frequencies and the derivatives of generating function by real actions,
				used in angle mapping.
				The three arrays of derivatives dS_i/dJ_{r,z,phi}, i=0..numParamsGenFnc-1,
				together with three frequencies Omega_{r,z,phi}, are the solutions of
				an overdetermined system of linear equations:
				\f$  M_{k,i} X_{i} = RHS_{k}, k=0..numPoints-1  \f$,
				where numPoints is the number of individual triplets of toy angles,
				\f$  X_i  \f$ is the solution vector {Omega, dS_i/dJ} for each direction (r,z,theta),
				\f$  RHS_k = dH/dJ(\theta_k)  \f$, again for three directions independently, and
				\f$  M_{k,i}  \f$ is the matrix of coefficients shared between all three equation systems:
				\f$  M_{k,0} = 1, M_{k,i+1} = -dH/dS_i(\theta_k)  \f$.
				The matrix M and three RHS vectors are filled using the same approach as during
				the Levenberg-Marquardt minimization, from the provided parameters of toy map and
				generating function; then the three linear systems are solved using
				the singular-value decomposition of the shared coefficient matrix,
				and the output frequencies and gen.fnc.derivatives are returned in corresponding arguments.
				The return value of this function is the same as `computeHamiltonianDisp()`.
			*/
			double fitAngleMap(const double params[],
					   double& Hbar, Frequencies& freqs, GenFncDerivs& dPdJ, 
					   bool& negJr, bool& negJz) const {
				unsigned int numPoints = GFFS.numPoints();
				unsigned int numParams = GFFS.numParams();
				// the matrix of coefficients shared between three linear systems
				math::Matrix<double> coefsdHdS(numPoints, numParams + 1);
				// tmp storage for dH/dS
				std::vector<double> dHdParams(numParams + 1);
				// derivs of Hamiltonian by toy actions (RHS vectors)
				std::vector<double> dHdJr(numPoints), dHdJz(numPoints), dHdJphi(numPoints);
				// accumulator for computing dispersion in H
				math::Averager Havg;
				negJr = false; negJz=false;
				// loop over grid of toy angles
				for (unsigned int indPoint = 0; indPoint < numPoints; indPoint++) {
					Actions JT(1, 1, 1), dHby;  // derivative of Hamiltonian by toy actions
					double H = new_computeHamiltonianAtPoint(&params[0], indPoint, JT,
						&dHby, &dHdParams.front());
					if (JT.Jr >= 0 && JT.Jz >= 0) Havg.add(H);
					else{
						if(!negJr && JT.Jr>0) negJr = true;
						if(!negJz && JT.Jz>0) negJz = true;
					}
					// fill the elements of each of three rhs vectors
					dHdJr[indPoint] = dHby.Jr;
					dHdJz[indPoint] = dHby.Jz;
					dHdJphi[indPoint] = dHby.Jphi;
					// fill the matrix row
					coefsdHdS(indPoint, 0) = 1;  // matrix coef for omega
					for (unsigned int p = 0; p < numParams; p++) // matrix coefs for dSk/dJ
						coefsdHdS(indPoint, p + 1) = -dHdParams[p];
				}
				Hbar = Havg.mean();
				// solve the overdetermined linear system in the least-square sense:
				// step 1: prepare the SVD of coefs matrix
				math::SVDecomp SVD(coefsdHdS);

				// step 2: solve three linear systems with the same matrix but different rhs
				std::vector<double> dSdJr(SVD.solve(dHdJr)), dSdJz(SVD.solve(dHdJz)),
				dSdJphi(SVD.solve(dHdJphi));
	// A diagnostic section
			/*	double resr=0,resz=0,resp=0;
				for(int i=0;i<numPoints;i++){
					double dfr=dHdJr[i];
					double dfz=dHdJz[i];
					double dfp=dHdJphi[i];
					for(int p=0;p<numParams+1;p++){
						dfr-=coefsdHdS(i,p)*dSdJr[p];
						dfz-=coefsdHdS(i,p)*dSdJz[p];
						dfp-=coefsdHdS(i,p)*dSdJphi[p];
					}
					resr+=(pow_2(dfr))/(numPoints);
					resz+=(pow_2(dfz))/(numPoints);
					resp+=(pow_2(dfp))/(numPoints);
				}
				printf("Omega: (%f,%f,%f) with residuals: (%f,%f,%f)\n",
				       dSdJr[0],dSdJz[0],dSdJphi[0],sqrt(resr),sqrt(resz),sqrt(resp));
				printf("number of points:%d with %d parameters\n",numPoints,numParams);*/
				// store output
				freqs.Omegar = dSdJr[0];
				freqs.Omegaz = dSdJz[0];
				freqs.Omegaphi = J.Jphi >= 0 ? dSdJphi[0] : -dSdJphi[0];
				dPdJ.resize(numParams);
				for (unsigned int p = 0; p < numParams; p++) {
					dPdJ[p].Jr = dSdJr[p + 1];
					dPdJ[p].Jz = dSdJz[p + 1];
					dPdJ[p].Jphi = dSdJphi[p + 1];
				}
				return sqrt(Havg.disp());
			}

		};
		//fit v(theta) and rn(r) map for isochrone
		class fitmap : public math::IFunctionNdimDeriv {
		private:
            int N, Nr;
            int nr, nz;
			int np;
			Actions J;
			Isochrone is;
			const potential::BasePotential& pot;
			math::ScalingInfTh sc;
			double delta;
		public:
			fitmap(const int _N, const int _Nr, Actions _J,
				Isochrone _is, const potential::BasePotential& _pot, double _R0, double _delta)
				: N(_N), Nr(_Nr), nr(100), nz(100), np(100 * 100),
				J(_J), is(_is), pot(_pot), sc(_R0), delta(_delta) {
			}
			virtual unsigned int numVars() const {
				return Nr + N;
			}
			virtual unsigned int numValues() const {
				return np;
			}
			double H(const double params[], int i, double* dHdParams = NULL) const {
				int ir = i / nr;int iz = i % nr;
				Angles aT(ir * M_PI / (double)nr, iz * M_PI / (double)nz, 0);
				ActionAngles aaT(J, aT);
				std::vector<double> p(N),pr(Nr);
				for(int i=0;i<Nr;i++){
					pr[i]=params[i];
				}
				for(int i=0;i<N;i++){
					p[i]=params[i+Nr];
				}
				PTIso PT(delta,sc,p,pr);
				ToyMapIso TMis(is,PT);
				std::vector<coord::PosMomCyl> dRzdP(N+Nr);
				coord::PosMomCyl Rz=TMis.from_aaT(aaT,NULL,NULL,dHdParams?&dRzdP[0]:NULL,dHdParams?&dRzdP[Nr]:NULL);
				coord::PosMomCyl dHdX;
				double E=H_dHdX(pot,Rz,dHdX);
				if(dHdParams){
					for(int i=0;i<N+Nr;i++){
						dHdParams[i]=dRzdP[i].R*dHdX.R+dRzdP[i].z*dHdX.z+dRzdP[i].pR*dHdX.pR+dRzdP[i].pz*dHdX.pz;
					}
				}
				return E;
			}
			virtual void evalDeriv(const double params[],
				double* dH, double* dfdParams)const
			{
				std::vector<double>Hval(np);
				int Nt = N + Nr;
				std::vector<double> dfdparams0(Nt * np);
				double Hav = 0.0;
				std::vector<double> dfdav(Nt, 0.0);
				std::vector<double> dfdp(Nt);
				for (int i = 0;i < np;i++) {
					//for (int k = 0;k < np;k++) {
						//double E = H(&params[0], actions::Angles(thetar[i], thetaz[k], 0.0), &dfdp[0]);
					double E = H(&params[0], i, &dfdp[0]);
					Hav += E;
					Hval[i] = E;
					if (dfdParams) {
						for (int j = 0;j < Nt;j++) {
							dfdparams0[i * Nt + j] = dfdp[j];
							dfdav[j] += dfdp[j];
						}
					}
				}
				Hav /= np;
				for (int i = 0;i < Nt;i++) {
					dfdav[i] /= np;
				}
				//double disp = 0.0;
				if (dH) {
					for (int i = 0;i < np;i++) {
						dH[i] = (Hval[i] - Hav) / Hav;
						//disp += pow_2(dH[i]) / np;
					}
				}
				if (dfdParams) {
					for (int i = 0;i < np;i++) {
						for (int j = 0;j < Nt;j++) {
							dfdParams[i * Nt + j] = dfdparams0[i * Nt + j] / Hav - Hval[i] / pow_2(Hav) * dfdav[j];
						}
					}
				}
			}
		};
		//fits v(z1) and R(R1) for Harmonic Oscillator Toy map
		class fitmapHarm : public math::IFunctionNdimDeriv {
			private:
				int N, Nr, N1;
				int nr, nz;
				int np;
				Actions J;
				HarmonicOscilattor os;
				const potential::BasePotential& pot;
				math::ScalingInfTh sc;
				math::ScalingInfTh scz;
				double delta;
			public:
				fitmapHarm(const int _N, const int _Nr, Actions _J,
					   HarmonicOscilattor _os,
					   const potential::BasePotential& _pot, double _R0, double _z0, double _delta)
						: N(_N),Nr(_Nr), N1(_Nr + _N), nr(100), nz(100), np(100 * 100),
				J(_J), os(_os), pot(_pot), sc(_R0), scz(_z0), delta(_delta)
				{
				}
				virtual unsigned int numVars() const {
					return N1;
				}
				virtual unsigned int numValues() const {
					return np;
				}
				double H(const double params[], int i, double* dHdParams = NULL) const {
					int ir = i / nr;int iz = i % nr;
					Angles aT(ir * M_PI / (double)nr, iz * M_PI / (double)nz, 0);
					ActionAngles aaT(J, aT);std::vector<double> p(N),pr(Nr);
					for(int i=0;i<Nr;i++){
						pr[i]=params[i];
					}
					for(int i=0;i<N;i++){
						p[i]=params[i+Nr];
					}
					PTHarm PT(delta,sc,scz,p,pr);
					ToyMapHarm TMH(os,PT);
					std::vector<coord::PosMomCyl> dRzdP(N+Nr);
					coord::PosMomCyl Rz=TMH.from_aaT(aaT,NULL,NULL,dHdParams?&dRzdP[0]:NULL,dHdParams?&dRzdP[Nr]:NULL);
					coord::PosMomCyl dHdX;
					double E=H_dHdX(pot,Rz,dHdX);
					if(dHdParams){
						for(int i=0;i<N+Nr;i++){
							dHdParams[i]=dRzdP[i].R*dHdX.R+dRzdP[i].z*dHdX.z+dRzdP[i].pR*dHdX.pR+dRzdP[i].pz*dHdX.pz;
						}
					}
					return E;
				}
				virtual void evalDeriv(const double params[],
					double* dH, double* dfdParams)const
				{
					std::vector<double>Hval(np);
					std::vector<double> dfdparams0(N1 * np);
					double Hav = 0.0;
					std::vector<double> dfdav(N1, 0.0);
					std::vector<double> dfdp(N1);
					for (int i = 0;i < np;i++) {
						double E = H(&params[0], i, &dfdp[0]);
						Hav += E;
						Hval[i] = E;
						if (dfdParams) {
							for (int j = 0;j < N1;j++) {
								dfdparams0[i * N1 + j] = dfdp[j];
								dfdav[j] += dfdp[j];
							}
						}
					}
					Hav /= np;
					for (int i = 0;i < N1;i++) {
						dfdav[i] /= np;
					}
					if (dH) {
						for (int i = 0;i < np;i++) {
							dH[i] = Hval[i] - Hav;
						}
					}
					if (dfdParams) {
						for (int i = 0;i < np;i++) {
							for (int j = 0;j < N1;j++) {
								dfdParams[i * N1 + j] = (dfdparams0[i * N1 + j] - dfdav[j]);
							}
						}
					}
				}
		};

		//Finds where x lies in xs[]. returns fractional distance from xs[top]
		double bot_top(const double x, const std::vector<double>& xs,
			int& top, int& bot) {
			if (xs[0] > xs[xs.size() - 1]) {
				top = 0; bot = xs.size() - 1;
			}
			else {
				top = xs.size() - 1; bot = 0;
			}//now top should point to largest x
			double f;
			if ((xs[top] - x) * (x - xs[bot]) < 0.) {//x lies out of grid
				if (x < xs[bot]) {
					top = bot + 1; f = 1;//T = Tgrid[bot];
				}
				else {
					bot = top - 1; f = 0;// T = Tgrid[top];
				}
			}
			else {
				while (abs(top - bot) > 1) {
					int n = (top + bot) / 2;
					if ((xs[top] - x) * (x - xs[n]) >= 0) bot = n;
					else top = n;
				}
				f = (xs[top] - x) / (xs[top] - xs[bot]);//distance from top
			}
			return f;
		}

	}//internal

	EXP Torus InterpTorus(const double x, const Torus& T0, const Torus& T1){
		if (x == 1) return T0;
		if (x == 0) return T1;
		const double xp = 1 - x;
		Actions J = T0.J * x + T1.J * xp;
		Frequencies freqs = T0.freqs * x + T1.freqs * xp;
		GenFnc newGF(interpGenFnc(x, T0.GF, T1.GF));
		PtrToyMap newptrTM(interpPtrToyMap(x, T0.ptrTM, T1.ptrTM));
		double newE = T0.E * x + T1.E * xp;
		bool negJr = T0.negJr || T1.negJr;
		bool negJz = T0.negJz || T1.negJz;
		return Torus(J, freqs, newGF, newptrTM, newE, negJr, negJz);
	}
		
	//interpolate between 2 tori
	EXP Torus TorusGenerator::interpTorus(const double x, const Torus& T0, const Torus& T1) const{
		Torus T0n=T0, T1n=T1;
		if(T0.ptrTM->getToyMapType() != T1.ptrTM->getToyMapType() ||
		   T0.ptrTM->getToyMapType() == ToyPotType::None){
			printf("Different map types\n");
			if(T0.ptrTM->getToyMapType() != ToyPotType::HO){
				printf("fitting HO to T0\n"); 
				T0n = fitTorus(T0.J,1,ToyPotType::HO);
			}
			if(T1.ptrTM->getToyMapType() != ToyPotType::HO){
				printf("fitting HO to T1\n"); 
				T1n = fitTorus(T1.J,1,ToyPotType::HO);
			}
		}
		return InterpTorus(x, T0n, T1n);
	}
	EXP Torus TorusGenerator::interpTorus(const double x, const std::vector<double>& xs,
					      const std::vector<Torus>& Tgrid) const{
		int bot, top;
		double f = bot_top(x, xs, bot, top);
		return interpTorus(f, Tgrid[bot], Tgrid[top]);
	}
	int TorusGrid1::botX(const double x) const {
		double bot = 0, top = nx - 1;
		while (top - bot > 1) {
			int n = .5 * (bot + top);
			if (x > xs[n]) bot = n; else top = n;
		}
		return bot;
	}

	Torus TorusGrid1::T(double x) const {
		int botx;
		if (x < xs.front()) botx = 0;
		else if (x > xs.back()) botx = nx - 2;
		else botx = botX(x);
		const double fx = (x - xs[botx]) / (xs[botx + 1] - xs[botx]);
		return TG.interpTorus(fx, Ts[botx], Ts[botx + 1]);
	}

	int TorusGrid3::botX(const double x) const {
		double bot = 0, top = nx - 1;
		while (top - bot > 1) {
			int n = .5 * (bot + top);
			if (x > xs[n]) bot = n; else top = n;
		}
		return bot;
	}
	int TorusGrid3::botY(const double x) const {
		double bot = 0, top = ny - 1;
		while (top - bot > 1) {
			int n = .5 * (bot + top);
			if (x > ys[n]) bot = n; else top = n;
		}
		return bot;
	}
	int TorusGrid3::botZ(const double x) const {
		double bot = 0, top = nz - 1;
		while (top - bot > 1) {
			int n = .5 * (bot + top);
			if (x > zs[n]) bot = n; else top = n;
		}
		return bot;
	}
	Torus TorusGrid3::T(const double x, const double y, const double z) const {
		int botx, boty, botz;
		if (x < xs.front()) botx = 0;
		else if (x > xs.back()) botx = nx - 2;
		else botx = botX(x);
		if (y < ys.front()) boty = 0;
		else if (y > ys.back()) boty = ny - 2;
		else boty = botY(y);
		if (z < zs.front()) botz = 0;
		else if (z > zs.back()) botz = nz - 2;
		else botz = botZ(z);
		double fx = (x - xs[botx]) / (xs[botx + 1] - xs[botx]);
		double fy = (y - ys[boty]) / (ys[boty + 1] - ys[boty]);
		double fz = (z - zs[botz]) / (zs[botz + 1] - zs[botz]);
		Torus T1(TG.interpTorus(fx, Tn(botx, boty, botz), Tn(botx + 1, boty, botz)));
		Torus T2(TG.interpTorus(fx, Tn(botx, boty + 1, botz), Tn(botx + 1, boty + 1, botz)));
		Torus T3(TG.interpTorus(fx, Tn(botx, boty, botz + 1), Tn(botx + 1, boty, botz + 1)));
		Torus T4(TG.interpTorus(fx, Tn(botx, boty + 1, botz + 1), Tn(botx + 1, boty + 1, botz + 1)));
		Torus T5(TG.interpTorus(fy, T1, T2));
		Torus T6(TG.interpTorus(fy, T3, T4));
		return TG.interpTorus(fz, T5, T6);
	}
	EXP PerturbingHamiltonian interpPerturbingHamiltonian(double x,
		const PerturbingHamiltonian& H0, const PerturbingHamiltonian& H1) {
		if (x == 1) return H0;
		if (x == 0) return H1;
		double xp = 1 - x;
		PerturbingHamiltonian H(H0), H2(H1);
		for (int i = 0; i < H0.indices.size(); i++) {//run over H0 indices
			int mr = H0.indices[i].mr, mz = H0.indices[i].mz, mphi = H0.indices[i].mphi;
			std::vector<std::complex<double> >::iterator jt = H2.values.begin();
			for (GenFncIndices::iterator it = H2.indices.begin(); it != H2.indices.end();) {
				if (mr == (*it).mr && mz == (*it).mz && mphi == (*it).mphi) {//this index matches mine
					H.values[i] = x * H.values[i] + xp * (*jt);
					it = H2.indices.erase(it);
					jt = H2.values.erase(jt);
					break;//my values updated for this index
				}
				else {
					it++; jt++;//move on to the next term in H2
				}
			}
		}
		if (H2.indices.size() > 0) {//something unaccounted for
			for (int j = 0; j < H2.indices.size(); j++) {//add them into my list
				H.indices.push_back(H2.indices[j]);
				H.values.push_back(H2.values[j]);
			}
		}
		return H;
	}

	EXP eTorus TorusGenerator::interpeTorus(const double x, const eTorus& eT0, const eTorus& eT1) const{
		if (x == 1) return eT0;
		if (x == 0) return eT1;
		return eTorus(interpTorus(x, Torus(eT0), Torus(eT1)),
			interpPerturbingHamiltonian(x, eT0.pH, eT1.pH));
	}
	EXP eTorus TorusGenerator::interpeTorus(const double x, const std::vector<double>& xs,
		const std::vector<eTorus>& Tgrid) const{
		int bot, top;
		double f = bot_top(x, xs, bot, top);
		return interpeTorus(f, Tgrid[bot], Tgrid[top]);
	}

	coord::PosMomCyl Torus::from_toy(const Angles& thetaT) const {
		ActionAngles aaT(GF.toyJ(J, thetaT), thetaT);
		return ptrTM->from_aaT(aaT);
	}
	coord::PosMomCyl Torus::from_true(const Angles& theta) const {//input true angles
		ActionAngles aaT(GF.true2toy(ActionAngles(J, theta)));//toy AAs computed from true
		return ptrTM->from_aaT(aaT);
	}
	coord::PosCyl Torus::new_PosDerivJ(const Angles& thetaT,
		DerivAct<coord::Cyl>& dRJ) const {
		ActionAngles aaT(GF.toyJ(J, thetaT), thetaT);
		// Compute derivs wrt JT
		DerivAct<coord::Cyl> dRJT;
		coord::PosMomCyl Rz(ptrTM->from_aaT(aaT, dRJT));
		//dtheta_i/dthetaT_j=dJT_j/dJ_i
		math::Matrix<double> dthetadthetaT(3, 3);
		GF.dtbydtT_Jacobian(thetaT, dthetadthetaT);

		math::Matrix<double> dxdJT(3, 3);
		dxdJT(0, 0) = dRJT.dbyJr.R; dxdJT(1, 0) = dRJT.dbyJr.z; dxdJT(2, 0) = dRJT.dbyJr.phi;
		dxdJT(0, 1) = dRJT.dbyJz.R; dxdJT(1, 1) = dRJT.dbyJz.z; dxdJT(2, 1) = dRJT.dbyJz.phi;
		dxdJT(0, 2) = dRJT.dbyJphi.R; dxdJT(1, 2) = dRJT.dbyJphi.z; dxdJT(2, 2) = dRJT.dbyJphi.phi;
		//dJTdJ_{ij} = dthetadthetaT_{ji}
		math::Matrix<double> Mat3(3, 3);
		math::blas_dgemm(math::CblasNoTrans, math::CblasTrans, 1.0, dxdJT, dthetadthetaT, 0.0, Mat3);

		math::Matrix<double> dpdJT(3, 3);
		dpdJT(0, 0) = dRJT.dbyJr.pR; dpdJT(1, 0) = dRJT.dbyJr.pz; dpdJT(2, 0) = dRJT.dbyJr.pphi;
		dpdJT(0, 1) = dRJT.dbyJz.pR; dpdJT(1, 1) = dRJT.dbyJz.pz; dpdJT(2, 1) = dRJT.dbyJz.pphi;
		dpdJT(0, 2) = dRJT.dbyJphi.pR; dpdJT(1, 2) = dRJT.dbyJphi.pz; dpdJT(2, 2) = dRJT.dbyJphi.pphi;

		math::Matrix<double> Mat4(3, 3);
		math::blas_dgemm(math::CblasNoTrans, math::CblasTrans, 1.0, dpdJT, dthetadthetaT, 0.0, Mat4);

		dRJ.dbyJr.R = Mat3(0, 0); dRJ.dbyJr.z = Mat3(1, 0); dRJ.dbyJr.phi = Mat3(2, 0);
		dRJ.dbyJz.R = Mat3(0, 1); dRJ.dbyJz.z = Mat3(1, 1); dRJ.dbyJz.phi = Mat3(2, 1);
		dRJ.dbyJphi.R = Mat3(0, 2); dRJ.dbyJphi.z = Mat3(1, 2); dRJ.dbyJphi.phi = Mat3(2, 2);
		dRJ.dbyJr.pR = Mat4(0, 0); dRJ.dbyJr.pz = Mat4(1, 0); dRJ.dbyJr.pphi = Mat4(2, 0);
		dRJ.dbyJz.pR = Mat4(0, 1); dRJ.dbyJz.pz = Mat4(1, 1); dRJ.dbyJz.pphi = Mat4(2, 1);
		dRJ.dbyJphi.pR = Mat4(0, 2); dRJ.dbyJphi.pz = Mat4(1, 2); dRJ.dbyJphi.pphi = Mat4(2, 2);
		return Rz;
	}
	/* Position from toy angle plus dR/dthetaT and dpR/dthetaT at fixed J (which causes JT
	 * to vary with thetaT)
	*/
	coord::PosCyl Torus::new_PosDerivs(const Angles& thetaT,
		DerivAngCyl& dRtT, double* det) const {
		ActionAngles aaT(GF.toyJ(J, thetaT), thetaT);
		coord::PosMomCyl Rz(ptrTM->from_aaT(aaT, dRtT));
		// Compute derivs wrt JT
		DerivAct<coord::Cyl> dRJT;
		ptrTM->from_aaT(aaT, dRJT);
		//dtheta_i/dthetaT_j=dJT_j/dJ_i
		math::Matrix<double> dthetadthetaT(3, 3);
		GF.dtbydtT_Jacobian(thetaT, dthetadthetaT);

		DerivAng<coord::Cyl> dJA = GF.dJdt(thetaT); //dJT/dthetaT

		dRtT.dbythetar.R += dRJT.dbyJr.R * dJA.dbythetar.R + dRJT.dbyJz.R * dJA.dbythetar.z;
		dRtT.dbythetaz.R += dRJT.dbyJr.R * dJA.dbythetaz.R + dRJT.dbyJz.R * dJA.dbythetaz.z;
		dRtT.dbythetar.z += dRJT.dbyJr.z * dJA.dbythetar.R + dRJT.dbyJz.z * dJA.dbythetar.z;
		dRtT.dbythetaz.z += dRJT.dbyJr.z * dJA.dbythetaz.R + dRJT.dbyJz.z * dJA.dbythetaz.z;
		dRtT.dbythetar.phi += dRJT.dbyJr.phi * dJA.dbythetar.R
			+ dRJT.dbyJz.phi * dJA.dbythetar.z
			+ dRJT.dbyJphi.phi * dJA.dbythetar.phi;
		dRtT.dbythetaz.phi += dRJT.dbyJr.phi * dJA.dbythetaz.R
			+ dRJT.dbyJz.phi * dJA.dbythetaz.z
			+ dRJT.dbyJphi.phi * dJA.dbythetaz.phi;;

		dRtT.dbythetar.pR += dRJT.dbyJr.pR * dJA.dbythetar.R + dRJT.dbyJz.pR * dJA.dbythetar.z;
		dRtT.dbythetaz.pR += dRJT.dbyJr.pR * dJA.dbythetaz.R + dRJT.dbyJz.pR * dJA.dbythetaz.z;
		dRtT.dbythetar.pz += dRJT.dbyJr.pz * dJA.dbythetar.R + dRJT.dbyJz.pz * dJA.dbythetar.z;
		dRtT.dbythetaz.pz += dRJT.dbyJr.pz * dJA.dbythetaz.R + dRJT.dbyJz.pz * dJA.dbythetaz.z;
		dRtT.dbythetar.pphi = dJA.dbythetar.phi;
		dRtT.dbythetaz.pphi = dJA.dbythetaz.phi;
		dRtT.dbythetaphi.pphi = dJA.dbythetaphi.phi;

		if (det) {//assume only index.mphi=0 non-zero
			(*det) = dRtT.dbythetar.R * dRtT.dbythetaz.z - dRtT.dbythetar.z * dRtT.dbythetaz.R;
		}
		return Rz;
	}

	void Torus::zSoS(std::vector<double>& Rs, std::vector<double>& vRs, const int N,
		double& Rmin, double& Rmax, double& Vmax, const double z0) const {
		const double tol = 1e-5;
		Rmin = 1e10, Rmax = 0, Vmax = 0;
		for (int i = 0; i < N; i++) {
			double thetaT_r = 2 * M_PI / (double)N * (-N / 2 + i);
			CrossingFinder CF(this, thetaT_r, z0);
			double thetaT_z, dtheta = .1, th_min = -.5 * M_PI,
				th_max = th_min + dtheta, z_min, z_max;
			CF.evalDeriv(th_min, &z_min); CF.evalDeriv(th_max, &z_max);
			while (z_min * z_max > 0) {// plod round looking for sign change
				th_min = th_max; z_min = z_max; th_max += dtheta;
				CF.evalDeriv(th_max, &z_max);
			}
			thetaT_z = math::findRoot(CF, th_min, th_max, tol);
			coord::PosMomCyl Rz(from_toy(Angles(thetaT_r, thetaT_z, 0)));
			if (Rz.pz < -.0001) {// keep going round
				do {
					th_min = th_max; z_min = z_max;
					th_max += dtheta; CF.evalDeriv(th_max, &z_max);
				} while (z_min * z_max > 0);
				thetaT_z = math::findRoot(CF, th_min, th_max, tol) - 2 * M_PI;
				Rz = from_toy(Angles(thetaT_r, thetaT_z, 0));
			}
			Rs.push_back(Rz.R); vRs.push_back(Rz.pR);
			Rmin = fmin(Rmin, Rz.R); Rmax = fmax(Rmax, Rz.R); Vmax = fmax(Vmax, fabs(Rz.pR));
		}
		Rs.push_back(Rs[0]); vRs.push_back(vRs[0]);
	}
	void Torus::rSoS(std::vector<double>& thetas, std::vector<double>& pthetas, const double rbar,
		const int N, double& thetamax, double& pmax, double delta) const {
		const double tol = 1e-5;
		thetamax = 0, pmax = 0;
		coord::ProlSph coordsys(delta);
		for (int i = 0; i < N; i++) {
			double thetaT_z = 2 * M_PI / (double)N * (-N / 2 + i);
			r_crossingFinder RCF(this, rbar, thetaT_z, delta);
			double thetaT_r, dtheta = .1, th_min = 0,
				th_max = th_min + dtheta, dr_min, dr_max;
			RCF.evalDeriv(th_min, &dr_min); RCF.evalDeriv(th_max, &dr_max);
			while (dr_min * dr_max > 0 && th_min < 2 * M_PI) {// plod round looking for sign change
				th_min = th_max; dr_min = dr_max; th_max += dtheta;
				RCF.evalDeriv(th_max, &dr_max);
			}
			if (th_min >= 2 * M_PI) continue;
			thetaT_r = math::findRoot(RCF, th_min, th_max, tol);

			coord::PosMomCyl Rv = from_toy(Angles(thetaT_r, thetaT_z, 0));
			coord::PosVelProlSph pps = coord::toPosVel<coord::Cyl, coord::ProlSph>(coord::toPosVelCyl(Rv), coordsys);
			double r = sqrt(pps.lambda - coordsys.Delta2);
			if (pps.nu < 0)pps.nu = -pps.nu;
			double v = (Rv.z == 0) ? .5 * M_PI : atan2(Rv.R / r, Rv.z / sqrt(pps.lambda));
			double pn = Rv.R * Rv.pR / (2 * (pps.nu - coordsys.Delta2)) + Rv.z * Rv.pz / (2 * pps.nu);
			double pv = -2 * pn * cos(v) * sin(v) * coordsys.Delta2;
			coord::PosMomSph rtheta(coord::toPosMomSph(from_toy(Angles(thetaT_r, thetaT_z, 0))));
			if (rtheta.pr < -.0001) {// keep going round
				do {
					th_min = th_max; dr_min = dr_max;
					th_max += dtheta; RCF.evalDeriv(th_max, &dr_max);
				} while (dr_min * dr_max > 0 && th_min < 2 * M_PI);
				thetaT_z = math::findRoot(RCF, th_min, th_max, tol) - 2 * M_PI;
				rtheta = coord::toPosMomSph(from_toy(Angles(thetaT_r, thetaT_z, 0)));
				Rv = from_toy(Angles(thetaT_r, thetaT_z, 0));
				coord::PosVelProlSph pps = coord::toPosVel<coord::Cyl, coord::ProlSph>(coord::toPosVelCyl(Rv), coordsys);
				r = sqrt(pps.lambda - coordsys.Delta2);
				v = (Rv.z == 0) ? .5 * M_PI : atan2(Rv.R / r, Rv.z / sqrt(pps.lambda));
				pn = Rv.R * Rv.pR / (2 * (pps.nu - coordsys.Delta2)) + Rv.z * Rv.pz / (2 * pps.nu);
				pv = -2 * pn * cos(v) * sin(v) * coordsys.Delta2;
			}
			thetas.push_back(v); pthetas.push_back(pv);
			thetamax = fmax(thetamax, v); pmax = fmax(pmax, fabs(pv));
		}
		thetas.push_back(thetas[0]); pthetas.push_back(pthetas[0]);
	}
	void Torus::SoSthetaz(std::vector<double>& X, std::vector<double>& pX,
		const double thetaz,
		const int N, double& Xmax, double& pXmax) const {
		X.clear(); pX.clear();
		Xmax = 0; pXmax = 0;
		for (int i = 0; i < N; i++) {
			double thetar = i * 2 * M_PI / (double)(N - 1);
			Angles thetaT(thetar, thetaz, 0);
			coord::PosMomCyl Rz(from_toy(thetaT));
			ActionAngles aaT(ptrTM->pq2aa(Rz));
			X.push_back(sqrt(aaT.Jr) * cos(aaT.thetar));
			pX.push_back(sqrt(aaT.Jr) * sin(aaT.thetar));
			Xmax = fmax(Xmax, fabs(X.back()));
			pXmax = fmax(pXmax, fabs(pX.back()));
		}
	}
	void Torus::SoSthetar(std::vector<double>& X, std::vector<double>& pX,
		const double thetar,
		const int N, double& Xmax, double& pXmax) const {
		X.clear(); pX.clear();
		Xmax = 0; pXmax = 0;
		for (int i = 0; i < N; i++) {
			double thetaz = i * 2 * M_PI / (double)(N - 1);
			Angles thetaT(thetar, thetaz, 0);
			coord::PosMomCyl Rz(from_toy(thetaT));
			ActionAngles aaT(ptrTM->pq2aa(Rz));
			X.push_back(sqrt(aaT.Jz) * cos(aaT.thetaz));
			pX.push_back(sqrt(aaT.Jz) * sin(aaT.thetaz));
			Xmax = fmax(Xmax, fabs(X.back()));
			pXmax = fmax(pXmax, fabs(pX.back()));
		}
	}

	std::vector<std::pair<coord::PosVelCyl, double> > Torus::orbit(const Angles& theta0, double dt, double T) const {
		std::vector<std::pair<coord::PosVelCyl, double> > traj;
		double t = 0;
		while (t < T) {
			Angles theta(theta0 + (Angles)(freqs * t));
			traj.push_back(std::pair<coord::PosVelCyl, double>(coord::toPosVelCyl(from_true(theta)), t));
			t += dt;
		}
		return traj;
	}

	// Returns true if (R,z,phi) is ever hit by the orbit, and false otherwise. If the 
	// torus passes through the point given, this happens four times, in each case
	// with a different velocity, but only two of these are independent:
	// theta_r -> -theta_r with theta_z -> Pi-theta_z leaves (R,z) and J^T fixed
	// but changes the sign of both velocities. |d(x,y,z)/d(theta_r,theta_z,theta_phi)|
	// is returned. The latter vanishes on the edge of the
	// orbit, such that its inverse, the density of the orbit, diverges there
	// (that's the reason why the density itself is not returned).

	bool Torus::containsPoint(const coord::PosCyl& p, std::vector<Angles>& As,
		std::vector<coord::VelCyl>& Vs,
		std::vector<double>& Jacobs,
		std::vector<DerivAngCyl>* dRdtheta,
		const double tol) const {
		coord::PosMomCyl peri(from_true(Angles(0, .5 * M_PI, 0))), apo(from_true(Angles(M_PI, 0, 0))),
			top(from_true(Angles(M_PI, .5 * M_PI, 0)));
		double Rmin = .95 * peri.R, Rmax = 1.05 * apo.R, zmax = 1.05 * fabs(top.z);
		if (p.R<Rmin || p.R>Rmax || fabs(p.z) > zmax) return false;
		locFinder LF(*this, p);
		double tolerance = 1e-8;
		double params[2] = { 1,1 }, result[2], det;
		int maxNumIter = 200;
		coord::PosMomCyl P1; Angles A1, Atrue;
		DerivAngCyl dA;
		int kmax = 60, nfail = 0;
		while (As.size() < 4 && nfail < kmax) {
			double kount = 0, distsq;
			for (int k = 0; k < kmax; k++) {
				int ntry = math::nonlinearMultiFit(LF, params, tolerance, maxNumIter, result, &distsq);
				if (sqrt(distsq) < 2 * tol) {
					A1 = Angles(math::wrapAngle(result[0]), math::wrapAngle(result[1]), 0.0);
					P1 = from_toy(A1);
					A1.thetaphi += p.phi - P1.phi;
					Atrue = GF.trueA(A1);
					if (is_new(Atrue, As)) break;
				} else nfail++;
				params[0] += .3; params[0] = math::wrapAngle(params[0]);
				params[1] += .7; params[1] = math::wrapAngle(params[1]);
				kount++;
				if (kount == nfail && kount > kmax / 4) return false;
				if (kount > kmax / 4 && As.size() == 2) return true;
				if (kount == kmax && As.size() == 0) return false;
			}
			if (kount < kmax) {
				math::Matrix<double> M(3, 3);//to hold dtheta/dthetaT
				As.push_back(Atrue); new_PosDerivs(A1, dA, &det);
				Vs.push_back(coord::VelCyl(P1.pR, P1.pz, P1.pphi / P1.R));
				Jacobs.push_back(fabs(det / GF.dtbydtT_Jacobian(A1, M)));
				if (dRdtheta) {
					new_assemble(dA, M); dRdtheta->push_back(dA);
				}
				A1.thetar = -A1.thetar; A1.thetaz = M_PI - A1.thetaz; A1.thetaphi = 0;
				P1 = from_toy(A1);
				A1.thetaphi += p.phi - P1.phi;
				Atrue = GF.trueA(A1);
				As.push_back(Atrue); new_PosDerivs(A1, dA, &det);
				Vs.push_back(coord::VelCyl(P1.pR, P1.pz, P1.pphi / P1.R));
				Jacobs.push_back(fabs(det / GF.dtbydtT_Jacobian(A1, M)));
				if (dRdtheta) {
					new_assemble(dA, M); dRdtheta->push_back(dA);
				}
			}
		}
		if (nfail >= kmax) printf("containsPoint error at Rz (%f %f) - %zd angles \n",
			p.R, p.z, As.size());
		return true;
	}
	double Torus::density(const coord::PosCyl& Rz) const {
		std::vector<Angles> As; std::vector<coord::VelCyl> Vs;
		std::vector<double> Jacobs;
		if (!containsPoint(Rz, As, Vs, Jacobs)) return 0;
		double rho = 0;
		for (int i = 0; i < As.size(); i++)
			rho += 1 / Jacobs[i];
		return rho;
	}
	void Torus::write(FILE* ofile) const {
		double val1,val2;
		PtrPointTransform PtrPT=ptrTM->getPointTrans();
		int N=PtrPT->numParams();
		std::vector<double> params(N);
		ptrTM->getParams(&val1,&val2);
		PtrPT->getParams(&params[0]);
		fprintf(ofile, "%g %g %g %g %g %g %g %g %g %g\n",
			J.Jr, J.Jz, J.Jphi, freqs.Omegar, freqs.Omegaz, freqs.Omegaphi,
			E, params[0], val1,val2);
		GF.write(ofile);
	}
	void Torus::read(FILE* ifile){
		double Delta, Js, b;
		fscanf_s(ifile, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
			&J.Jr, &J.Jz, &J.Jphi, &freqs.Omegar, &freqs.Omegaz, &freqs.Omegaphi,
			&E, &Delta, &Js, &b);
		PTIso newPT(Delta);
		Isochrone newIs(Js, b);
		ptrTM=PtrToyMap(new ToyMapIso(newIs,newPT));
		GF.read(ifile);
	}

	/* Gather all terms in the pH that are harmonics of the specified line
	 */
	std::vector<std::complex<double> > PerturbingHamiltonian::get_hn(const GenFncIndex& I,
		std::vector<float>& multiples) const {
		std::vector<std::complex<double> > Hs;
		for (int i = 0; i < indices.size(); i++) {
			GenFncIndex In(indices[i]);
			if ((I.mr == 0 && In.mr != 0) || (I.mr != 0 && In.mr == 0)) continue;
			if ((I.mz == 0 && In.mz != 0) || (I.mz != 0 && In.mz == 0)) continue;
			if ((I.mphi == 0 && In.mphi != 0) || (I.mphi != 0 && In.mphi == 0)) continue;
			//now either matching zero indices or both non-zero
			std::vector<double> Rats;
			if (I.mr != 0) Rats.push_back((double)I.mr / (double)In.mr);
			if (I.mz != 0) Rats.push_back((double)I.mz / (double)In.mz);
			if (I.mphi != 0) Rats.push_back((double)I.mphi / (double)In.mphi);
			if (Rats[1] != Rats[0]) continue;
			if ((Rats.size() == 3) && (Rats[2] != Rats[0])) continue;
			Hs.push_back(values[i]);
			if (Rats[0] >= 1) multiples.push_back(Rats[0]);
			else multiples.push_back(1 / Rats[0]);
		}
		return Hs;
	}
	TMfitter::TMfitter(const potential::BasePotential& pot,
		std::vector<std::pair<coord::PosMomCyl, double> >& _traj,
		double _pphi) : pphi(_pphi),traj(_traj) {
		xmin = 1e6; ymax = 0;
		double ymin = 1e6, xmax = 0, pxmin = 1e6, pxmax = 0, pymin = 1e6, pymax = 0;
		for (int i = 1; i < traj.size(); i++) {//Find where crosses axes
			if (traj[i].first.R * traj[i - 1].first.R <= 0) {//x axis
				double f = traj[i].first.R / (traj[i].first.R - traj[i - 1].first.R);
				double y = (1 - f) * traj[i].first.z + f * traj[i - 1].first.z;
				double px = (1 - f) * traj[i].first.pR + f * traj[i - 1].first.pR;
				ymin = fmin(ymin, fabs(y)); ymax = fmax(ymax, fabs(y));
				pxmin = fmin(pxmin, fabs(px)); pxmax = fmax(pxmax, fabs(px));
			}
			if (traj[i].first.z * traj[i - 1].first.z <= 0) {//z axis
				double f = traj[i].first.z / (traj[i].first.z - traj[i - 1].first.z);
				double x = (1 - f) * traj[i].first.R + f * traj[i - 1].first.R;
				double py = (1 - f) * traj[i].first.pz + f * traj[i - 1].first.pz;
				xmin = fmin(xmin, fabs(x)); xmax = fmax(xmax, fabs(x));
				pymin = fmin(pymin, fabs(py)); pymax = fmax(pymax, fabs(py));
			}
		}
		//xbar, ybar estimated axes of underlying loop orbit
		xbar = .5 * (xmin + xmax); double ybar = .5 * (ymin + ymax);
		Delta2 = (ybar * ybar - xbar * xbar);
		double Phi; coord::GradCyl grad;
		pot.eval(coord::PosCyl(0, ymax, 0), &Phi, &grad);
		Frat = grad.dz;
		pot.eval(coord::PosCyl(xmin, 0, 0), &Phi, &grad);
		Frat /= grad.dR;
		double pxbar = pphi > 0 ? -.5 * (pxmin + pxmax) : .5 * (pxmin + pxmax),
			pybar = pphi > 0 ? .5 * (pymin + pymax) : -.5 * (pymin + pymax);
		aPT = 0.25 * pphi * (1 / (ybar * pybar) + 1 / (xbar * pxbar));
		bPT = 0.125 * pphi * (1 / (ybar * pybar) - 1 / (xbar * pxbar)) - 0.25;
		printf("xbar %f ybar %f ", xbar, ybar);
		printf("pxbar %f pybar %f\n", pxbar, pybar);
		printf("Delta2 %f Frat %f a %f b %f\n", Delta2, Frat, aPT, bPT);
	}
	//Function with root where we have the right force ratio
	double TMfitter::value(double Js) const {
		double g2 = pow((fabs(pphi) + sqrt(pphi * pphi + 4 * Js * Js)) / (2 * Js), 4) - 1;
		double bIso = xbar / sqrt(g2);
		double ap = sqrt(ymax * ymax + bIso * bIso), am = sqrt(xmin * xmin + bIso * bIso);
		//	printf("Js pphi g2 b %f %g %g %g\n",Js,pphi,g2,b);
		return Frat - (ymax * am * pow_2(bIso + am)) / (xmin * ap * pow_2(bIso + ap));
	}
	//Solve for Js and bIso
	std::vector<double> TMfitter::fitTM() const {
		const double Jsmin = .01, Jsmax = 10;
		double Js = math::findRoot(*this, Jsmin, Jsmax, 1e-5);
		double g2 = (pow((fabs(pphi) + sqrt(pphi * pphi + 4 * Js * Js)) / (2 * Js), 4) - 1);
		double bIso = xbar / sqrt(g2);
		std::vector<double> ans(5);
		ans[0] = sqrt(Delta2); ans[1] = aPT; ans[2] = bPT;
		ans[3] = Js; ans[4] = bIso;
		return ans;
	}

	TorusGenerator::TorusGenerator(const potential::BasePotential& _pot,
		const double _tol, std::string _logfname) :
		pot(_pot), defaultTol(_tol),
	    invPhi0(1. / _pot.value(coord::PosCyl(0, 0, 0))),
	    logfname(_logfname), tmax(250) {
		printf("inside TG\n");
		FILE* logfile;
		if (fopen_s(&logfile, logfname.c_str(), "w"))
			printf("I can't open logfile %s\n", logfname.c_str());
		fclose(logfile);
		printf("Preparing TorusGenerator...");
		math::QuinticSpline2d interpEJr;
		//We use mapHJr defined in actions_spherical to create
		//2d splines relating E, Jr, L
		mapHJr(pot,interpEJr,interpJrE);
		printf("done\n");
	}
	double TorusGenerator::Hamilton(const Torus& T, const potential::BasePotential* ePot, const Angles& theta)
	{
		coord::PosMomCyl Rz(T.from_true(theta));
		double H = .5 * (pow_2(Rz.pR) + pow_2(Rz.pz) + pow_2(Rz.pphi / Rz.R)) + pot.value(Rz);
		coord::PosCyl pos(Rz.R, Rz.z, Rz.phi);
		return ePot ? H + ePot->value(pos) : H;
	}
	PerturbingHamiltonian TorusGenerator::get_pH(const Torus& T, int nf, bool ifp,
		const potential::BasePotential* ePot) {//Fourier analyses H
		int nfr = nf, nfz = nf;
		int nfp = ePot ? nf / 4 : 1;
		double N = (nfp * nfr * nfz);
		double* h = new double[nfp * nfr * nfz];
		Angles thetas;
		double dtr = 2 * M_PI / (double)nfr, dtz = 2 * M_PI / (double)nfz, dtp = M_PI / (double)nfp;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) 
#endif
		for (int k = 0; k < nfp; k++) {
			thetas.thetaphi = k * dtp;//if !ePot sticks at 0
			for (int i = 0; i < nfr; i++) {
				thetas.thetar = i * dtr;
				for (int j = 0; j < nfz; j++) {
					thetas.thetaz = j * dtz;
					double a = Hamilton(T, ePot, thetas);
					h[nfr * nfz * k + nfz * i + j] = a;
				}
			}
		}
		//need speq to return cpts at Nyquist vals m_phi, which we won't use 
		math::Matrix<double> speq(nfr, 2 * nfz);
		rlft3(h, speq, nfp, nfr, nfz, 1);
		std::vector<double> Hmods;
		GenFncIndices Hindices;
		int ntop = 0;
		double hmax = 0;
		for (int k = 0; k < nfp; k++) {//find largest perturbing terms
			for (int i = 0; i < nfr; i++) {
				for (int j = 0; j < nfz / 2; j++) {
					double s = sqrt(pow_2(h[nfr * nfz * k + nfz * i + 2 * j])
						+ pow_2(h[nfr * nfz * k + nfz * i + 2 * j + 1])) / N;
					if (ntop<tmax || s>hmax) {
						hmax = insertLine(ntop, tmax, s, GenFncIndex(i, j, k),
							Hmods, Hindices);
					}
				}
			}
		}
		ifp = false;
		std::vector<std::complex<double> > Hvalues;
		for (int i = 0; i < ntop; i++) {
			Hvalues.push_back(std::complex<double>(h[nfr * nfz * Hindices[i].mphi + nfz * Hindices[i].mr + 2 * Hindices[i].mz] / N,
				h[nfr * nfz * Hindices[i].mphi + nfz * Hindices[i].mr + 2 * Hindices[i].mz + 1] / N));
		}
		if (ifp) printf("Terms in perturbing H (*100)\n");
		for (int i = 0; i < ntop; i++) {
			if (Hindices[i].mphi > nfp / 2) Hindices[i].mphi -= nfp;
			if (Hindices[i].mr > nfr / 2)   Hindices[i].mr -= nfr;
			//		Hindices[i].mz*=2;
			//		if(Hindices[i].mz>nfz/2)   Hindices[i].mz   -= nfz;//won't happen!
			if (ifp) printf("(%3d %3d %3d) (%g %g)\n",
				Hindices[i].mr, Hindices[i].mz, Hindices[i].mphi,
				100 * math::modulus(Hvalues[i]), math::arg(Hvalues[i]));
		}
		delete[] h;
		return PerturbingHamiltonian(Hindices, Hvalues);
	}

	void writeClosed(std::vector<coord::PosVelCyl>& shell, Actions J, PtrToyMap& TM, PtrToyMap& baldTM) {
		FILE* ofile; fopen_s(&ofile, "writeClosed.dat", "w");
		int npt = 200;
		fprintf(ofile, "%d\n", npt);
		std::vector<coord::PosMomCyl> Rzr, Rzt;
		for (int i = 0; i < npt; i++) {//draw images of shell isochrone 
			double xr, yr, xt, yt, theta = 2 * M_PI * i / (double)npt;
			coord::PosMomCyl Rzr(TM->from_aaT(J, Angles(M_PI * theta, theta, 0)));
			xr = Rzr.R; yr = Rzr.pR;
			coord::PosMomCyl Rzt(TM->from_aaT(J, Angles(1, theta, 0)));
			xt = Rzt.z; yt = Rzt.pz;
			fprintf(ofile, "%f %f %f %f\n", xr, yr, xt, yt);
		};
		fprintf(ofile, "%zd\n", shell.size());
		for (int i = 0; i < shell.size(); i++) {
			fprintf(ofile, "%f %f %f %f\n",
				shell[i].R, shell[i].vR, shell[i].z, shell[i].vz);
		}
		fclose(ofile);
	}
	PtrToyMap TorusGenerator::chooseTM(GenFncFit& GFFS0,
					   std::vector<double>& params,
					   const Actions& J, double& Jscale,
					   double& freqScale, double& Rsh,
					   ToyPotType ToyMapType, FILE* logfile) const {
		const double L = fabs(J.Jphi) + J.Jz, Xi = J.Jz / L, Xip=fabs(J.Jphi)/L;
		const double Jtot = L + J.Jr;
		Jscale = J.Jr + J.Jz;
		double Delta; pot.getRshDelta(L, Xip, Rsh, Delta);
		// Estimate E
		double Jtot1 = J.Jr + fabs(J.Jphi), Lrel = fabs(J.Jphi)/Jtot1;
		double scaledE = interpJrE.value(log(Jtot1), Lrel);
		double E1=unscaleE(scaledE, invPhi0), E2=potential::E_circ(pot,L);
		double fac = (J.Jz==0)? 0 : exp(-pow(J.Jr / J.Jz, 1.5));//fac -> 1 for shell
		double E = fac*E2 + (1-fac)*E1;
		// Determine if Jphi is low so a HO should be considered
		double Phi0; coord::HessCyl d2Phi;
		pot.eval(coord::PosCyl(0, 0, 0), &Phi0, NULL, &d2Phi);
		math::ScalingSemiInf Sc;
		if (!std::isnan(Phi0) && ToyMapType == ToyPotType::None) {
			double wv2 = 2 * (E-Phi0) + pow_2(J.Jphi / Delta);
			if (!std::isnan(d2Phi.dz2) && d2Phi.dz2 > 0)
				wv2 += d2Phi.dz2 * pow_2(Delta);
			double Jcr = sqrt(wv2) * Delta, k = 0.3;
			double Jphilow = k * Jcr;
			if (J.Jphi < Jphilow) {
				//Choose Is unless Jz<Jzcrit & Jr>0
				if (J.Jr == 0) ToyMapType = ToyPotType::Is;
				else {
					double Jzcrit = pot.getJzcrit(2*J.Jr+J.Jz);//math::evalPoly(JzcritSpl, scale(Sc,2*J.Jr+J.Jz));
					ToyMapType = J.Jz>Jzcrit? ToyPotType::Is :
						ToyPotType::HO;
					if (logfile) {
						if (ToyMapType == ToyPotType::Is)
							fprintf(logfile, "Using Isochrone\n");
						else fprintf(logfile, "Using Harmonic Oscillator\n");
						fprintf(logfile, "box loop orbit transition at this Jr is Jz(Jr)=%f\n", Jzcrit);
					if (fabs(Jzcrit - J.Jz) < 0.2 * Jzcrit)
						fprintf(logfile, "warning:very close to box loop orbit transition. This will lead to inaccurate results\n");
					}
				}
			}
		}
		if(ToyMapType==ToyPotType::None) ToyMapType=ToyPotType::Is;
		double Rmin, Rmax;
		potential::findPlanarOrbitExtent(pot, E1, J.Jphi, Rmin, Rmax);
		freqScale = potential::v_circ(pot, Rsh) / Rsh; //frequency scale set
		int Nn = 5;
		int Nnr = 5;
		double tolerance = 1e-9;//controls optimisation
		PtrToyMap TM;
		if (ToyMapType == ToyPotType::Is) {
			//Now choose isochrone
			//For any Js b=Rsh/ISO.g(Js) f_apo_peri=ISO.f(b.Js,e) & where
			//this matches F_apo_peri in pot we pick Js and b
			double Rs = (1.1*Rmax * (1 - fac) + 2*fac*Rsh);
			Iso ISO(L, J.Jr);
			double Jsmax = 1.1 * Jtot, Jsmin = .1 * Jsmax;
			JsFinder JF(pot, ISO,Rsh);
			double val1, val2;
			JF.evalDeriv(Jsmin, &val1); JF.evalDeriv(Jsmax, &val2);
			while (val1 > 0) {
				Jsmin *= .75; JF.evalDeriv(Jsmin, &val1);
			}
			while (val2 < 0) {
				Jsmax *= 1.5; JF.evalDeriv(Jsmax, &val2);
			}
			const double relToler = 1e-4;
			double Js_iso = math::findRoot(JF, Jsmin, Jsmax, relToler);
			double b_iso = Rsh / ISO.g(Js_iso);
			std::vector<double> params3(Nnr + Nn, 0);
			Isochrone Is(Js_iso,b_iso);
			fitmap fm(Nn, Nnr, J, Is, pot, Rs, Delta);
			if(Nnr+Nn>0)math::nonlinearMultiFit(fm, &params3[0], tolerance, 20, &params3[0]);
			std::vector<double> pr(Nnr), p(Nn);
			for (int i = 0;i < Nnr;i++) {
				pr[i] = params3[i];
				if (logfile) fprintf(logfile, "paramsr: %f\n", pr[i]);
			}
			for (int i = 0;i < Nn;i++){
				p[i] = params3[i + Nnr];
				if (logfile) fprintf(logfile, "paramstheta: %f\n", p[i]);
			}
			if (logfile) fprintf(logfile, "Js, b: %f %f Delta Rs: %f %f\n",
				Js_iso, b_iso, Delta, Rs);
			math::ScalingInfTh sc(Rs);
			PTIso PT(Delta,sc,p,pr);
			TM=PtrToyMap(new ToyMapIso(Is,PT));
		} else {
			// Here we choose HO
			double Rs = 2*(Rmax * (1 - fac) + fac * Rsh);
			//double z1 = sqrt(pow_2(Rsh) + pow_2(Delta));
			double J1 = 2 * J.Jr + fabs(J.Jphi);
			double e = 2 * sqrt(J.Jr * (J.Jr + fabs(J.Jphi))) / J1;
			double omegaR = J1 * (1 + e) * ((1 - fac) / pow_2(Rmax)
				+ fac / pow_2(Rsh));
			double X = 2 / M_PI * Rsh * J.Jz / (J.Jz + fabs(J.Jphi));
			double omegaz = 2 * J.Jz / pow_2(X) * fac + (1 - fac) * omegaR;
			HarmonicOscilattor os(omegaR, omegaz);
			fitmapHarm fm(Nn, Nnr, J, os, pot, Rs, Rs, Delta);
			std::vector<double> params2(Nn + Nnr, 0.0);
			if(Nn+Nnr>0)math::nonlinearMultiFit(fm, &params2[0], tolerance, 20, &params2[0]);
			if (logfile) fprintf(logfile, "HO freqs: %f %f\n", omegaR, omegaz);
			std::vector<double> p(Nnr, 0), pzv(Nn, 0);
			for (int i = 0;i < Nnr;i++) {
				p[i] = params2[i];
				if (logfile) fprintf(logfile, "paramsr: %f\n", p[i]);
			}
			for (int i = 0;i < Nn;i++) {
				pzv[i] = params2[i+Nnr];
				if (logfile) fprintf(logfile, "paramsz: %f\n", pzv[i]);
			}
			math::ScalingInfTh sc(Rs);
			math::ScalingInfTh scz(Rs);
			PTHarm PT(Delta,sc,scz,pzv,p);
			TM=PtrToyMap(new ToyMapHarm(os,PT));
		}
		return TM;
	}
	Torus TorusGenerator::giveBaseTorus(const Actions& J, const PtrToyMap& ptrTM) const {
		std::vector<double> params;
		int nrmax = 0, nzmax = 0;// nzmax must be even
		double Hbar;
		GenFncIndices indices = makeGridIndices(nrmax, nzmax);
		//GenFncFit GFFS(indices, 8, 6, J);
		GenFncFit GFFS(indices, 1, 1, J);
		torusFitter TF(J, pot, 1, ptrTM, GFFS);
		Frequencies freqs;
		GenFncDerivs dPdJ;
		bool negJr, negJz;
		TF.fitAngleMap(&params[0], Hbar, freqs, dPdJ, negJr, negJz);
		GenFnc G(indices, params, dPdJ);
		return Torus(J, freqs, G, ptrTM, Hbar, negJr, negJz);
	}

	Torus TorusGenerator::fitTorus(const Actions& J, const double tighten, const ToyPotType ToyMapType) const {
		FILE* logfile; fopen_s(&logfile, logfname.c_str(), "a");
		fprintf(logfile, "\nActions (%f %f %f)\n", J.Jr, J.Jz, J.Jphi);
		int nrmax = 6, nzmax = 4;// nzmax must be even
		GenFncIndices indices = makeGridIndices(nrmax, nzmax);
		std::vector<double> params(indices.size(), 0);
		double timesr = 1.5, timesz = 1.5;
		GenFncFit GFFS0(indices, timesr, timesz, J);
		double Jscale, freqScale, Rsh;
		PtrToyMap ptrTM(chooseTM(GFFS0, params, J, Jscale, freqScale, Rsh,
					 ToyMapType, logfile));
		double tolerance = 1e-9;//controls optimisation of the given Sn
		double tol = defaultTol * tighten;
		double Hbar, Hdisp = 1e20, Htarget = tol * freqScale * Jscale;
		bool converged = false;
		int Loop = 0, MaxLoop = 12, maxNumIter = 10;
		bool stuck = false;
		std::vector<double> rep;
		do {
			double Hdisp_old = Hdisp;
			int numTerms_old = GFFS0.numTerms();
			GenFncFit GFFS(indices, timesr, timesz, J);
			torusFitter TF(J, pot, freqScale, ptrTM, GFFS);
			if (std::isnan(Hdisp)) exit(0);
			try {
				NANbar = 0;
				int numIter = math::nonlinearMultiFit(TF, &params[0], tolerance, maxNumIter, &params[0], &Hdisp);
				Hdisp = sqrt(Hdisp);
				fprintf(logfile, "%d terms %d points -> dH %7.3e",
					GFFS.numTerms(), GFFS.numPoints(), Hdisp);
				rep.push_back(Hdisp);
				if (NANfrac + NANbar > 0)
					fprintf(logfile, " NANfrac %f NANbar %f\n", NANfrac, NANbar);
				else
					fprintf(logfile, "\n");
				converged = (Hdisp < tol * freqScale * Jscale);
			}
			catch (std::exception& e) {
				std::cout << "Exception in fitTorus: " << e.what() << '\n';
			}
			if(Loop>8)
				stuck = fabs(1 - rep[Loop] / rep[Loop - 1]) < 1e-2 && GFFS.numTerms() > numTerms_old;
			if (converged || (Loop > 1 && stuck)) break;
			if (GFFS.numTerms() > numTerms_old && Hdisp > Hdisp_old) {//something wrong: back up 
				indices = GFFS0.indices;
				params.resize(indices.size());
				Hdisp = Hdisp_old;
				timesr *= 1.2; timesz *= 1.2;
			}
			else {
				GFFS0 = GFFS;
				indices = GFFS.expand(params);
			}
			Loop++;
		} while (Loop < MaxLoop);
		if(stuck){
			printf("Quit early because Hvar not decreasing\n");
			for(int i=0;i<Loop;i++) printf("%g ",rep[i]);
		}
		GenFncFit GFFS(indices, 2*timesr, 2*timesz, J);
		printf("Hdisp %e after %d expansions\n", Hdisp, Loop);
		if (!converged) {
			fprintf(logfile, "\nfitTorus failed to converge: %7.3e vs %7.3e target. ",
				Hdisp, Htarget);
			fprintf(logfile, "%zd terms", indices.size());
			fprintf(logfile, "\nNANfracs: %f %f, resids:\n", NANfrac, NANbar);
			for (int i = 0; i < rep.size(); i++) fprintf(logfile, "%7.2e ", rep[i]); printf("\n");
		}
		torusFitter TF(J, pot, freqScale, ptrTM, GFFS);
		Frequencies freqs;
		GenFncDerivs dPdJ;
		bool negJr,negJz;
		Hdisp = TF.fitAngleMap(&params[0], Hbar, freqs, dPdJ,negJr,negJz);
		params.resize(indices.size()); dPdJ.resize(indices.size());
		GenFnc G(indices, params, dPdJ);
		fclose(logfile);
		return Torus(J, freqs, G, ptrTM, Hbar, negJr, negJz);
	}

	Torus TorusGenerator::fitBaseTorus(const Actions& J, const double tighten, const ToyPotType ToyMapType) const {
		FILE* logfile; fopen_s(&logfile, logfname.c_str(), "a");
		fprintf(logfile, "\nActions (%f %f %f)\n", J.Jr, J.Jz, J.Jphi);
		int nrmax = 2, nzmax = 6;// nzmax must be even
		GenFncIndices indices = makeGridIndices(nrmax, nzmax);
		std::vector<double> params(indices.size(), 0);
//		GenFncFit GFFS0(indices, 8, 6, J);
		GenFncFit GFFS0(indices, 1, 1, J);
		double Jscale, freqScale, Rsh;
		PtrToyMap ptrTM(chooseTM(GFFS0, params, J, Jscale, freqScale, Rsh,ToyMapType));
		fclose(logfile);
		return giveBaseTorus(J, ptrTM);
	}

	Torus TorusGenerator::fitFrom(const Actions& J, const Torus& T, const double tighten) const {
		GenFncIndices indices(T.GF.indices);
		std::vector<double> params(T.GF.values);
		double freqScale = T.freqs.Omegaz;
		PtrToyMap ptrTM(T.ptrTM);
		double tolerance = 1e-9;//controls optimisation of the given Sn
		double Hbar, Hdisp = 1e20;
		int maxNumIter = 10;
		std::vector<double> rep;
		double timesr=1.6, timesz=1.6;
		GenFncFit GFF(indices, timesr, timesz, J);
		torusFitter TF(J, pot, freqScale, ptrTM, GFF);
		try {
			int numIter = math::nonlinearMultiFit(TF, &params[0], tolerance, maxNumIter, &params[0], &Hdisp);
			Hdisp = sqrt(Hdisp);
			rep.push_back(Hdisp);
		}
		catch (std::exception& e) {
			std::cout << "Exception in fitFrom: " << e.what() << '\n';
		}
		Frequencies freqs;
		GenFncDerivs dPdJ;
		bool negJr, negJz;
		GenFncFit GFFS(indices, 2*timesr, 2*timesz, J);
		torusFitter TFS(J, pot, freqScale, ptrTM, GFFS);
		Hdisp = TFS.fitAngleMap(&params[0], Hbar, freqs, dPdJ, negJr, negJz);
		if(negJr || negJz) printf("actions<0 in fitAngleMap: %d %d\n",negJr,negJz);
		GenFnc G(indices, params, dPdJ);
		return Torus(J, freqs, G, ptrTM, Hbar, negJr, negJz);
	}
	eTorus TorusGenerator::fiteTorus(const Actions& J, const double tighten,
		const potential::BasePotential* ePhi) {
		Torus T = fitTorus(J, tighten);
		int nf = 128;
		PerturbingHamiltonian pH(get_pH(T, nf, true, ePhi));
		return eTorus(T, pH);
	}
	eTorus TorusGenerator::fiteTorus(const Actions& J, const potential::BasePotential* ePhi) {
		Torus T = fitTorus(J);
		int nf = 128;
		PerturbingHamiltonian pH(get_pH(T, nf, true, ePhi));
		return eTorus(T, pH);
	}
	double TorusGenerator::getRsh(const Actions& J) const {
		const double L = J.Jr + J.Jz + fabs(J.Jphi), Xip = fabs(J.Jphi)/L;
		double Rsh, Delta;
		pot.getRshDelta(L, Xip, Rsh,Delta);
		return Rsh;
//		return interpR.value(L, (J.Jr + J.Jz) / L);
	}
	double TorusGenerator::getDelta(const Actions& J) const {
		double Rsh, Delta, L = fabs(J.Jphi) + J.Jz, Xip = fabs(J.Jphi) / L;
		pot.getRshDelta(L, Xip, Rsh, Delta);
		return Delta;
	}
	std::vector<Torus> TorusGenerator::constE(const double Jrmin, const Actions& Jstart, const int Nsteps) {
		double fac = exp(-log(Jstart.Jr / Jrmin) / (Nsteps - 1));
		Actions Jnext(Jstart);
		std::vector<Torus> Tgrd;
		Torus T(fitTorus(Jnext));
		Tgrd.push_back(T);
		for (int i = 1; i < Nsteps; i++) {
			double frat = Tgrd.back().freqs.Omegar / Tgrd.back().freqs.Omegaz;
			double dJr = (1 - fac) * Jnext.Jr;
			Jnext.Jr -= dJr; Jnext.Jz += dJr * frat;
			Torus T1(fitTorus(Jnext));
			Torus T2 = interpTorus(0.5, Tgrd.back(), T1);
			double frat1 = T2.freqs.Omegar / T2.freqs.Omegaz;
			Jnext.Jz = Jnext.Jz + dJr * frat - dJr * frat1;
			printf("Next Jr/Jz: %f Rsh: %f\n", Jnext.Jr / Jnext.Jz, getRsh(Jnext));
			Tgrd.push_back(fitTorus(Jnext));
		}
		return Tgrd;
	}

	EXP ActionAngles ActionFinderTG::actionAngles(const coord::PosVelCyl& point, Frequencies* freq) const {
		Torus T;
		ActionAngles aa(actionAnglesTorus(point, T));
		if (freq) *freq = T.freqs;
		return aa;
	}
	EXP ActionAngles ActionFinderTG::actionAnglesTorus(
		const coord::PosVelCyl& point, Torus& T) const {
		double E = potential::totalEnergy(*pot, point);
		if (E > 0) {
			printf("Energy is: %f and positive so no orbit\n", E);
			return ActionAngles(Actions(NAN, NAN, NAN), Angles(NAN, NAN, NAN));
		}
		const double tol = 1e-5;
		double phi0 = point.phi;
		while (fabs(phi0) > M_PI) phi0 += phi0 > M_PI ? -M_PI : M_PI;
		coord::PosMomCyl P0(coord::toPosMomCyl(point));
		Angles trueA;
		Actions J(AF->actions(point));
		if (std::isnan(J.Jr) || std::isnan(J.Jz)) {
			printf("Fudge has returned a NaN: Jr,Jz,Jphi (%f %f %f)\n",
				J.Jr, J.Jz, J.Jphi);
			exit(0);
		}
		T = TG.fitTorus(J);
		PtrPointTransform PtrPT = T.ptrTM->getPointTrans();
		std::vector<double> dJt_old;
		ActionAngles aaT = T.ptrTM->pq2aa(P0);
		Angles tT(aaT);	Actions JT(aaT);
		printf("thetaT: %f %f %f\n",tT.thetar,tT.thetaz,tT.thetaphi);
		double diff = 1e20;
		int kount = 0;
		while (kount < 10) {
			Torus T1 = kount>0? TG.fitFrom(J, T) : T;
			Actions JT1 = T1.GF.toyJ(J, tT);
			std::vector<double> df = {JT.Jr - JT1.Jr, JT.Jz - JT1.Jz, 0.0};
			math::Matrix<double> dthetadthetaT(3, 3);
			T1.GF.dtbydtT_Jacobian(tT, dthetadthetaT);
			math::LUDecomp LUM(dthetadthetaT);
			math::Matrix<double> inv = LUM.inverse(3);
			std::vector<double> dJt(3, 0.0);
			math::blas_dgemv(math::CblasTrans, 1.0, inv, df, 0.0, dJt);
			double old_diff = diff;
			diff = sqrt(pow_2(dJt[0]) + pow_2(dJt[1]));
			if (diff >= old_diff) {//We've gone backwards!
				J.Jr = T.J.Jr; J.Jz = T.J.Jz;
				trueA = T.GF.trueA(tT);// T=T1;
				break;
			}
			T=T1;
			dJt_old = dJt;
			J.Jr += dJt[0]; J.Jz += dJt[1];
			printf("kount:%d (%f,%f) %f\n", kount, J.Jr, J.Jz, diff);
			if (diff < tol && kount>0) {
				trueA = T.GF.trueA(tT);	break;
			}
			kount++;
		}
		return ActionAngles(J, trueA);
	}


}//namespace
