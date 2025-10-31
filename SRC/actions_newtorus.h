#pragma once
#include "map.h"
#include "actions_base.h"
#include "actions_staeckel.h"
#include "actions_newgenfnc.h"
#include "actions_newisochrone.h"
#include "actions_harmonicoscillator.h"
#include "actions_focal_distance_finder.h"
#include "actions_spherical.h"
#include "potential_base.h"
#include "potential_utils.h"
#include "math_core.h"
#include "math_specfunc.h"
#include "math_fit.h"
#include "math_linalg.h"
#include "math_fourier.h"
#include "coord.h"
#include "particles_base.h"
#include <iostream>
#include <map>
#include <complex>
#include <cassert>


namespace actions {
	/* Class to hold the Fourier decomposition of the residual Hamiltonian
	 * after torus fitting
	*/
	class EXP PerturbingHamiltonian {
	private:
	public:
		GenFncIndices indices;
		std::vector<std::complex<double> > values;
		PerturbingHamiltonian(const GenFncIndices& _indices,
			const std::vector<std::complex<double> >& _values) :
			indices(_indices), values(_values) {
		}
		GenFncIndex index(const int i) const {
			return indices[i];
		}
		std::complex<double> value(const int i) const {
			return values[i];
		}
		int numTerms() {
			return indices.size();
		}
		std::vector<std::complex<double> > get_hn(const GenFncIndex&, std::vector<float>&) const;
	};
	EXP PerturbingHamiltonian interpPerturbingHamiltonian(const double x,
		const PerturbingHamiltonian& H0, const PerturbingHamiltonian& H1);

	/*
	 * Base of all torus classes
	 */
	class EXP Torus {
	public:
		Actions J;
		Frequencies freqs;
		GenFnc GF;
		PtrToyMap ptrTM;
		double E;
		bool negJr, negJz;
		Torus() {}
		/* Creator called by TorusGenerator rather than users */
		Torus(const Actions& _J, const Frequencies& _freqs, const GenFnc& _GF,
			const PtrToyMap _ptrTM, double _E, bool _negJr, bool _negJz) :
		    J(_J), freqs(_freqs), GF(_GF), ptrTM(_ptrTM), E(_E), negJr(_negJr), negJz(_negJz) {
		}
		Torus(const Torus& T) :
		    J(T.J), freqs(T.freqs), GF(T.GF), ptrTM(T.ptrTM), E(T.E), negJr(T.negJr), negJz(T.negJz) {
		}
		void printGF(void) {
			GF.print();
		}
		coord::PosMomCyl from_toy(const Angles&) const;
		coord::PosMomCyl from_true(const Angles&) const;
		coord::PosCyl new_PosDerivJ(const Angles&, actions::DerivAct<coord::Cyl>&) const;
		coord::PosCyl new_PosDerivs(const Angles&, actions::DerivAng<coord::Cyl>&,
			double* det = NULL) const;/*dR/dtheta etx*/
		/* compute the surface of section z=0 pz>0 */
		void zSoS(std::vector<double>& R, std::vector<double>& vR, const int N,
			double& Rmin, double& Rmax, double& Vmax, const double z0 = 0) const;
		/* compute the surface of section R=Rbar pR>0 */
		void rSoS(std::vector<double>& z, std::vector<double>& vz, const double Rbar, const int N,
			double& zmax, double& Vmax, double delta = 1e-3) const;
		void SoSthetaz(std::vector<double>& X, std::vector<double>& pX, const double thetaz,
			const int N, double& Xmax, double& pXmax) const;
		void SoSthetar(std::vector<double>& X, std::vector<double>& pX, const double thetar,
			const int N, double& Xmax, double& pXmax) const;
		Frequencies Omega(void) const {
			return freqs;
		}
		/* Compute the orbit from given angles using computed frequencies */
		std::vector<std::pair<coord::PosVelCyl, double> > orbit(const Angles& theta0, double dt, double T) const;
		/* Does the torus pass through given point? If so for
		 * what true thetas? */
		bool containsPoint(const coord::PosCyl& pt, std::vector<Angles>& As,
			std::vector<coord::VelCyl>& Vs,
			std::vector<double>& Jacobs,
			std::vector<actions::DerivAngCyl>* dA = NULL,
			const double tol = 1e-5) const;
		/* returns density contributed at location */
		double density(const coord::PosCyl&) const;
		void write(FILE*) const;
		void read(FILE*);
	};

	/* The class of tori that include perturbing Hamiltonians. */
	class EXP eTorus : public Torus {
	private:
	public:
		PerturbingHamiltonian pH;
		eTorus(const Actions& _J, const Frequencies& _freqs, const GenFnc& _GF,
			const PtrToyMap _ptrTM, double _E, bool _negJr, bool _negJz,
			const PerturbingHamiltonian& _pH) :
			Torus(_J, _freqs, _GF, _ptrTM, _E, _negJr, _negJz), pH(_pH) {
		}
		eTorus(const Torus& T, const PerturbingHamiltonian& _pH) :
			Torus(T), pH(_pH) {
			//printf("eTorus at E = %f created: %d terms in pH\n",T.E, pH.numTerms());
		}
		PerturbingHamiltonian Hns() const {
			return pH;
		}
		std::vector<std::complex<double> > get_hn(const GenFncIndex& Indx, std::vector<float>& multiples) const {
			return pH.get_hn(Indx, multiples);
		}
	};

	/*
	 * Class for fitting torus to an orbit
	*/
	class EXP TMfitter : public math::IFunctionNoDeriv {
	private:
		double pphi, Delta2, xmin, ymax, xbar, Frat, aPT, bPT;
		std::vector<std::pair<coord::PosMomCyl, double> >& traj;
	public:
		TMfitter(const potential::BasePotential&,
			std::vector<std::pair<coord::PosMomCyl, double> >&, double);
		std::vector<double> fitTM() const;
		virtual double value(double) const;
		virtual unsigned int numVars() const { return 4; }
		virtual unsigned int numValues() const { return 1; }
	};
	/*
	 * Class for generators of tori.
	*/
	class EXP TorusGenerator {
	private:
		const potential::BasePotential& pot;
		const double defaultTol, invPhi0;
		math::LinearInterpolator2d interpD; //for Delta values
		math::LinearInterpolator2d interpR; //for Rshell(L,Xi) values
		math::LinearInterpolator2d interpV; //for Vshell(L,Xi) values
		math::LinearInterpolator2d interpRE; //for Rshell(E,Xi) values
		math::LinearInterpolator interpJz;//Jr(Jr) for box loop orbit transition
		math::QuinticSpline2d interpJrE;//Esc(Q,Y) for planar orbit, Q=log(Lz+Jr),Y=Lz/(Lz+Jr). Esc=log(1/Phi(0)-1/E)
		std::string logfname;
#ifdef TEST
		/* Test_it compares analytic and numerical derivatives */
		void test_it(actions::Actions& J, std::vector<double>&);
#endif
		/* Hamilton computes residual H for Phi+ePhi where
		 * ePhi is an optional additional potential not used in fitting */
		double Hamilton(const Torus&, const potential::BasePotential*, const Angles&);
		/* PerturbingHamiltonian computes Fourier decomposition
		 * of the residual H */
		PerturbingHamiltonian get_pH(const Torus&,
			int nf, bool ifp, const potential::BasePotential*);
		//		void setConsts(actions::Actions, double, double&, double&, double&, Isochrone&, coord::UVSph&) const;
				//ToyMap chooseTM(actions::Actions, double&, double&, double&) const;
		int tmax;// Max number of terms retained in residual H
	public:
		/* Creator of tori in given potential. GF deemed ok if
		 * dispersion in H < tol*freqScale*Jtotal */
		TorusGenerator(const potential::BasePotential& _pot,
			const double _tol = 1e-9, std::string _logfname = "TG.log");
		PtrToyMap chooseTM(GenFncFit&, std::vector<double>&, const Actions&,
			double&, double&, double&, ToyPotType ToyMapType=ToyPotType::None, FILE* logfile = NULL) const;
		Torus fitTorus(const Actions& J, const double tighten = 1, const ToyPotType ToyMapType=ToyPotType::None) const;
		/* Fit a torus with all Sn=0 */
		Torus fitBaseTorus(const Actions&, const double tighten = 1, const ToyPotType ToyMapType=ToyPotType::None) const;
		/* Build a torus with all Sn=0 around the given ToyMap */
		Torus giveBaseTorus(const Actions&, const PtrToyMap&) const;
		/* fitFrom uses given TM & varies Sn starting from given values
		*/
		Torus fitFrom(const Actions&, const Torus&, const double tighten = 1) const;
		eTorus fiteTorus(const Actions&, const potential::BasePotential* _addPhi = NULL);
		eTorus fiteTorus(const Actions&, const double, const potential::BasePotential* _addPhi = NULL);
		double getDelta(double&, double&);
		double getDelta(Actions&);
		double getRsh(Actions&);
		void test_it(const Actions&, std::vector<double>&);
		std::vector<Torus> constE(const double Jrmin, const Actions& Jstart, const int Nstep);
		void old_getHn(const Torus&, int);
		//interpolate between 2 tori
		Torus interpTorus(const double x, const Torus& T0, const Torus& T1) const;
		//interpolae on an indexed array of tori 
		Torus interpTorus(const double x, const std::vector<double>&, const std::vector<Torus>&) const;
		//interpolate between 2 tori
		eTorus interpeTorus(const double x, const eTorus& T0, const eTorus& T1) const;
		//interpolate on an indexed array of tori 
		eTorus interpeTorus(const double x, const std::vector<double>&, const std::vector<eTorus>&) const;

	};

	class EXP TorusGrid1 {
	private:
		const std::vector<double>& xs;
		const::std::vector<Torus> &Ts;
		const int nx;
		const TorusGenerator& TG;
		int botX(const double) const;
	public:
		TorusGrid1(std::vector<double>& _xs,
			std::vector<actions::Torus> &_Ts, const TorusGenerator& _TG) :
			xs(_xs), Ts(_Ts), nx(xs.size()), TG(_TG) {
		}
		Torus T(const double x) const;
	};
	class EXP TorusGrid3 {
	private:
		const std::vector<double>& xs, ys, zs;
		const std::vector<Torus>& Ts;
		const int nx, ny, nz;
		const TorusGenerator& TG;
		int botX(const double) const;
		int botY(const double) const;
		int botZ(const double) const;
		Torus Tn(int ix, int iy, int iz) const {
			return Ts[(ix * ny + iy) * ny + iz];
		}
	public:
		TorusGrid3(std::vector<double>& _xs, std::vector<double>& _ys,
			std::vector<double>& _zs, std::vector<actions::Torus>& _Ts, const TorusGenerator& _TG) :
			xs(_xs), ys(_ys), zs(_zs), Ts(_Ts), nx(xs.size()), ny(ys.size()), nz(zs.size()), TG(_TG) {
		}
		Torus T(const double x, const double y, const double z) const;
		Torus T(const Actions J) const {
			return T(J.Jr, J.Jz, J.Jphi);
		}
	};
	class EXP ActionFinderTG : public BaseActionFinder {
	private:
		const potential::PtrPotential pot;
		const TorusGenerator& TG;
		const ActionFinderAxisymFudge* AF;

	public:
		ActionFinderTG(const potential::PtrPotential& _pot,
			const TorusGenerator& _TG,
			const ActionFinderAxisymFudge* _AF) :
			pot(_pot), TG(_TG), AF(_AF) {
		}
		ActionFinderTG(const potential::PtrPotential& _pot,
			const TorusGenerator& _TG) :
			pot(_pot), TG(_TG), AF(new ActionFinderAxisymFudge(_pot)) {
		}
		//virtual Actions actions(const coord::PosMomCyl& point) const;
		virtual ActionAngles actionAngles(const coord::PosVelCyl& point,
			Frequencies* freq = NULL) const;
		virtual ActionAngles actionAngles(const coord::PosVelCar(point)) const {
			return actionAngles(coord::toPosVelCyl(point));
		}
		virtual ActionAngles actionAnglesTorus(const coord::PosVelCyl& point,
			Torus& T) const;
		virtual ActionAngles actionAngles(const coord::PosCar& x, const coord::VelCar& v)  const {
			return actionAngles(coord::PosVelCar(x, v));
		}
		virtual Actions actions(const coord::PosVelCyl& point) const {
			return Actions(actionAngles(point));
		}
	};
	EXP void getGridBoxLoop(const potential::BasePotential& pot,
				std::vector<double>& gridE,
				std::vector<double>& gridJr,
				std::vector<double>& gridJz);
	EXP std::vector<double> mapJcrit(const potential::BasePotential& pot);
}//namespace actions