#include "potential_interpolators.h"

/// number of sampling points for a shell orbit (equally spaced in time)
static const unsigned int NUM_STEPS_TRAJ = 64;
/// accuracy of root-finding for the radius of thin (shell) orbit
static const double ACCURACY_RSHELL = 1e-6;
/// accuracy of orbit integration for shell orbit
static const double ACCURACY_INTEGR = 1e-8;
/// accuracy parameter determining the spacing of the interpolation grid along the energy axis
static const double ACCURACY_INTERP2 = 1e-5;
/// upper limit on the number of timesteps in ODE solver (should be enough to track half of the orbit)
static const unsigned int MAX_NUM_STEPS_ODE = 2000;

namespace potential{

namespace{

// Class to find the largest FD that yields a centrifugal barrier even
// with Lz=0. Rsh is the radius of the shell orbit and vR is vR is
// what's required to create the orbit of interest shooting from
// (Rsh,0). Delta is an estimate of FD. The method bestFD() returns the
// FD with its 1st & 2nd argunents the value of u at which
// 0=p_u^2=dp_u^2/du and the value there of d^2p_u^2/du^2
//
class EXP FDfinder{
	private:
		double E,Rsh,vR,P0,Delta0,EmP0;
		const potential::BasePotential& pot;
	public:
		FDfinder(const double _E,const double _Rsh,const double _vR,
			 const double _Delta,
			 const potential::BasePotential& _pot) :
		    E(_E), Rsh(_Rsh), vR(_vR), Delta0(_Delta), pot(_pot){
			P0=pot.value(coord::PosCyl(Rsh,0,0));
			EmP0=E-P0;
		}
		math::Matrix<double> derivs(double u,double Delta,double& d2p2du2,
					    double* p2=NULL,double* p2prime=NULL);
		double bestFD(double&, double&);
};

EXP math::Matrix<double> FDfinder::derivs(double u,double Delta, double& d2p2du2,
					  double* p2,double* p2prime){
	double sh=sinh(u), ch=cosh(u), u0=asinh(Rsh/Delta);
	double R=Delta*sh, R2=R*R, Delta2=pow_2(Delta);
	double P;
	coord::GradCyl gradP;
	coord::HessCyl hessP;
	pot.eval(coord::PosCyl(R,0,0),&P,&gradP,&hessP);
	double EmP=E-P;
	double X=2*(EmP*R2-Delta2*P);
	double dXdR=4*EmP*R-2*(R2+Delta2)*gradP.dR;
	double dXdDelta=-4*Delta*P;
	double d2XdR2=4*EmP-8*R*gradP.dR-2*(R2+Delta2)*hessP.dR2;
	double d2XdDeltadR=-4*Delta*gradP.dR;
	double pu0=Delta*cosh(u0)*vR;
	if(p2!=NULL) *p2=pow_2(pu0)+X-2*(EmP0*pow_2(Rsh)-Delta2*P0);
	double dp2du=dXdR*Delta*ch;
	d2p2du2 = d2XdR2*pow_2(Delta*ch)+dXdR*Delta*sh;
	if(p2prime!=NULL) *p2prime=dp2du;
	double dp2dDelta=2*pu0*vR/cosh(u0) + dXdR*sh + dXdDelta + 4*Delta*P0;
	double dp2primedu=d2XdR2*pow_2(Delta*ch)+dXdR*Delta*sh;
	double dp2primedDelta=d2XdR2*Delta*ch*sh+dXdR*ch+d2XdDeltadR*Delta*ch;
	math::Matrix<double> M(2,2);
	M(0,0)=dp2du;      M(0,1)=dp2dDelta;
	M(1,0)=dp2primedu; M(1,1)=dp2primedDelta;
	return M;
}

EXP double FDfinder::bestFD(double& umin, double& d2p2du2){//implements N-R search for p_u^2=dp_u^2/du=0
	double u=asinh(.8*Rsh/Delta0), u0=asinh(Rsh/Delta0), pu0=vR*Delta0*cosh(u0);
	double p2=1, p2prime=1, Delta=Delta0, fac=1;
	int i=0;
	while(i<20 && u>0 && (fabs(p2)>1e-4 || fabs(p2prime)>1e-4)){
		math::Matrix<double> M0(derivs(u,Delta,d2p2du2,&p2,&p2prime));
		double det=M0(0,0)*M0(1,1)-M0(1,0)*M0(0,1);
		double dY0 = ( M0(1,1)*p2-M0(0,1)*p2prime)/det;
		double dY1 = (-M0(1,0)*p2+M0(0,0)*p2prime)/det;
		//printf("%d %g %g %g %g %g\n",i,det,p2,p2prime,u,Delta);
		while(fabs(dY0)>.5*u){
			dY0*=.5; dY1*=.5;
		}
		while(fabs(dY1)>.5*Delta){
			dY0*=.5; dY1*=.5;
		}
		u-=dY0; Delta-=dY1;
		i++;
	}
	umin=u;
	return Delta;
}	

/* The next three classs are used to compute Jr & Jz at the box/loop
 * transition.
 * pzf.evalDeriv returns pz for motion along the z axis
 */
class pzf : public math::IFunction {
	private:
		const potential::BasePotential& pot;
		const double E;
	public:
		pzf(const potential::BasePotential& _pot,double _E) : pot(_pot), E(_E) {};
		virtual void evalDeriv(const double z,
				       double* value = NULL, double* deriv = NULL,
				       double* deriv2 = NULL) const {
			coord::PosCar X(0, 0, z);
			double Phi;
			coord::GradCar dx;
			pot.eval(X, &Phi, &dx);
				//std::cout << "Value:" << x << " " << pot.value(X) << "\n";
			double p2 = 2 * (E - Phi);
			if (p2 <= 0) {
				if (value)*value = 0;
			}
			else {
				double p = sqrt(p2);
				if (value)*value = p;
				if (deriv)*deriv = -dx.dz / p;
			}
		};
		virtual unsigned int numDerivs() const { return 1; }
};
/* Now some code used to compute (Jr,Jz) at the box/loop transition
 * when Jphi=0. We compute Jf=2*Jr+Jz by integrating along the
 * (unstable) short-axis orbit and then get Jr from the area of the SoS
 * when the instability of the orbit is triggered
 */
//pxf.evalDeriv returns px for motion along x axis
class pxf : public math::IFunction {
	private:
		const potential::BasePotential& pot;
		const double E;
	public:
		pxf(const potential::BasePotential& _pot, double _E) :pot(_pot), E(_E){};
		virtual void evalDeriv(const double x,
				       double* value = NULL, double* deriv = NULL,
				       double* deriv2 = NULL) const {
			coord::PosCar X(x, 0, 0);
			double Phi;
			coord::GradCar dx;
			pot.eval(X, &Phi, &dx);
				//std::cout << "Value:" << x << " " << pot.value(X) << "\n";
			double p2 = 2 * (E - Phi);
			if (p2 <= 0) {
				if (value) *value = 0;
			}
			else {
				double p = sqrt(p2);
				if (value) *value = p;
				if (deriv) *deriv = -dx.dx / p;
			}
		};
		virtual unsigned int numDerivs() const { return 1; }
};
//finds the z1 s.t. Iz=int pzdz between 0 and z1, pz is the z momentum at that energy.
class Ifind : public math::IFunction {
	private:
		const double I;
		pzf pzfunc;
	public:
		Ifind(double _Iz, const potential::BasePotential& _pot, double _E)
				:I(_Iz), pzfunc(_pot,_E){
		};
		virtual void evalDeriv(const double z,
				       double* value = 0, double* deriv = 0, double* deriv2 = 0)  const {
			if (value)*value = math::integrateGL(math::ScaledIntegrand<math::ScalingCub>
				(math::ScalingCub(0, z), pzfunc), 0, 1, math::MAX_GL_ORDER) - I;
			if (deriv) {
				double pz;
				pzfunc.evalDeriv(z, &pz);
				*deriv = pz;
			}
		}
		virtual unsigned int numDerivs() const { return 1; }
};
		//Bubble sort (x,y) in ascending order of x.
void sort(std::vector<double>& x, std::vector<double>& y) {
	int n = x.size();
	for (int i = 0; i < n - 1; i++) {
		bool swap = false;
		for (int j = 0; j < n - i - 1; j++) {
			if (x[j] > x[j + 1]) {
				std::swap(x[j], x[j + 1]);
				std::swap(y[j], y[j + 1]);
				swap = true;
			}
		}
		if (!swap) break;
	}
}
// computes area inside the curve defined by  the points (x,y) - used to compute Jr at box/loop
// transition. On return biggest values of x and y in x2max & y2max
// and in ysh the y-value of the curve when x=Rsh
double Area(std::vector<double> x, std::vector<double> y, const double Rsh,
	    double& x2max, double& ymv2, double& ysh) {
	std::vector<double> xn, yn;
	for (int i = 0; i < x.size(); i++) {//Assume 4-fold symmetry
		xn.push_back(fabs(x[i])); yn.push_back(fabs(y[i]));
	}
	sort(xn, yn);
	int n = yn.size();
	double I = yn[0] * (xn[1] - xn[0])
		   + yn[n - 1] * (xn[n - 1] - xn[n - 2]);
	for (int i = 1; i < n - 1; i++) 
		I +=  yn[i] * (xn[i + 1] - xn[i - 1]);
	I *= .5;
	for(int i=1; i<n; i++){
		if((Rsh-x[i-1]) * (x[i]-Rsh) < 0) continue;
		double f = (Rsh-x[i-1]) / (x[i]-x[i-1]);
		ysh = (1-f)*y[i-1] + f*y[i];
	}
	x2max = xn[xn.size()-1];
	ymv2 = yn[yn.size()-1];
	return I;
}

// computes Actions J of the orbit with enegy E at the box/loop
// transition. On return vR has Rdot(Rsh)
actions::Actions BoxLoopTrAct(const potential::BasePotential& pot, double E,
			      const double Rsh, double& vR) {
	pzf pzt(pot,E);
	double zmax = potential::z_max(pot, E);
			//Angle to z axis at origin
	double t = 0;
	double Rmax0 = potential::R_max(pot,E);
	double x0 = 1e-3*Rmax0;
	double Phi1 = pot.value(coord::PosCyl(x0, 0, 0));
	double v = sqrt(2 * (E - Phi1));
	double vz = v * cos(t), vx = v * sin(t);
	coord::PosVelCyl xv0(x0, 0, 0, vx, vz, 0);
	std::vector<double> R, pR;
	orbit::makeSoS(xv0, pot, R, pR, 1000);
	double Rmax, pRm;
	double Jr = Area(R, pR, Rsh, Rmax, pRm, vR) / M_PI;
	double v1 = sqrt(2*(E-pot.value(coord::PosCyl(Rmax,0,0))));
	if(pRm/v1>0.8){//From Rmax to Rmax0 approx motion as on x axis 
		pxf pxfunc(pot,E);
		Jr+=math::integrateGL(math::ScaledIntegrand<math::ScalingCub>
				      (math::ScalingCub(Rmax, Rmax0), pxfunc), 0, 1, math::MAX_GL_ORDER)/M_PI;
	}
	//Get Jfast=Jx+Jz=2*Jr+Jz from motion along the z axis
	double Jfast = 2 * math::integrateGL(math::ScaledIntegrand<math::ScalingCub>
					     (math::ScalingCub(0, zmax), pzt), 0, 1, math::MAX_GL_ORDER) / M_PI;
	double Jz = Jfast - 2 * Jr;
	if(Jz<0){
		Jz=Jfast;
		Jr=0;
	}
	return actions::Actions(Jr, Jz, 0);
}

/// function to use in ODE integrator
class OrbitIntegratorMeridionalPlane: public math::IOdeSystem {
	public:
		OrbitIntegratorMeridionalPlane(const potential::BasePotential& p, double Lz) :
		    poten(p), Lz2(Lz*Lz) {};

    /** apply the equations of motion in R,z plane without tracking the azimuthal motion.
        Integration variables are: R, z, vR, vz, dR, dz, dvR, dvz
        R here can have a negative sign (this happens for Lz=0, when the orbit flips to x<0
        and crosses the z=0 plane at negative 'R', but we compute the potential derivatives at |R|,
        and multiply by sign(R) when necessary.
    */
		virtual void eval(const double /*t*/, const double x[], double dxdt[]) const
		{
			coord::GradCyl grad;
			coord::HessCyl hess;
			double signR = x[0]>=0 ? 1 : -1;
			coord::PosCyl pos(fabs(x[0]), x[1], 0);
			poten.eval(pos, NULL, &grad, &hess);
			double Lz2ovR4 = Lz2>0 ? Lz2/pow_2(pow_2(pos.R)) : 0;
			dxdt[0] = x[2];
			dxdt[1] = x[3];
			dxdt[2] = -(grad.dR - Lz2ovR4 * pos.R) * signR;
			dxdt[3] = - grad.dz;
			dxdt[4] = x[6];
			dxdt[5] = x[7];
			dxdt[6] = -(hess.dR2 + 3*Lz2ovR4) * x[4] - hess.dRdz * signR * x[5];
			dxdt[7] = - hess.dRdz * signR * x[4] - hess.dz2 * x[5];
			dxdt[8] = pow_2(x[3]) + pow_2(x[2]);//dJz = vz*dvz + vR*dvR 
		}

		virtual unsigned int size() const { return 9; }  // two coordinates and two velocities
	private:
		const potential::BasePotential& poten;
		const double Lz2;
};

/// function to use in locating the exact time of the x-y plane crossing
class FindCrossingPointZequal0: public math::IFunction {
	public:
		FindCrossingPointZequal0(const math::BaseOdeSolver& _solver) :
		    solver(_solver) {};
    /** used in root-finder to locate the root z(t)=0 */
		virtual void evalDeriv(const double time, double* val, double* der, double*) const
		{
			if(val)
				*val = solver.getSol(time, 1);  // z
			if(der)
				*der = solver.getSol(time, 3);  // vz
		}
		virtual unsigned int numDerivs() const { return 1; }
	private:
		const math::BaseOdeSolver& solver;
};

/// function to use in locating the exact time vz goes -ve
class FindCrossingPointVZequal0: public math::IFunction {
	public:
		FindCrossingPointVZequal0(const math::BaseOdeSolver& _solver) :
		    solver(_solver) {};
    /** used in root-finder to locate the root z(t)=0 */
		virtual void evalDeriv(const double time, double* val, double* der, double*) const
		{
			if(val)	*val = solver.getSol(time, 3);  // Vz
		}
		virtual unsigned int numDerivs() const { return 0; }
	private:
		const math::BaseOdeSolver& solver;
};

/** launch an orbit perpendicularly to x-y plane from radius R0 with vz>0,
    and record the radius at which it crosses this plane downward (vz<0).
    \param[in]  poten  is the potential;
    \param[in]  E  is the orbit energy;
    \param[in]  Lz  is the z-component of angular momentum;
    \param[in]  R0  is the radius of the starting point;
    \param[out] timeCross stores the time required to complete the half-oscillation in z;
    \param[out] traj stores the trajectory recorded at equal intervals of time;
    \param[out] Rcross stores the radius of the crossing point;
    \param[out] dRcrossdR0  stores the derivative dRcross/dR0, computed from the variational equation.
*/
void findCrossingPointR(
			const potential::BasePotential& poten, double E, double Lz, double R0,
			double& timeCross,
			std::vector<std::pair<coord::PosVelCyl, double> >& traj, double& Rcross,
			double& dRcrossdR0, double& Jz)
{
	double Phi;
	coord::GradCyl grad;
	poten.eval(coord::PosCyl(R0, 0, 0), &Phi, &grad);
    // initial vertical velocity
	double vz0 = sqrt(fmax( 2 * (E-Phi) - (Lz>0 ? pow_2(Lz/R0) : 0), 0));
    // initial R-component of the deviation vector
	double dR0 = 1.;
    // initial vz-component (assigned from the requirement that E=const)
	double dvz0= vz0>0 ? ((Lz>0 ? pow_2(Lz) / pow_3(R0) : 0) - grad.dR) / vz0 * dR0 : 0;
	double vars[9] = {R0, 0, 0, vz0, dR0, 0, 0, dvz0, 0};
	OrbitIntegratorMeridionalPlane odeSystem(poten, Lz);
	math::OdeSolverDOP853 solver(odeSystem, ACCURACY_INTEGR);
	solver.init(vars);
	bool finished = false;
	unsigned int numStepsODE = 0;
	double timeCurr = 0;
	double timeTraj = 0;
	//Store the rising quarter of the orbit
	const double timeStepTraj = timeCross*0.5/(NUM_STEPS_TRAJ-1);
	traj.clear();
	while(!finished) {
		if(solver.doStep() <= 0 || numStepsODE >= MAX_NUM_STEPS_ODE) { // signal of error
			utils::msg(utils::VL_WARNING, FUNCNAME,
				   "Failed to compute shell orbit for E="+utils::toString(E,16)+
				   ", Lz="+utils::toString(Lz,16)+", R="+utils::toString(R0,16));
			timeCross  = 0;
			Rcross     = R0;   // this would terminate the root-finder, but we have no better option..
			dRcrossdR0 = NAN;
			printf("Failed to compute shell orbit for E=%g, Lz=%g, R=%g\n",E,Lz,R0);
			return;
		} else {
			numStepsODE++;
			double timePrev = timeCurr;
			timeCurr = solver.getTime();
			if(timeStepTraj!=INFINITY)
			{   // store first part of trajectory
				while(timeTraj <= timeCurr && traj.size() < NUM_STEPS_TRAJ) {
		    // store R, z, vR, vz at equal intervals of time
					double R = solver.getSol(timeTraj, 0);
					double z = solver.getSol(timeTraj, 1);
					double vR = solver.getSol(timeTraj, 2);
					double vz = solver.getSol(timeTraj, 3);
					traj.push_back(std::make_pair(coord::PosVelCyl(fabs(R), z, 0, vR, vz, 0),timeTraj));
					timeTraj += timeStepTraj;
				}
			}
			if(solver.getSol(timeCurr, 1) <= 0) {  // z<=0 - we're done
				finished = true;
				timeCurr = math::findRoot(FindCrossingPointZequal0(solver),
					timePrev, timeCurr, ACCURACY_RSHELL);
			}
		}
	}
	timeCross = timeCurr;    // the moment of crossing of the equatorial plane
	Rcross    = solver.getSol(timeCurr, 0);
	double vR = solver.getSol(timeCurr, 2);
	double vz = solver.getSol(timeCurr, 3);
	double dR = solver.getSol(timeCurr, 4);  // component of the deviation vector dR at the crossing
	double dz = solver.getSol(timeCurr, 5);  // -"- dz
	//Action obtained from entire solution over half period
	Jz = solver.getSol(timeCurr, 8)/M_PI;
	dRcrossdR0= dR - dz * vR / vz;
	if(Rcross < 0) {  // this happens for Lz=0, when the orbit crosses the x axis at negative x
		Rcross     = -Rcross;
		dRcrossdR0 = -dRcrossdR0;
	}
}

}//namespace anon

/* We integrate shell orbits on grid in E,Xi=Jphi/Jc(E) (energy and inclination)
 * For each of D,Rsh we produce 2 types of LinearInterpolator2d:
 * For actionFinder x,y are scaledE and scaledXi
 * For torusGenerator x,y are L,Xi
 * Note: Xi has different meaning in the 2 cases:
    for actionFinder Xi=Jphi/Lc(E)
    for TorusGenerator Xi=Jphi/(Jz+Jphi)
*/
EXP ShellInterpolator::ShellInterpolator(const BasePotential& pot){
	const double invPhi0 = 1/pot.value(coord::PosCyl(0,0,0));
	std::vector<double> gridR = potential::createInterpolationGrid
				    (pot, ACCURACY_INTERP2);
	const int sizeE=gridR.size(), sizeXi = 25;//sizeE/2;
	std::vector<double> gridE(sizeE), gridEscaled(sizeE);
	std::vector<double> gridL(sizeE), gridLscaled(sizeE);
	const math::ScalingSemiInf scalingSemi;
	for (int i = 0; i < sizeE; i++) {
		gridE[i] = pot.value(coord::PosCyl(gridR[i],0,0));
		gridEscaled[i] = scaleE(gridE[i], invPhi0);
		gridL[i] = L_circ(pot, gridE[i]);
		gridLscaled[i] = scale(scalingSemi,gridL[i]);
	}
	std::vector<double> gridXi(sizeXi), gridXiscaled(sizeXi);
	math::ScalingCub scalingXi(0, 1);//fix Xi grid so it's densest near ends
	for(int i=0; i<sizeXi; i++){
		gridXiscaled[i] = math::unscale(scalingXi, i/(double)(sizeXi-1));
		gridXi[i] = math::unscale(scalingXi, gridXiscaled[i]);  // Lzrel = u(chi)
	}
	math::Matrix<double> grid2dDE(sizeE,sizeXi);  // focal distance
	math::Matrix<double> grid2dRE(sizeE,sizeXi);  // Rshell / Rcirc(E)
	math::Matrix<double> grid2dL(sizeE,sizeXi);  // L = Jz+Jphi values
	math::Matrix<double> grid2dDL(sizeE,sizeXi);  // focal distance
	math::Matrix<double> grid2dRL(sizeE,sizeXi);  // Rshell / Rcirc(E)
	std::string errorMessage;  // store the error text in case of an exception in the openmp block
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
	for(int iE = 0; iE < sizeE-1; iE++) {//avoid almost free orbits
		double E  = gridE[iE];
		double Rc = R_circ(pot, E);
		double Jz, Lc = gridL[iE];
		std::vector<double> L_vals(sizeXi);
		std::vector<double> Xi_vals(sizeXi);
		std::vector<double> D_vals(sizeXi);
		std::vector<double> R_vals(sizeXi);
		for(int iXi=0; iXi < sizeXi-1; iXi++){//omit circ orbits
			try{
				double Jphi = Lc * gridXi[iXi];
				std::vector<coord::PosVelCyl> shell;
				double Rsh, FD = estimateFocalDistanceShellOrbit
					(pot, E, Jphi, &Rsh, &Jz, &shell);
				grid2dDE(iE, iXi) = FD;
				grid2dRE(iE, iXi) = Rsh / Rc;
				double L = Jphi + Jz;
				L_vals[iXi] = L; Xi_vals[iXi] = Jphi / L;
				D_vals[iXi] = FD; R_vals[iXi] = Rsh;
			}
			catch(std::exception& ex) {
				std::cout << ex.what() << "\n";
			}
		}
	// values for circular (planar) orbits
		grid2dDE(iE, 0) = grid2dDE(iE, 1);//don't trust Jphi=0 result
		grid2dDE(iE,sizeXi-1) = grid2dDE(iE, sizeXi-2);
		grid2dRE(iE,sizeXi-1) = 1;// Rsh = Rc
		D_vals[sizeXi - 1] = D_vals[sizeXi - 2];
		R_vals[sizeXi - 1] = R_vals[sizeXi - 2];
		L_vals[sizeXi - 1] = Lc;
		Xi_vals[sizeXi - 1] = 1;
	//Now interpolate values L, D, R at given E onto
	//regular grid in Jz/L 
		math::LinearInterpolator interpL(Xi_vals, L_vals);
		math::LinearInterpolator interpD(Xi_vals, D_vals);
		math::LinearInterpolator interpR(Xi_vals, R_vals);
		for (int iXi = 0; iXi < sizeXi; iXi++) {
			interpL.evalDeriv(gridXi[iXi], &grid2dL(iE, iXi));
			interpD.evalDeriv(gridXi[iXi], &grid2dDL(iE, iXi));
			interpR.evalDeriv(gridXi[iXi], &grid2dRL(iE, iXi));
		}
	}
    // limiting case of E=0 - assume a Keplerian potential at large radii
	for(int iXi=0; iXi<sizeXi; iXi++) {
		grid2dDE(sizeE-1, iXi) = grid2dDE(sizeE-2, iXi);
		grid2dRE(sizeE-1, iXi) = 1;  // Rshell = Rcirc
		grid2dL (sizeE-1, iXi) = L_circ(pot, gridE[sizeE-1]);
		grid2dDL(sizeE-1, iXi) = grid2dDL(sizeE-2,iXi);
		grid2dRL(sizeE-1, iXi) = R_circ(pot, gridE[sizeE-1]);// Rshell = Rcirc
	}
//	for(int i=0; i<sizeE; i+=2) printf("%g %g %g %g\n",
//					   gridE[i],grid2dDE(i,0),grid2dRE(i,0),gridR[i]);
	//grid2dD contains D on regular grid in Xi but irregular
	//values of L that are stored in grid2dL
	for (int iXi = 0; iXi < sizeXi; iXi++) {
		std::vector<double> L_vals(sizeE);
		std::vector<double> D_vals(sizeE);
		std::vector<double> R_vals(sizeE);
		for (int iE = 0; iE < sizeE; iE++) {
			L_vals[iE] = grid2dL(iE, iXi);
			D_vals[iE] = grid2dDL(iE, iXi);
			R_vals[iE] = grid2dRL(iE, iXi);
		}
		//We now have D and R at series of L values
		math::LinearInterpolator DL(L_vals, D_vals);
		math::LinearInterpolator RL(L_vals, R_vals);
		for (int iE = 0; iE < sizeE; iE++) {
			DL.evalDeriv(gridL[iE], &grid2dDL(iE, iXi));
			RL.evalDeriv(gridL[iE], &grid2dRL(iE, iXi));
		}

	}
	if(!errorMessage.empty())
		throw std::runtime_error(errorMessage);
	std::vector<double> gridJr, gridJz, gridI3, gridFD, scaledJ;
	interpDE = math::LinearInterpolator2d(gridEscaled, gridXiscaled, grid2dDE);
	interpRE = math::LinearInterpolator2d(gridEscaled, gridXiscaled, grid2dRE);
	interpDL = math::LinearInterpolator2d(gridLscaled, gridXiscaled, grid2dDL);
	interpRL = math::LinearInterpolator2d(gridLscaled, gridXiscaled, grid2dRL);
}
/* The main job of PolarInterpolator is to hold the curve Jz(Jf) along
 * which the box/loop transition lies. In adition it holds the values
 * of Delta(E) that cause the I3 centrifugal barrier to vanish on this
 * curve and te associaed I3(E), where I3 is computed from the velocity
 * of the transition orbit at (Rsh,0) 
*/
EXP PolarInterpolator::PolarInterpolator(const potential::BasePotential& pot,
					 const PtrShellInterpolator PtrShellI) {
			//spherical potential only has loop orbits by conservation of angular momentum.
			// Also if potential infinite as centre also always loop orbits.
	double Phi0 = pot.value(coord::PosCar(0, 0, 0));//potential at centre
	std::vector<double> gridR = potential::createInterpolationGrid
				    (pot, ACCURACY_INTERP2);
	const int sizeE=gridR.size();
	std::vector<double> gridE(sizeE), gridEscaled(sizeE);
	for (int i = 0; i < sizeE; i++) {
		gridE[i] = pot.value(coord::PosCyl(gridR[i],0,0));
		gridEscaled[i] = scaleE(gridE[i], 1/Phi0);
	}
	math::ScalingSemiInf Sc;
	std::vector<double> gridJr(sizeE), gridJz(sizeE), gridJfScaled(sizeE),
	gridI3, gridFD;
	if (potential::isSpherical(pot)|| std::isnan(Phi0) || std::isinf(Phi0)) {
		double Jr0 = 0;
		for (int i = 0; i < sizeE; i++) {
			gridJr[i] = Jr0;
			gridJz[i] = 0;
			gridFD.push_back(0);
			gridI3.push_back(0);//change this
			gridJfScaled[i]=scale(Sc,2*Jr0);
			Jr0 += 1.;
		}
	} else {
		int N = 5;//number of points to be fitted
		bool fitted = false;//gets if fitted straight point in E,z1 plane
		bool interp = false;
		double E0, z0;
		double b;//z=b(E-E0)+z0
		int i = 0;
		double Rsh, vR, umin, d2pu2du2;
		while(i<sizeE) {
			Rsh = PtrShellI->getRsh(gridE[i], 0.5, 1/Phi0) * R_circ(pot, gridE[i]);
			;
			if (gridE[i] / Phi0 > 0.2 && !interp) {
				actions::Actions Jcrit = BoxLoopTrAct(pot, gridE[i], Rsh, vR);
				double Delta0 = PtrShellI->getDelta(gridE[i], 0, 1/Phi0);
				FDfinder FDf(gridE[i], Rsh, vR, Delta0, pot);
				gridFD.push_back(FDf.bestFD(umin, d2pu2du2));
				const coord::ProlSph coordsys(gridFD[i]);
				coord::UVSph cs(gridFD[i]);
				double vz2 = 2*(gridE[i]-pot.value(coord::PosCyl(Rsh,0,0)))-pow_2(vR); 
				double vz = vz2>0? sqrt(vz2) : 0;
				coord::PosVelCyl point(Rsh,0,0,vR,vz,0);
				coord::PosMomUVSph uvp(coord::toPosMom(coord::PosMomCyl(Rsh,0,0,vR,vz,0),cs));
				gridI3.push_back(gridE[i]*pow_2(Rsh/cs.Delta) - .5*pow_2(uvp.pu)/cs.Delta2);
				gridJr[i] = Jcrit.Jr; gridJz[i] = Jcrit.Jz;
				if ((i > N+1) && (1-gridE[i]/Phi0) > 1e-4) {
					if (gridJz[i] < gridJz[i-1]) {
						interp = true;
					}
				}
			}
			else {
				if (!fitted) {
					interp = true;
					int index = i - 1;
					if (N > index)N = index;
					E0 = gridE[index];
					Ifind fI(.5 * M_PI * gridJz[index], pot, gridE[index]);
					double zmax0 = potential::z_max(pot, gridE[index]);
					z0 = math::findRoot(fI, 0, zmax0, 1e-6);
					std::vector<double> x(N - 1), y(N - 1);
					for (int j = 0;j < N - 1;j++) {
						int j1 = index + j - N + 1;
						Ifind fI(.5 * M_PI * gridJz[j1], pot, gridE[j1]);
						double zmax = potential::z_max(pot, gridE[j1]);
						x[j] = (gridE[j1] - E0);
						y[j] = (math::findRoot(fI, 0, zmax, 1e-6) - z0);
					}
					b = math::linearFitZero(x, y, NULL);
					fitted = true;
				}
				pzf pzfunc(pot,gridE[i]);
				double zmax = potential::z_max(pot,gridE[i]);
				double z1 = b * (gridE[i] - E0) + z0;
				gridJz[i] = 2 * math::integrateGL(pzfunc, 0, z1, math::MAX_GL_ORDER) / M_PI;
				gridJr[i] = math::integrateGL(pzfunc, z1, zmax, math::MAX_GL_ORDER) / M_PI;
			}
			i++;
		}
		int sizeFD = gridFD.size();
		for(int i=0; i<sizeE; i++){
			gridJfScaled[i]=scale(Sc,2*gridJr[i]+gridJz[i]);
			if(i>=sizeFD){
				gridFD.push_back(0);
				gridI3.push_back(0);
			}
		}
	}
	interpI3 = math::LinearInterpolator(gridEscaled, gridI3);
	interpFD = math::LinearInterpolator(gridEscaled, gridFD);
	coeffsJz = math::fitPoly(15, gridJfScaled, gridJz);
}

void findCrossingPointV(
			const potential::BasePotential& poten, double R0, double Jphi, double V0,
			double& timeCross, std::vector<std::pair<coord::PosVelCyl,double> >& traj,
			double& Rcross, double& dRcrossdV0, double& Jz)
{
	double Phi;
	coord::GradCyl grad;
	poten.eval(coord::PosCyl(R0, 0, 0), &Phi, &grad);
    // initial vertical velocity
	double vz0 = V0;
    // initial V-component of the deviation vector
	double dV0 = 1.;
	double vars[9] = {R0, 0, 0, vz0, 0, 0, 0, dV0, 0};
	OrbitIntegratorMeridionalPlane odeSystem(poten, Jphi);
	math::OdeSolverDOP853 solver(odeSystem, ACCURACY_INTEGR);
	solver.init(vars);
	bool finished = false;
	unsigned int numStepsODE = 0;
	double timeCurr = 0;
	double timeTraj = 0;
	//We only remember the rising quarter of the orbit 
	const double timeStepTraj = timeCross*0.5/(NUM_STEPS_TRAJ-1);
	traj.clear();
	while(!finished) {
		if(solver.doStep() <= 0 || numStepsODE >= MAX_NUM_STEPS_ODE) { // signal of error
			std::cout << "findCrossingPointV: Failed to compute orbit for Jphi="
					+utils::toString(Jphi,16)+", R="+utils::toString(R0,16)+
					": doStep() numStepsODE:"<<solver.doStep()<<" "<<numStepsODE<<"\n";
			timeCross  = 0;
			Rcross     = R0;   // this would terminate the root-finder, but we have no better option..
			dRcrossdV0 = NAN;
			return;
		} else {
			numStepsODE++;
			double timePrev = timeCurr;
			timeCurr = solver.getTime();
			if(timeStepTraj!=INFINITY)
			{   // store first part of trajectory
				while(timeTraj <= timeCurr && traj.size() < NUM_STEPS_TRAJ) {
		    // store R and z at equal intervals of time
					double R = solver.getSol(timeTraj, 0);
					double z = solver.getSol(timeTraj, 1);
					double vR = solver.getSol(timeTraj, 2);
					double vz = solver.getSol(timeTraj, 3);
					traj.push_back(std::make_pair(coord::PosVelCyl(R, z, 0, vR, vz, 0),timeTraj));
					timeTraj += timeStepTraj;
				}
			}
			if(solver.getSol(timeCurr, 1) <= 0) {  // z<=0 - we're done
				finished = true;
				timeCurr = math::findRoot(FindCrossingPointZequal0(solver),
					timePrev, timeCurr, ACCURACY_RSHELL);
			}
		}
	}
	timeCross = timeCurr;    // the moment of crossing of the equatorial plane
	Rcross    = solver.getSol(timeCurr, 0);
	double vR = solver.getSol(timeCurr, 2);
	double vz = solver.getSol(timeCurr, 3);
	double dR = solver.getSol(timeCurr, 4);  // component of the deviation vector dR at the crossing
	double dz = solver.getSol(timeCurr, 5);  // -"- dz
	//Action obtained from entire solution over half period
	Jz = solver.getSol(timeCurr, 8)/M_PI;
	dRcrossdV0= dR - dz * vR / vz;
	if(Rcross < 0) {  // this happens for Jphi=0, when the orbit crosses the x axis at negative x
		Rcross     = -Rcross;
		dRcrossdV0 = -dRcrossdV0;
	}
}


EXP void FindClosedOrbitRZplane::evalDeriv(const double R0, double* val, double* der, double*) const {
	// first two calls in root-finder are for the boundary points, we already know the answer
	if(R0==Rmin || R0==Rmax) {
		if(val) *val = R0==Rmin ? Rmax-Rmin : Rmin-Rmax;
		if(der) *der = NAN;
		return;
	}
	double Rcross, dRcrossdR=NAN;
	findCrossingPointR(poten, E, Jphi, R0, timeCross, traj, Rcross, dRcrossdR, Jz);
	if(val)
		*val = Rcross-R0;
	if(der)
		*der = dRcrossdR-1;
}

EXP void FindRzClosedOrbitV::evalDeriv(const double V0, double* val, double* der, double*) const {
	double Rcross, dRcrossdV=NAN;
	findCrossingPointV(poten, R0, Jphi, V0, timeCross, traj, Rcross, dRcrossdV, Jz);
	if(val)
		*val = Rcross-R0;
	if(der)
		*der = dRcrossdV;
}

EXP double estimateFocalDistanceShellOrbit(
					   const potential::BasePotential& poten, double E, double Jphi,
					   double* _Rshell, double* _Jz, std::vector<coord::PosVelCyl>* shell)
{
	double Rmin, Rmax, FD;
	findPlanarOrbitExtent(poten, E, Jphi, Rmin, Rmax);
	double timeCross = 5*pow_2(Rmin)/Jphi, Jz;
	std::vector<std::pair<coord::PosVelCyl,double> > traj;
	FindClosedOrbitRZplane finder(poten, E, Jphi, Rmin, Rmax, timeCross, traj, Jz);
    // locate the radius of a shell orbit;  as a by-product, store the orbit in 'traj'
	double Rshell = math::findRoot(finder, Rmin, Rmax, ACCURACY_RSHELL);
	if(shell){//return shell orbit up to p_theta=0
		for(int i=0; i<traj.size(); i++)
			shell->push_back(traj[i].first);
	}

#ifdef OLD_METHOD
	if(traj.size() >= 2)
	// now find the best-fit value of delta for this orbit
		FD = fitFocalDistanceShellOrbit(traj);
	else {
	// something went wrong; use a backup solution
		if(!isFinite(Rshell))
			Rshell = 0.5 * (Rmin+Rmax);
		utils::msg(utils::VL_WARNING, FUNCNAME,
			   "Could not find a thin orbit for E="+utils::toString(E,16)+", Jphi="+utils::toString(Jphi,16)+
			   " - assuming Rthin="+utils::toString(Rshell,16));
	// if we don't have a proper orbit, make a short vertical step out of the z=0 plane
	// and estimate the focal distance from the mixed derivative at this single point
		FD = estimateFocalDistancePoints(poten, std::vector<coord::PosCyl>(1,
			coord::PosCyl(Rshell, /* z=very small number */ Rshell*ACCURACY_RSHELL, 0)));
	}
#else
	double Phi;
	coord::GradCyl grad;
	poten.eval(coord::PosCyl(Rshell,0,0), &Phi, &grad);
	double vphi = Jphi!=0 ? Jphi / Rshell : 0;
	FD = Rshell * sqrt( math::clip((2 * (E-Phi) - Rshell * grad.dR) / ( Rshell * grad.dR - vphi*vphi), 0., 1e6) );
#endif
	if(_Jz) *_Jz = Jz;
	if(_Rshell) *_Rshell = Rshell; 
	return FD;
}

}//namespace