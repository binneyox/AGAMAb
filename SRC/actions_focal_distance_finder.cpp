#include "actions_focal_distance_finder.h"
#include <iostream>
#include <cassert>
#include <stdexcept>
#include <cmath>

namespace actions{

/// number of sampling points for a shell orbit (equally spaced in time)
static const unsigned int NUM_STEPS_TRAJ = 64;
/// accuracy of root-finding for the radius of thin (shell) orbit
static const double ACCURACY_RSHELL = 1e-6;
/// accuracy of orbit integration for shell orbit
static const double ACCURACY_INTEGR = 1e-8;
/// upper limit on the number of timesteps in ODE solver (should be enough to track half of the orbit)
static const unsigned int MAX_NUM_STEPS_ODE = 2000;


// estimate focal distance for a series of points in R-z plane
template<typename PointT>
EXP double estimateFocalDistancePoints(
    const potential::BasePotential& potential, const std::vector<PointT>& traj)
{
    if(traj.size()==0)
        throw std::invalid_argument("Error in finding focal distance: empty array of points");
    std::vector<double> x(traj.size()), y(traj.size());
    double sumsq = 0;
    for(unsigned int i=0; i<traj.size(); i++) {
        const coord::PosCyl p = coord::toPosCyl(traj[i]);
        coord::GradCyl grad;
        coord::HessCyl hess;
        potential.eval(p, NULL, &grad, &hess);
        x[i] = hess.dRdz;
        y[i] = 3*p.z * grad.dR - 3*p.R * grad.dz + p.R*p.z * (hess.dR2-hess.dz2)
             + (p.z*p.z - p.R*p.R) * hess.dRdz;
        sumsq += pow_2(x[i]);
    }
    double result = sumsq>0 ? math::linearFitZero(x, y, NULL) : 0;
    return sqrt( fmax( result, 0) );  // ensure that the computed value is non-negative
}

template EXP double estimateFocalDistancePoints(
    const potential::BasePotential& potential, const std::vector<coord::PosCar>& traj);
template EXP double estimateFocalDistancePoints(
    const potential::BasePotential& potential, const std::vector<coord::PosVelCar>& traj);
template EXP double estimateFocalDistancePoints(
    const potential::BasePotential& potential, const std::vector<coord::PosCyl>& traj);
template EXP double estimateFocalDistancePoints(
    const potential::BasePotential& potential, const std::vector<coord::PosVelCyl>& traj);
template EXP double estimateFocalDistancePoints(
    const potential::BasePotential& potential, const std::vector<coord::PosSph>& traj);
template EXP double estimateFocalDistancePoints(
    const potential::BasePotential& potential, const std::vector<coord::PosVelSph>& traj);

/** find the best-fit value of focal distance for a shell orbit.
    \param[in] traj  contains the trajectory of this orbit in R-z plane,
    \return  the parameter `delta` of a prolate spheroidal coordinate system which minimizes
    the variation of `lambda` coordinate for this orbit
    If the best-fit value is negative, it is replaced with zero.
*/
EXP double fitFocalDistanceShellOrbit(const std::vector<coord::PosCyl>& traj)
{
	if(traj.size()==0)
		throw std::invalid_argument("Error in finding focal distance for a shell orbit: empty array");
	math::Matrix<double> coefs(traj.size(), 2);
	std::vector<double> rhs(traj.size());
	std::vector<double> result;  // regression parameters:  lambda/(lambda-delta), lambda
	for(unsigned int i=0; i<traj.size(); i++) {
		coefs(i, 0) = -pow_2(traj[i].R);
		coefs(i, 1) = 1.;
		rhs[i] = pow_2(traj[i].z);
	}
	math::linearMultiFit(coefs, rhs, NULL, result);
	return sqrt( fmax( result[1] * (1 - 1/result[0]), 0) );
}

namespace{

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

}  // namespace

EXP std::vector<coord::PosVelCyl> toTop(const potential::BasePotential& poten,
		double Rsh, double Vsh, double Jphi, double& riseTime, double& Jz)
{
	double dV0 = 1.;
	double vars[9] = {Rsh, 0, 0, Vsh, 0, 0, 0, dV0, 0};
	OrbitIntegratorMeridionalPlane odeSystem(poten, Jphi);
	math::OdeSolverDOP853 solver(odeSystem, ACCURACY_INTEGR);
	solver.init(vars);
	double timeCurr = 0;
	double timeTraj = 0;
	const double timeStepTraj = riseTime/(NUM_STEPS_TRAJ-1);
	std::vector<coord::PosVelCyl> shell;
	bool finished = false;
	while(!finished) {
		solver.doStep();
		double timePrev = timeCurr;
		timeCurr = solver.getTime();
		while(timeTraj <= timeCurr) {
		    // store R and z at equal intervals of time
			double R  = solver.getSol(timeTraj, 0);
			double z  = solver.getSol(timeTraj, 1);
			double vR = solver.getSol(timeTraj, 2);
			double vz = solver.getSol(timeTraj, 3);
			if(vz>0) shell.push_back(coord::PosVelCyl(R, z, 0, vR, vz, 0));
			timeTraj += timeStepTraj;
		}
		if(solver.getSol(timeCurr, 3) <= 0) {  // Vz<=0 - we're done
			finished = true;
			timeCurr = math::findRoot(FindCrossingPointVZequal0(solver),
				timePrev, timeCurr, ACCURACY_RSHELL);
		}
	}
	riseTime=timeCurr;
	double R  = solver.getSol(timeCurr, 0);
	double z  = solver.getSol(timeCurr, 1);
	double vR = solver.getSol(timeCurr, 2);
	double vz = solver.getSol(timeCurr, 3);
	shell.push_back(coord::PosVelCyl(R, z, 0, vR, vz, 0));
	//Action obtained from entire solution over quarter period
	Jz = 2*solver.getSol(timeCurr, 8)/M_PI;
	return shell;
}

/** launch an orbit perpendicularly to x-y plane from radius R0 with vz>0,
    and record the radius at which it crosses this plane downward (vz<0).
    We vary vz rather than R and don't hold E=const.
    \param[in]  poten  is the potential;
    \param[in]  Lz  is the z-component of angular momentum;
    \param[in]  R0  is the radius of the starting point;
    \param[out] timeCross stores the time required to complete the half-oscillation in z;
    \param[out] traj stores the trajectory recorded at equal intervals of time;
    \param[out] Rcross stores the radius of the crossing point;
    \param[out] dRcrossdV0  stores the derivative dRcross/dV0, computed from the variational equation.
*/
void findCrossingPointV(
			const potential::BasePotential& poten, double R0, double Lz, double V0,
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
	OrbitIntegratorMeridionalPlane odeSystem(poten, Lz);
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
			std::cout << "findCrossingPointV: Failed to compute orbit for Lz="
					+utils::toString(Lz,16)+", R="+utils::toString(R0,16)+
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
	if(Rcross < 0) {  // this happens for Lz=0, when the orbit crosses the x axis at negative x
		Rcross     = -Rcross;
		dRcrossdV0 = -dRcrossdV0;
	}
}


void FindClosedOrbitRZplane::evalDeriv(const double R0, double* val, double* der, double*) const {
	// first two calls in root-finder are for the boundary points, we already know the answer
	if(R0==Rmin || R0==Rmax) {
		if(val) *val = R0==Rmin ? Rmax-Rmin : Rmin-Rmax;
		if(der) *der = NAN;
		return;
	}
	double Rcross, dRcrossdR=NAN;
	findCrossingPointR(poten, E, Lz, R0, timeCross, traj, Rcross, dRcrossdR, Jz);
	if(val)
		*val = Rcross-R0;
	if(der)
		*der = dRcrossdR-1;
}

void FindRzClosedOrbitV::evalDeriv(const double V0, double* val, double* der, double*) const {
	double Rcross, dRcrossdV=NAN;
	findCrossingPointV(poten, R0, Lz, V0, timeCross, traj, Rcross, dRcrossdV, Jz);
	if(val)
		*val = Rcross-R0;
	if(der)
		*der = dRcrossdV;
}

EXP double estimateFocalDistanceShellOrbit(
					   const potential::BasePotential& poten, double E, double Lz,
					   double* _Rshell, double* _Jz, std::vector<coord::PosVelCyl>* shell)
{
    double Rmin, Rmax, FD;
    findPlanarOrbitExtent(poten, E, Lz, Rmin, Rmax);
    double timeCross = 5*pow_2(Rmin)/Lz, Jz;
    std::vector<std::pair<coord::PosVelCyl,double> > traj;
    FindClosedOrbitRZplane finder(poten, E, Lz, Rmin, Rmax, timeCross, traj, Jz);
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
            "Could not find a thin orbit for E="+utils::toString(E,16)+", Lz="+utils::toString(Lz,16)+
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
    double vphi = Lz!=0 ? Lz / Rshell : 0;
    FD = Rshell * sqrt( math::clip((2 * (E-Phi) - Rshell * grad.dR) / ( Rshell * grad.dR - vphi*vphi), 0., 1e6) );
#endif
    if(_Jz) *_Jz = Jz;
    if(_Rshell) *_Rshell = Rshell; 
    return FD;
}

}  // namespace actions
