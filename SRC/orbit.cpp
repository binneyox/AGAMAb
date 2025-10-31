#include "orbit.h"
#include "potential_base.h"
#include "utils.h"
#include "math_core.h"
#include <stdexcept>
#include <cmath>

namespace orbit{

namespace{
// normalize the position-velocity returned by orbit integrator in case of r<0 or R<0
template<typename CoordT>
inline EXP coord::PosVelT<CoordT> getPosVel(const double data[6]) { return coord::PosVelT<CoordT>(data); }

template<>
inline EXP coord::PosVelCyl getPosVel(const double data[6])
{
    if(data[0] >= 0)
        return coord::PosVelCyl(data);
    else
        return coord::PosVelCyl(-data[0], data[1], data[2]+M_PI, -data[3], data[4], -data[5]);
}

template<>
inline EXP coord::PosVelSph getPosVel(const double data[6])
{
    double r = data[0];
    double phi = data[2];
    int signr = 1, signt = 1;
    double theta = fmod(data[1], 2*M_PI);
    // normalize theta to the range 0..pi, and r to >=0
    if(theta<-M_PI) {
        theta += 2*M_PI;
    } else if(theta<0) {
        theta = -theta;
        signt = -1;
    } else if(theta>M_PI) {
        theta = 2*M_PI-theta;
        signt = -1;
    }
    if(r<0) {
        r = -r;
        theta = M_PI-theta;
        signr = -1;
    }
    if((signr == -1) ^ (signt == -1))
        phi += M_PI;
    phi = math::wrapAngle(phi);
    return coord::PosVelSph(r, theta, phi, data[3] * signr, data[4] * signt, data[5] * signr * signt);
}
/// function to use in locating the exact time of crossing the x-y plane
class FindCrossingPointZequal0: public math::IFunction {
	public:
		FindCrossingPointZequal0(const math::BaseOdeSolver& _solver, const double _z0=0) :
		    solver(_solver), z0(_z0) {};
    /** used in root-finder to locate the root z(t)=0 */
		virtual void evalDeriv(const double time, double* val, double* der, double*) const
		{
			if(val)
				*val = solver.getSol(time, 1) - z0;  // z
			if(der)
				*der = solver.getSol(time, 3);  // vz
		}
		virtual unsigned int numDerivs() const { return 1; }
	private:
		const math::BaseOdeSolver& solver;
		const double z0;
};
/// function to use in locating the exact time of crossing r=rbar
class FindCrossingPoint_rbar: public math::IFunction {
	public:
		FindCrossingPoint_rbar(const math::BaseOdeSolver& _solver, const double _rbar) :
		    solver(_solver), rbar(_rbar) {};
    /** used in root-finder to locate the root r-rbar=0 */
		virtual void evalDeriv(const double time, double* val, double* der, double*) const
		{
			double R=solver.getSol(time, 0), z=solver.getSol(time, 1);
			double r = sqrt(R*R + z*z);
			if(val)
				*val = r - rbar;
			if(der)
				*der = (R*solver.getSol(time, 2) + z*solver.getSol(time, 3))/r;
		}
		virtual unsigned int numDerivs() const { return 1; }
	private:
		const math::BaseOdeSolver& solver;
		const double rbar;
};

/// function to use in ODE integrator
class EXP OrbitIntegratorRzPlane: public math::IOdeSystem {
	public:
		OrbitIntegratorRzPlane(const potential::BasePotential& p, double Lz) :
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
			double signR = x[0]>=0 ? 1 : -1;
			coord::PosCyl pos(fabs(x[0]), x[1], 0);
			poten.eval(pos, NULL, &grad, NULL);
			double Lz2ovR4 = Lz2>0 ? Lz2/pow_2(pow_2(pos.R)) : 0;
			dxdt[0] = x[2];
			dxdt[1] = x[3];
			dxdt[2] = -(grad.dR - Lz2ovR4 * pos.R) * signR;
			dxdt[3] = - grad.dz;
		}

		virtual unsigned int size() const { return 4; }  // two coordinates and two velocities
	private:
		const potential::BasePotential& poten;
		const double Lz2;
};
/// function to use in ODE integrator for HenonHeiles etc
class EXP OrbitIntegratorXzPlane: public math::IOdeSystem {
	public:
		OrbitIntegratorXzPlane(const potential::BasePotential& p) :
		    poten(p) {};

    /** apply the equations of motion in x,z plane when Lz=0.
        Integration variables are: x, z, vx, vz, dR, dz, dvR, dvz
    */
		virtual void eval(const double /*t*/, const double x[], double dxdt[]) const
		{
			coord::GradCyl grad;
			coord::PosCyl pos(x[0], x[1], 0);
			poten.eval(pos, NULL, &grad, NULL);
			dxdt[0] = x[2];
			dxdt[1] = x[3];
			dxdt[2] = - grad.dR;
			dxdt[3] = - grad.dz;
		}

		virtual unsigned int size() const { return 4; }  // two coordinates and two velocities
	private:
		const potential::BasePotential& poten;
};

}

static const double ACCURACY_INTEGR = 1e-8;

static const double ACCURACY_Zcoord = 1e-5;
/// upper limit on the number of timesteps in ODE solver
static const unsigned int MAX_NUM_STEPS_ODE = 20000;

template<typename CoordT>
StepResult RuntimeTrajectory<CoordT>::processTimestep(
    const math::BaseOdeSolver& solver, const double tbegin, const double tend, double vars[])
{
    if(samplingInterval == INFINITY) {
        // store just one point from the trajectory, which always contains the current orbital state and time
        if(trajectory.empty())
            trajectory.resize(1);
        trajectory[0].first  = getPosVel<CoordT>(vars);
        trajectory[0].second = tend;
    } else if(samplingInterval > 0) {
        // store trajectory at regular intervals of time
        size_t ibegin = static_cast<size_t>(tbegin / samplingInterval);
        size_t iend   = static_cast<size_t>(tend   / samplingInterval + 1e-10);
        trajectory.resize(iend + 1);
        for(size_t iout=ibegin; iout<=iend; iout++) {
            double tout = samplingInterval * iout;
            if(tout >= tbegin && tout <= tend) {
                double data[6];
                for(int d=0; d<6; d++)
                    data[d] = solver.getSol(tout, d);
                trajectory[iout].first  = getPosVel<CoordT>(data);
                trajectory[iout].second = tout;
            }
        }
    } else {
        // store trajectory at every integration timestep
        if(trajectory.empty()) {
            // add the initial point
            double data[6];
            for(int d=0; d<6; d++)
                data[d] = solver.getSol(tbegin, d);
            trajectory.push_back(
                std::pair<coord::PosVelT<CoordT>, double>(getPosVel<CoordT>(data), tbegin));
        }
        // add the current point (at the end of the timestep)
        trajectory.push_back(
            std::pair<coord::PosVelT<CoordT>, double>(getPosVel<CoordT>(vars), tend));
    }
    return SR_CONTINUE;
}


void OrbitIntegratorRot::eval(const double /*t*/, const double x[], double dxdt[]) const
{
    coord::GradCar grad;
    potential.eval(coord::PosCar(x[0], x[1], x[2]), NULL, &grad);
    // time derivative of position
    dxdt[0] = x[3] + Omega * x[1];
    dxdt[1] = x[4] - Omega * x[0];
    dxdt[2] = x[5];
    // time derivative of velocity
    dxdt[3] = -grad.dx + Omega * x[4];
    dxdt[4] = -grad.dy - Omega * x[3];
    dxdt[5] = -grad.dz;
}

double OrbitIntegratorRot::getAccuracyFactor(const double /*t*/, const double x[]) const
{
    double Epot = potential.value(coord::PosCar(x[0], x[1], x[2]));
    double Ekin = 0.5 * (x[3]*x[3] + x[4]*x[4] + x[5]*x[5]);
    return fmin(1, fabs(Epot + Ekin) / fmax(fabs(Epot), Ekin));
}

template<>
void OrbitIntegrator<coord::Car>::eval(const double /*t*/, const double x[], double dxdt[]) const
{
    coord::GradCar grad;
    potential.eval(coord::PosCar(x[0], x[1], x[2]), NULL, &grad);
    // time derivative of position
    dxdt[0] = x[3];
    dxdt[1] = x[4];
    dxdt[2] = x[5];
    // time derivative of velocity
    dxdt[3] = -grad.dx;
    dxdt[4] = -grad.dy;
    dxdt[5] = -grad.dz;
}

template<>
void OrbitIntegrator<coord::Cyl>::eval(const double /*t*/, const double x[], double dxdt[]) const
{
    coord::PosVelCyl p(x);
    if(x[0]<0) {    // R<0
        p.R = -p.R; // apply reflection
        p.phi += M_PI;
    }
    coord::GradCyl grad;
    potential.eval(p, NULL, &grad);
    double Rinv = p.R!=0 ? 1/p.R : 0;  // avoid NAN in degenerate cases
    dxdt[0] = p.vR;
    dxdt[1] = p.vz;
    dxdt[2] = p.vphi * Rinv;
    dxdt[3] = (x[0]<0 ? grad.dR : -grad.dR) + pow_2(p.vphi) * Rinv;
    dxdt[4] = -grad.dz;
    dxdt[5] = -(grad.dphi + p.vR*p.vphi) * Rinv;
}

template<>
void OrbitIntegrator<coord::Sph>::eval(const double /*t*/, const double x[], double dxdt[]) const
{
    double r = x[0];
    double phi = x[2];
    int signr = 1, signt = 1;
    double theta = fmod(x[1], 2*M_PI);
    // normalize theta to the range 0..pi, and r to >=0
    if(theta<-M_PI) {
        theta += 2*M_PI;
    } else if(theta<0) {
        theta = -theta;
        signt = -1;
    } else if(theta>M_PI) {
        theta = 2*M_PI-theta;
        signt = -1;
    }
    if(r<0) {
        r = -r;
        theta = M_PI-theta;
        signr = -1;
    }
    if((signr == -1) ^ (signt == -1))
        phi += M_PI;
    const coord::PosVelSph p(r, theta, phi, x[3], x[4], x[5]);
    coord::GradSph grad;
    potential.eval(p, NULL, &grad);
    double rinv = r!=0 ? 1/r : 0, sintheta, costheta;
    math::sincos(theta, sintheta, costheta);
    double sinthinv = sintheta!=0 ? 1./sintheta : 0;
    double cottheta = costheta * sinthinv;
    dxdt[0] = p.vr;
    dxdt[1] = p.vtheta * rinv;
    dxdt[2] = p.vphi * rinv * sinthinv;
    dxdt[3] = -grad.dr*signr + (pow_2(p.vtheta) + pow_2(p.vphi)) * rinv;
    dxdt[4] = (-grad.dtheta*signt + pow_2(p.vphi)*cottheta - p.vr*p.vtheta) * rinv;
    dxdt[5] = (-grad.dphi * sinthinv - (p.vr+p.vtheta*cottheta) * p.vphi) * rinv;
}

template<typename CoordT>
double OrbitIntegrator<CoordT>::getAccuracyFactor(const double /*t*/, const double x[]) const
{
    double Epot = potential.value(coord::PosT<CoordT>(x[0], x[1], x[2]));
    double Ekin = 0.5 * (x[3]*x[3] + x[4]*x[4] + x[5]*x[5]);
    return fmin(1, fabs(Epot + Ekin) / fmax(fabs(Epot), Ekin));
}


template<typename CoordT>
coord::PosVelT<CoordT> integrate(
    const coord::PosVelT<CoordT>& initialConditions,
    const double totalTime,
    const math::IOdeSystem& orbitIntegrator,
    const RuntimeFncArray& runtimeFncs,
    const OrbitIntParams& params)
{
    math::OdeSolverDOP853 solver(orbitIntegrator, params.accuracy);
    int NDIM = orbitIntegrator.size();
    if(NDIM < 6)
        throw std::runtime_error("orbit::integrate() needs at least 6 variables");
    std::vector<double> state(NDIM);
    double* vars = &state.front();
    initialConditions.unpack_to(vars);  // first 6 variables are always position/velocity
    solver.init(vars);
    size_t numSteps = 0;
    double timeCurr = 0;
    while(timeCurr < totalTime) {
        if(!(solver.doStep() > 0.)) {
            // signal of error
            utils::msg(utils::VL_WARNING, "orbit::integrate", "terminated at t="+utils::toString(timeCurr));
            break;
        }
        double timePrev = timeCurr;
        timeCurr = std::min(solver.getTime(), totalTime);
        for(int d=0; d<NDIM; d++)
            vars[d] = solver.getSol(timeCurr, d);
        bool reinit = false, finish = timeCurr >= totalTime;
        for(size_t i=0; i<runtimeFncs.size(); i++) {
            switch(runtimeFncs[i]->processTimestep(solver, timePrev, timeCurr, vars))
            {
                case orbit::SR_TERMINATE: finish = true; break;
                case orbit::SR_REINIT:    reinit = true; break;
                default: /*nothing*/;
            }
        }
        if(reinit)
            solver.init(vars);
        if(finish || ++numSteps >= params.maxNumSteps)
            break;
    }
    return coord::PosVelT<CoordT>(vars);
}

EXP coord::PosVelCyl makeSoS(const coord::PosVelCyl& Rz, const potential::BasePotential& poten,
			     std::vector<double>& Rs,
			     std::vector<double>& pRs, int Npt, const double z0)
{
	const double Lz=Rz.R*Rz.vphi;
	double vars[4] = {Rz.R, Rz.z, Rz.vR, Rz.vz};
	OrbitIntegratorRzPlane odeSystem(poten, Lz);
	math::OdeSolverDOP853 solver(odeSystem, ACCURACY_INTEGR);
	solver.init(vars);
	bool finished = false;
	unsigned int numStepsODE = 0;
	double timeCurr = 0;
	double z_last = Rz.z;
	Rs.clear(); pRs.clear();
	while(!finished) {
		if(solver.doStep() <= 0 || numStepsODE >= MAX_NUM_STEPS_ODE) { // signal of error
			return Rz;
		} else {
			numStepsODE++;
			double timePrev = timeCurr;
			timeCurr = solver.getTime();
			// Check for upward crossing of z=z0
			if((z_last-z0) * (solver.getSol(timeCurr, 1)-z0) <= 0 && solver.getSol(timeCurr,3)>0) {
				double timeCross = math::findRoot(FindCrossingPointZequal0(solver,z0),
					timePrev, timeCurr, ACCURACY_Zcoord);
				Rs.push_back( solver.getSol(timeCross, 0));
				pRs.push_back(solver.getSol(timeCross, 2));
			}
			z_last = solver.getSol(timeCurr, 1);
			finished = Rs.size()>=Npt;
		}
	}
	return coord::PosVelCyl(solver.getSol(timeCurr,0),solver.getSol(timeCurr,1),0,
				solver.getSol(timeCurr,2),solver.getSol(timeCurr,3),
				Lz/solver.getSol(timeCurr,0));
}

EXP coord::PosVelCyl makeSoS(const coord::PosVelCyl& Rz, const potential::BasePotential& poten,
			     std::vector<double>& Rs, std::vector<double>& pRs, const double rbar,
			     std::vector<double>& thetas, std::vector<double>& pthetas,
			     int Npt, const double z0)
{
	const double Lz=Rz.R*Rz.vphi;
	double vars[4] = {Rz.R, Rz.z, Rz.vR, Rz.vz};
	OrbitIntegratorRzPlane odeSystem(poten, Lz);
	math::OdeSolverDOP853 solver(odeSystem, ACCURACY_INTEGR);
	solver.init(vars);
	bool finished = false;
	unsigned int numStepsODE = 0;
	double timeCurr = 0;
	double z_last = Rz.z, dr_last = sqrt(Rz.R*Rz.R + Rz.z*Rz.z) - rbar;
	Rs.clear(); pRs.clear();
	while(!finished) {
		if(solver.doStep() <= 0 || numStepsODE >= MAX_NUM_STEPS_ODE) { // signal of error
			return Rz;
		} else {
			numStepsODE++;
			double timePrev = timeCurr;
			timeCurr = solver.getTime();
			// Check for upward crossing of z=0
			if((z_last-z0) * (solver.getSol(timeCurr, 1)-z0) <= 0 && solver.getSol(timeCurr,3)>0) {
				double timeCross = math::findRoot(FindCrossingPointZequal0(solver,z0),
					timePrev, timeCurr, ACCURACY_Zcoord);
				Rs.push_back( solver.getSol(timeCross, 0));
				pRs.push_back(solver.getSol(timeCross, 2));
			}
			// Check for outward crossing of r=rbar
			double r=sqrt(pow_2(solver.getSol(timeCurr,0)) + pow_2(solver.getSol(timeCurr,1)));
			if(dr_last * (r - rbar) <= 0 &&
			   solver.getSol(timeCurr,2)>0) {
				double timeCross = math::findRoot(FindCrossingPoint_rbar(solver, rbar),
					timePrev, timeCurr, ACCURACY_Zcoord);
				coord::PosMomCyl Rzp(solver.getSol(timeCross,0),solver.getSol(timeCross,1),0,
					solver.getSol(timeCross,2),solver.getSol(timeCross,3),Lz);
				coord::PosMomSph rtheta(coord::toPosMomSph(Rzp));
				thetas.push_back(rtheta.theta);
				pthetas.push_back(rtheta.ptheta);
			}
			z_last = solver.getSol(timeCurr, 1);
			dr_last = r - rbar;
			finished = Rs.size()>=Npt;
		}
	}
	return coord::PosVelCyl(solver.getSol(timeCurr,0),solver.getSol(timeCurr,1),0,
				solver.getSol(timeCurr,2),solver.getSol(timeCurr,3),
				Lz/solver.getSol(timeCurr,0));
}

EXP coord::PosVelCyl makeSoSXz(const coord::PosVelCyl& Rz, const potential::BasePotential& poten,
			     std::vector<double>& Rs,
			     std::vector<double>& pRs, int Npt)
{
	double vars[4] = {Rz.R, Rz.z, Rz.vR, Rz.vz};
	OrbitIntegratorXzPlane odeSystem(poten);
	math::OdeSolverDOP853 solver(odeSystem, ACCURACY_INTEGR);
	solver.init(vars);
	bool finished = false;
	unsigned int numStepsODE = 0;
	double timeCurr = 0;
	double z_last = Rz.z;
	Rs.clear(); pRs.clear();
	while(!finished) {
		if(solver.doStep() <= 0 || numStepsODE >= MAX_NUM_STEPS_ODE) { // signal of error
			printf("Error in orbit: %f %d\n",solver.doStep(),numStepsODE);
			return Rz;
		} else {
			numStepsODE++;
			double timePrev = timeCurr;
			timeCurr = solver.getTime();
			if(z_last * solver.getSol(timeCurr, 1) <= 0 && solver.getSol(timeCurr,3)>0) {
				double timeCross = math::findRoot(FindCrossingPointZequal0(solver),
					timePrev, timeCurr, ACCURACY_Zcoord);
				Rs.push_back( solver.getSol(timeCross, 0));
				pRs.push_back(solver.getSol(timeCross, 2));
			}
			z_last = solver.getSol(timeCurr, 1);
			finished = Rs.size()>=Npt;
		}
	}
	return coord::PosVelCyl(solver.getSol(timeCurr,0),solver.getSol(timeCurr,1),0,
				solver.getSol(timeCurr,2),solver.getSol(timeCurr,3),0);
}

// explicit template instantiations to make sure all of them get compiled
template class EXP OrbitIntegrator<coord::Car>;
template class EXP OrbitIntegrator<coord::Cyl>;
template class EXP OrbitIntegrator<coord::Sph>;
template class EXP RuntimeTrajectory<coord::Car>;
template class EXP RuntimeTrajectory<coord::Cyl>;
template class EXP RuntimeTrajectory<coord::Sph>;
template EXP coord::PosVelCar integrate(const coord::PosVelCar&,
    const double, const math::IOdeSystem&, const RuntimeFncArray&, const OrbitIntParams&);
template EXP coord::PosVelCyl integrate(const coord::PosVelCyl&,
    const double, const math::IOdeSystem&, const RuntimeFncArray&, const OrbitIntParams&);
template EXP coord::PosVelSph integrate(const coord::PosVelSph&,
    const double, const math::IOdeSystem&, const RuntimeFncArray&, const OrbitIntParams&);

}  // namespace orbit