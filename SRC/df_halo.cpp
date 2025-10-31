#include "df_halo.h"
#include "math_core.h"
#include "math_ode.h"
#include <cmath>
#include <stdexcept>

//#define HENON

namespace df{

namespace {  // internal ns

/* Class used by NewOxford DF for anisotropic spherical systems */
class NewOxford_r_syst: public math::IOdeSystem{
	private:
		const double r_rat, beta;
	public:
		NewOxford_r_syst(double _r_rat, double _beta) :
		    r_rat(_r_rat), beta(_beta) {};
		void eval(const double t,const double s[],double dsdt[]) const{
			//returns Om_r/Om_phi as function of circularity
			double c = s[1]>0? s[0]/(s[0]+s[1]) : 1;//circularity
			double psi = c*.5*M_PI;
			double top = beta!=0? exp(-beta*sin(psi)) : 1; 
			dsdt[0]=top/(.51+c*(r_rat-.51));//r_rat=Om/kappa
			dsdt[1]=-1;//t=const-Jr
		}
		unsigned int size()  const {return 2;}
};

/*Class used to infer H(Jr,L) for a Plummer potential */
class y_getter : public math::IFunction {
	const double Jr, L, L0;
	double cL, Ebot;
	const math::LinearInterpolator& jrmax;
	const math::LinearInterpolator& cls;
	public:
		y_getter(const double _Jr,const double _L,const double _L0,
			 const math::LinearInterpolator& _jrmax,
			 const math::LinearInterpolator& _cls):
		    Jr(_Jr), L(_L), L0(_L0), jrmax(_jrmax), cls(_cls){
			cL=(.5*(L+sqrt(L*L+L0*L0)));
			Ebot=(fmax(jrmax.evalReverse(Jr),cls.evalReverse(cL)));
			//printf("cL,Ebot %g %g\n",cL,Ebot);
		}
		double Eb(){
			return Ebot;
		}
		virtual void evalDeriv(const double E,double* val,double* deriv=NULL, double* deriv2=NULL) const{
			double Jrmax,cLm;
			jrmax.evalDeriv(E,&Jrmax); cls.evalDeriv(E,&cLm);
			if(val) *val = Jr/Jrmax + (cL-cLm)/(cLm-.5*L0);
		}
		virtual unsigned int numDerivs() const{
			return 0;
		}
};

/* helper class used in the root-finder to determine the auxiliary
 coefficient beta for a cored halo of DoublePowerLaw DF*/
class BetaFinderNDPL: public math::IFunctionNoDeriv{
	const DoublePowerLawParam& par;

    // return the difference between the non-modified and modified DF as a function of beta and
    // the appropriately scaled action variable (t -> hJ), weighted by d(hJ)/dt for the integration in t
	double deltaf(const double t, const double beta) const
	{
	// integration is performed in a scaled variable t, ranging from 0 to 1,
	// which is remapped to hJ ranging from 0 to infinity as follows:
		double hJ    = par.Jcore * t*t*(3-2*t) / pow_2(1-t) / (1+2*t); 
		double dhJdt = par.Jcore * 6*t / pow_3(1-t) / pow_2(1+2*t);
		return hJ * hJ * dhJdt *
				math::pow(1 + math::pow(par.J0 / hJ, par.steepness),  par.slopeIn  / par.steepness) *
				math::pow(1 + math::pow(hJ / par.J0, par.steepness), -par.slopeOut / par.steepness) *
				(math::pow(1 + par.Jcore/hJ * (par.Jcore/hJ - beta), -0.5*par.slopeIn) - 1);
	}

	public:
		BetaFinderNDPL(const DoublePowerLawParam& _par) : par(_par) {}

		virtual double value(const double beta) const
		{
			double result = 0;
	// use a fixed-order GL quadrature to compute the integrated difference in normalization between
	// unmodified and core-modified DF, which is sought to be zero by choosing an appropriate beta
			static const int GLORDER = 20;  // should be even, to avoid singularity in the integrand at t=0.5
			for(int i=0; i<GLORDER; i++)
				result += math::GLWEIGHTS[GLORDER][i] * deltaf(math::GLPOINTS[GLORDER][i], beta);
			return result;
		}
};

/* helper function to compute the auxiliary coefficient beta in the case of a
 * central core */
double computeBetaNDPL(const DoublePowerLawParam &par)
{
	if(par.Jcore<=0)
		return 0;
	return math::findRoot(BetaFinderNDPL(par), 0.0, 2.0, /*root-finder tolerance*/ SQRT_DBL_EPSILON);
}

/* helper class used in the root-finder to determine the auxiliary
 coefficient beta for a cored halo of dSph DF*/
class BetaFinderNdSph: public math::IFunctionNoDeriv{
	const DoublePowerLawParam& par;

    // return the difference between the non-modified and modified DF as a function of beta and
    // the appropriately scaled action variable (t -> hJ), weighted by d(hJ)/dt for the integration in t
	double deltaf(const double t, const double beta) const
	{
	// integration is performed in a scaled variable t, ranging from 0 to 1,
	// which is remapped to hJ ranging from 0 to infinity as follows:
		double hJ    = par.Jcore * t*t*(3-2*t) / pow_2(1-t) / (1+2*t); 
		double dhJdt = par.Jcore * 6*t / pow_3(1-t) / pow_2(1+2*t);
		return hJ * hJ * dhJdt *
				math::pow(1 + math::pow(par.J0 / hJ, par.steepness),  par.slopeIn  / par.steepness) *
				math::pow(1 + math::pow(hJ / par.J0, par.steepness), -par.slopeOut / par.steepness) *
				(math::pow(1 + par.Jcore/hJ * (par.Jcore/hJ - beta), -0.5*par.slopeIn) - 1);
	}

	public:
		BetaFinderNdSph(const DoublePowerLawParam& _par) : par(_par) {}

		virtual double value(const double beta) const
		{
			double result = 0;
	// use a fixed-order GL quadrature to compute the integrated difference in normalization between
	// unmodified and core-modified DF, which is sought to be zero by choosing an appropriate beta
			static const int GLORDER = 20;  // should be even, to avoid singularity in the integrand at t=0.5
			for(int i=0; i<GLORDER; i++)
				result += math::GLWEIGHTS[GLORDER][i] * deltaf(math::GLPOINTS[GLORDER][i], beta);
			return result;
		}
};

/* helper function to compute the auxiliary coefficient beta in the case of a
 * central core */
double computeBetaNdSph(const DoublePowerLawParam &par)
{
	if(par.Jcore<=0)
		return 0;
	return math::findRoot(BetaFinderNdSph(par), 0.0, 2.0, /*root-finder tolerance*/ SQRT_DBL_EPSILON);
}


}  // internal ns

EXP DoublePowerLaw::DoublePowerLaw(const DoublePowerLawParam &inparams) :
    par(inparams), beta(computeBetaNDPL(par))
{
    // sanity checks on parameters
	if(!(par.norm>0))
		throw std::invalid_argument("DoublePowerLaw: normalization must be positive");
	if(!(par.J0>0))
		throw std::invalid_argument("DoublePowerLaw: break action J0 must be positive");
	if(!(par.Jcore>=0 && beta>=0))
		throw std::invalid_argument("DoublePowerLaw: core action Jcore is invalid");
	if(!(par.Jcutoff>=0))
		throw std::invalid_argument("DoublePowerLaw: truncation action Jcutoff must be non-negative");
	if(!(par.slopeOut>3) && par.Jcutoff==0)
		throw std::invalid_argument(
					    "DoublePowerLaw: mass diverges at large J (outer slope must be > 3)");
	if(!(par.slopeIn<3))
		throw std::invalid_argument(
					    "DoublePowerLaw: mass diverges at J->0 (inner slope must be < 3)");
	if(!(par.steepness>0))
		throw std::invalid_argument("DoublePowerLaw: transition steepness parameter must be positive");
	if(!(par.cutoffStrength>0))
		throw std::invalid_argument("DoublePowerLaw: cutoff strength parameter must be positive");
	if(!(par.coefJrIn>0 && par.coefJzIn >0 && par.coefJrIn +par.coefJzIn <3 &&
	     par.coefJrOut>0 && par.coefJzOut>0 && par.coefJrOut+par.coefJzOut<3) )
		throw std::invalid_argument(
					    "DoublePowerLaw: invalid weights in the linear combination of actions");
	if(!(fabs(par.rotFrac)<=1))
		throw std::invalid_argument(
					    "DoublePowerLaw: amplitude of odd-Jphi component must be between -1 and 1");
	if(par.Fname!=""){
		readBrighterThan(par.Fname);
	}
}

EXP double DoublePowerLaw::value(const actions::Actions &J, const double Jzc) const
{
	double L=J.Jz+fabs(J.Jphi);
    // linear combination of actions in the inner part of the model (for J<J0)
	double hJ  = fmax(0, par.coefJrIn * J.Jr +par.coefJzIn * J.Jz +fabs(J.Jphi));
    // linear combination of actions in the outer part of the model (for J>J0)
	double gJ  = fmax(0, par.coefJrOut * J.Jr + par.coefJzOut * J.Jz + fabs(J.Jphi));
	double val = par.norm / pow_3(2*M_PI * par.J0) *
		     math::pow(1 + math::pow(par.J0 / hJ, par.steepness),  par.slopeIn  / par.steepness) *
		     math::pow(1 + math::pow(gJ / par.J0, par.steepness), -par.slopeOut / par.steepness);
	if(par.rotFrac!=0)  // add the odd part
		val *= 1 + par.rotFrac * tanh(J.Jphi / par.Jphi0);
	if(par.Jcutoff>0)   // exponential cutoff at large J
		val *= exp(-math::pow(gJ / par.Jcutoff, par.cutoffStrength));
	if(par.Jcore>0) {   // central core of nearly-constant f(J) at small J
		if(hJ==0) return par.norm / pow_3(2*M_PI * par.J0);
		val *= math::pow(1 + par.Jcore/hJ * (par.Jcore/hJ - beta), -0.5*par.slopeIn);
	}
	return val;
}

NewDoublePowerLaw::NewDoublePowerLaw(const DoublePowerLawParam& _params,
				     const potential::BasePotential& _pot) :
    epsilonJ(_params.epsilonJ), DF0(_params) {}
double NewDoublePowerLaw::wt(const actions::Actions& J) const{
	return 1/(1+pow_2(J.Jphi)/(epsilonJ*(epsilonJ+J.Jr+J.Jz)));
}
double NewDoublePowerLaw::value(const actions::Actions& J, const double Jzc) const{
	double w = wt(J), sgn = J.Jphi>0? 1 : -1;
	actions::Actions J1;
	double eps = .2*epsilonJ;
	if(J.Jz < Jzc){
		J1=actions::Actions(J.Jr+.5*(fabs(J.Jphi)-eps), J.Jz, eps);
	} else {
		J1=actions::Actions(J.Jr, J.Jz+(fabs(J.Jphi)-eps), eps);
	}
	double F0=DF0.value(J, Jzc), F1=DF0.value(J1, Jzc);
	return w*F1 + (1-w)*F0;
}

DwarfSpheroid::DwarfSpheroid(const DwarfSpheroidParam& _par): par(_par){}
	
double DwarfSpheroid::value(const actions::Actions &J, const double Jzc) const {
	double fJphi=fabs(J.Jphi);
	double jt = fmax(0, (par.coefJr*J.Jr + par.coefJz*J.Jz + fJphi)/par.J0);
	double val = par.norm * exp(-pow( jt, par.alpha));
	if(par.rotFrac!=0)  // add the odd part
		val *= 1 + par.rotFrac * tanh(J.Jphi / par.Jphi0);
	return val; 
}

NewDwarfSpheroid::NewDwarfSpheroid(const DwarfSpheroidParam& _par,
				     const potential::BasePotential& _pot) :
    epsilonJ(_par.epsilonJ), DF0(_par) {}
double NewDwarfSpheroid::wt(const actions::Actions& J) const{
	return 1/(1+pow_2(J.Jphi)/(epsilonJ*(epsilonJ+J.Jr+J.Jz)));
}
double NewDwarfSpheroid::value(const actions::Actions& J, const double Jzc) const{
	double w = wt(J), sgn = J.Jphi>0? 1 : -1;
	double eps = .2*epsilonJ;
	actions::Actions J1;
	if(J.Jz < Jzc){
		J1=actions::Actions(J.Jr+.5*(fabs(J.Jphi)-eps), J.Jz, eps);
	} else {
		J1=actions::Actions(J.Jr, J.Jz+(fabs(J.Jphi)-eps), eps);
	}
	double F0 = DF0.value(J, Jzc), F1 = DF0.value(J1, Jzc);
	return w*F1 + (1-w)*F0;
}

EXP double NewOxford::h(const actions::Actions &J) const{
	//Move over E surface to circular orbit then to circular orbit &
	//return Lc
	const double L=J.Jz+fabs(J.Jphi);
	double rats[2];	freq.epicycle_ratios(1.4*J.Jr+L, rats);//rats{Om/kappa,nu/Om} at this E
	double ic[2]={L/par.J0, J.Jr/par.J0};
	NewOxford_r_syst r_system(rats[0], par.beta);
	math::OdeSolverDOP853 r_solver(r_system);
	r_solver.init(ic);
	double t=0; int j=0; double jr=ic[1];
	const int nmax=30;
	while(j<nmax && jr>0){//until Jr=0
		r_solver.doStep(0); t=r_solver.getTime(); jr=r_solver.getSol(t,1);
		j++;
	}
	double jp = j==0? ic[0] : r_solver.getSol(ic[1],0);// Jphi/L0 when Jr=0
	if(j==nmax){
		double s[2]={r_solver.getSol(t,0),r_solver.getSol(t,1)};
		printf("r_solver error: %g (%g %g) (%g %g)\n",t,ic[0],ic[1],s[0],s[1]);
		exit(0);
	}
	return par.J0*jp;//return Jphi when Jr=0
}

EXP double NewOxford::value(const actions::Actions &J, const double Jzc) const{
	double hJ=h(J);
	double val =  par.norm / pow_3(2*M_PI * par.J0) *
		      math::pow(1 + math::pow(par.J0 / hJ, par.steepness),  par.slopeIn  / par.steepness) *
		      math::pow(1 + math::pow(hJ / par.J0, par.steepness), -par.slopeOut / par.steepness);
	if(par.rotFrac!=0)  // add the odd part
		val *= 1 + par.rotFrac * tanh(J.Jphi / par.Jphi0);
	if(par.Jcutoff>0)   // exponential cutoff at large J
		val *= exp(-math::pow(hJ / par.Jcutoff, par.cutoffStrength));
	return val;
}

EXP PlummerDF::PlummerDF(const PlummerParam& params): par(params){
	if(par.mass<=0) printf("PlummerDF: mass must be >0\n");
	if(par.scaleRadius <= 0) printf("PlummerDF: scaleRadius must be >0\n");
	FILE* ifile;
	if(fopen_s(&ifile,"/u/gmb/c/PlummerEs.dat","r"))
		printf("PlummerDF: I can't find jrmax, cL data\n");
	std::vector<double> Es,Jrmax,cLs;
	while(!feof(ifile)){	//load Lc,Jrmax
		double x,y,z;
		if(3!=fscanf_s(ifile,"%lf %lf %lf",&x,&y,&z)) break;
		Es.push_back(x); cLs.push_back(y); Jrmax.push_back(z);
		//printf("%f %f %f\n",Es.back(),cLs.back(),Jrmax.back());
	}
	fclose(ifile);
	jrmax = math::LinearInterpolator(Es,Jrmax);
	cls =   math::LinearInterpolator(Es,cLs);
	Etop=Es.back(); Ebot=Es[0];
	Jrtop=Jrmax.back(); cLtop=cLs.back();
	Jrbot=Jrmax[0]; cLbot=cLs[0];
	norm =1;
	norm = par.mass/totalMass();	printf("norm %g\n",norm);
	if(norm<1e-10) exit(0);
}
EXP void PlummerDF::write_params(std::ofstream& strm, const units::InternalUnits& intUnits) const{
	strm << "type\t = PlummerDF\n";
	strm << "mass\t = " << par.mass << "\n";
	strm << "scaleRadius\t = " << par.scaleRadius << "\n";
	strm << "scaleAction\t = " << par.scaleAction << "\n";
	strm << "mu = " << par.mu << " nu = " << par.nu << "\n";
}
EXP double PlummerDF::value(const actions::Actions& J, const double Jzc) const{
	double L=J.Jz+fabs(J.Jphi),Jmod=J.Jr+L;
	double sq0=L*L+2*pow_2(par.scaleAction), muL=par.mu*L*L/sq0, emu=exp(muL);
	double Jrp=J.Jr/emu, Lp=L*emu, H;
	if(Jrp>Jrtop || Lp>cLtop)//Kepler limit
		H=-.5*pow_2(par.mass/(Jrp+Lp+.25*pow_2(par.scaleAction)/Lp));
	else if(Jrp<Jrbot || Lp<cLbot)//HO limit
		H=-pow_2(par.mass*.5*par.scaleAction/(.5*par.scaleAction+Jrp+.5*Lp*(1+.5*Lp/par.scaleAction)));
	else{
		y_getter yg(Jrp, Lp, par.scaleAction, jrmax, cls);
		H = math::findRoot(yg,yg.Eb(),Etop,1e-5);
		if(Jmod>4 && H<-.1){
			double tr,Jrm,cLm;  yg.evalDeriv(H,&tr);
			jrmax.evalDeriv(H,&Jrm); cls.evalDeriv(H,&cLm);
//			double dum=Jrp/Jrm+(cL-cLm)/(cLm-.5*par.scaleAction);
			printf("PlummerDF error: %g %g %g %g %g\n",H,tr,Jrp,L,Lp);// exit(0);
			printf("%g %g %g\n",Jrp/Jrm,cLm,(.5*(Lp+sqrt(Lp*Lp+pow_2(par.scaleAction)))-cLm)/(cLm-.5*par.scaleAction));
		}
		if(std::isnan(H) || H>=-1e-3){
			if(J.Jr+L<.1) H=-1;
			else return 0;
		}
	}
	return norm*pow(-H,3.5);
}
EXP IsochroneDF::IsochroneDF(const IsochroneParam& params) : par(params){
	if(par.mass<=0) printf("IsochronDF: mass must be >0\n");
	if(par.scaleRadius <= 0) printf("IsochronDF: scaleRadius must be >0\n");
	norm =1; norm=par.mass/totalMass();
}
EXP void IsochroneDF::write_params(std::ofstream& strm, const units::InternalUnits& intUnits) const{
	strm << "type\t = IsochroneDF\n";
	strm << "mass\t = " << par.mass << "\n";
	strm << "scaleRadius\t = " << par.scaleRadius << "\n";
	strm << "mu\t = " << par.mu << " nu = " << par.nu << "\n";
}
#ifdef HENON //DF as given in Henon's paper with obscure scaling
EXP double IsochroneDF::value(const actions::Actions& J, const double Jzc) const{
	const double fac=sqrt(2)/(pow_3(2*M_PI));
	double L=J.Jz+fabs(J.Jphi), J0=J.Jr+.5*(L+sqrt(L*L+4*par.mass*par.scaleRadius));
	double Et=-.5*pow_2(par.mass/J0)*par.scaleRadius/(par.mass);
	if(Et>=0) return 0;
	double U=(1+Et);
	return fac*((48+U*(-53+U*(14+U*(14+U*4))))*sqrt(1-U)
		    +3*(9+U*(-22+U*4))*acos(U)/sqrt(U+1))/pow(U+1,4);
}
#else
EXP double IsochroneDF::value(const actions::Actions& J, const double Jzc) const{
	const double fac=norm/(sqrt(2)*pow_3(2*M_PI));
	double L=J.Jz+fabs(J.Jphi), L02=4*par.mass*par.scaleRadius;
	double sq0=L*L+1.5*L02, muL=par.mu*L*L/sq0, emu=exp(muL);
	double Lp=L*emu;
	double J0 = J.Jr/emu + .5*(Lp + sqrt(Lp*Lp+L02));
	double Et=-.5*pow_2(par.mass/J0)*par.scaleRadius/par.mass;
	if(Et>=0) return 0;
	double fac1=fac/(pow(par.mass*par.scaleRadius,1.5)*pow(2*(1+Et),4));
	return fac1*(sqrt(-Et)*(27+Et*(66+Et*(320+Et*(240+Et*64))))+3*((16*Et-28)*Et-9)*
		     asin(sqrt(-Et))/sqrt(1+Et));
}
#endif
}  // namespace df
