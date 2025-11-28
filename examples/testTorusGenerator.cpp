/*Code to generate Fig. 7 of Binney et al 2026: surface of section of a
 *perfect ellipsoid. Data written to testTG.dat. See aloso TG.log for
 *diagnostics.
 */

#include "potential_perfect_ellipsoid.h"
#include "potential_composite.h"
#include "potential_interpolators.h"
#include "potential_utils.h"
#include "actions_staeckel.h"
#include "actions_newtorus.h"

using namespace potential;

double I3Stack(OblatePerfectEllipsoid pot, coord::PosVelCyl p){
	double E = totalEnergy(pot,p), Lz2 = pow_2(p.R*p.vphi);
	double D2 = pot.coordsys().Delta2, Glambda;
	coord::PosVelProlSph lnu(coord::toPosVel(p, pot.coordsys()));
	double absnu=fabs(lnu.nu);
	pot.evalDeriv(lnu.lambda, &Glambda);
	return lnu.lambda*(E-.5*Lz2/(lnu.lambda-D2)+Glambda
		-.125*pow_2(lnu.lambdadot*(lnu.lambda-absnu))/((lnu.lambda-D2)*pow_2(lnu.lambda)));
}
double I3val(const OblatePerfectEllipsoid& pot, const double E,
	     const double R, const double Lz, double& Phi){
	pot.eval(coord::PosCyl(R,0,0), &Phi);
	double vR=0, vphi=Lz/R, vz=sqrt(2*(E-Phi)-pow_2(vphi)); 
	coord::PosVelCyl Rzp(R,0,0,vR,vz,vphi);
	return I3Stack(pot,Rzp);
}
double lambdadot(OblatePerfectEllipsoid pot, coord::PosVelCyl p){
	double E = totalEnergy(pot,p), Lz2 = pow_2(p.R*p.vphi);
	double D2 = pot.coordsys().Delta2, Glambda;
	coord::PosVelProlSph lnu(coord::toPosVel(p, pot.coordsys()));
	double absnu = fabs(lnu.nu);
	pot.evalDeriv(lnu.lambda, &Glambda);
	double I3 = (lnu.lambda-D2)*(E-.5*Lz2/(lnu.lambda-D2)+Glambda
	-.125*lnu.lambdadot*(lnu.lambda-absnu)/((lnu.lambda-D2)*pow_2(lnu.lambda)));	pot.evalDeriv(lnu.lambda, &Glambda);
	return 8*(lnu.lambda-D2)*lnu.lambda/(lnu.lambda-absnu)*(lnu.lambda*(E-.5*Lz2/(lnu.lambda-D2)
		+Glambda)-I3);
}
// helper class to find the initial conditions corresponding to the given triplet of actions;
// the point is placed at a radius R with vertical velocity vz, and should correspond to the pericentre.
class RFinder: public math::IFunctionNoDeriv {
	const OblatePerfectEllipsoid& pot;
	const double E, I3, Lz;
	public:
		RFinder(const OblatePerfectEllipsoid& _pot, 
			 const double _E, const double _I3, const double _Lz):
		    pot(_pot), E(_E), I3(_I3), Lz(_Lz) {
		}
		virtual double value(const double R) const
		{
			double Phi; pot.eval(coord::PosCyl(R,0,0), &Phi);
			double vz = 2*(E-Phi) - pow_2(Lz/R);
			vz = vz>0? sqrt(vz) : 0;
			coord::PosVelCyl point(R, 0, 0, 0, vz, Lz/R);
			return I3 - I3Stack(pot, point);
		}
};
class vRFinder: public math::IFunctionNoDeriv {
	const OblatePerfectEllipsoid& pot;
	const double R, EmP, I3, Lz;
	public:
		vRFinder(const OblatePerfectEllipsoid& _pot, 
			 const double _R, const double _EmP, const double _I3, const double _Lz):
		    pot(_pot), R(_R), EmP(_EmP), I3(_I3), Lz(_Lz) {
		}
		virtual double value(const double vR) const
		{
			double vz = 2*EmP - pow_2(vR)-pow_2(Lz/R);
			vz = vz>0? sqrt(vz) : 0;
			coord::PosVelCyl point(R, 0, 0, vR, vz, Lz/R);
			return I3 - I3Stack(pot, point);
		}
};
void getCurve(const OblatePerfectEllipsoid& pot,int np,
	      double E, double Lz, double I3, double Rmin, double Rmax,
	      std::vector<double>& Rs, std::vector<double>& vRs){
	double tol=1e-6;
	for(int i=0; i<np; i++){
		double R=.5*(Rmin+Rmax) + .5*(Rmax-Rmin)*cos(i*M_PI/(double)(np-1));
		if(i==0) R -= 1e-10;
		double Phi; pot.eval(coord::PosCyl(R,0,0), &Phi);
		vRFinder vF(pot, R, E-Phi, I3, Lz);
		double vmax = sqrt(2*(E-Phi)-pow_2(Lz/R));
		double bot = vF.value(0), top=vF.value(vmax);
		if(bot*top>0 || std::isnan(vmax)) break;
		double vR=findRoot(vF, -1e-10, vmax, tol);
		Rs.push_back(R); vRs.push_back(vR);
	}
}

int main(void) {
	FILE* ofile; fopen_s(&ofile,"testTG.dat","w");
	PtrOblatePerfectEllipsoid
			PEpot(new OblatePerfectEllipsoid(1, 1, 0.6));
	PtrPotential ptr(new OblatePerfectEllipsoid(1, 1, 0.6));
	PtrPotential pot(wrapPotential(ptr));
	double tol=1e-6;
	int np=150;
	std::vector<double> Rs,vRs;
	actions::TorusGenerator tg(*pot,tol,"TG.log");
	double Phi0 = pot->value(coord::PosCyl(0,0,0));
	double E = .5*Phi0, Lz = 0.05;
	double Rpmin, Rpmax;
	findPlanarOrbitExtent(*pot,E,Lz,Rpmin,Rpmax);
	printf("Rmin/max %8.5f %8.5f\n",Rpmin,Rpmax);
	double Phi, I3zero = I3val(*PEpot,E,Rpmax-1e-10,Lz,Phi);//for zero-v curve
	getCurve(*PEpot,np,E,Lz,I3zero,Rpmin,Rpmax-1e-10,Rs,vRs);
	double vRmax = 0, I3trans;
	for(int i=0; i<vRs.size(); i++) vRmax=fmax(vRmax,vRs[i]);
	fprintf(ofile,"(R,vR) of zero-velocity curve at Jphi = %8.4f)\n", Lz);
	for(int i=0; i<Rs.size(); i++)
		fprintf(ofile,"%8.5f\t%8.5f\n",Rs[i],vRs[i]);
	int nc=11;
	for(int j=0; j<nc; j++){
		double Rmin=1.1*Rpmin, Rmax=.95*Rpmax, R=Rmax;
		double I3 = I3val(*PEpot, E, R, Lz, Phi);
		RFinder RF(*PEpot, E, I3, Lz);
		double Rout=R-1e-8, Rin=.9*Rout, outer=RF.value(.95*Rout);
		while(RF.value(Rin)*outer > 0 && Rin>1e-5) Rin *= .9;
		if(Rin>1.2e-5) Rin = math::findRoot(RF,Rin,Rout,tol);
		getCurve(*PEpot,np,E,Lz,I3,Rin,Rout,Rs,vRs);
		bool onlyTM=false;//true;//
		double vR=0, vphi=Lz/R, vz=sqrt(2*(E-Phi)-pow_2(vphi)); 
		coord::PosVelCyl Rzp(R,0,0,vR,vz,vphi);
		actions::Frequencies omegas;
		actions::ActionAngles aa(actions::actionAnglesAxisymStaeckel(*PEpot, Rzp, &omegas));
		printf("Doing J: (%8.5f %8.5f) omegas %8.5f %8.5f %8.5f %8.5f\n",aa.Jr,aa.Jz,
		       omegas.Omegar,omegas.Omegaz,omegas.Omegaphi,omegas.Omegaz/omegas.Omegaphi);
		actions::Actions J(aa.Jr,aa.Jz,aa.Jphi);
		actions::Torus T, BT;
		T=tg.fitTorus(J);
		printf("E %8.5f Jf %8.5f,Jzcrit %8.5f\n",T.E,2*J.Jr+J.Jz,pot->getJzcrit(2*J.Jr+J.Jz));
		double params[5];
		T.ptrTM->getPointTrans()->getParams(params);
		BT=tg.giveBaseTorus(J, T.ptrTM);
		std::vector<double> RsT, vRsT, RsBT, vRsBT;
		double rmin,rmax,vzmax;
		T.zSoS(RsT,vRsT,2*np,rmin,rmax,vzmax);
		BT.zSoS(RsBT,vRsBT,2*np,rmin,rmax,vzmax);
		fprintf(ofile,"(R,vR) for torus and toy map at (Jr,Jz) = (%8.4f, %8.2f)\n",
			  J.Jr, J.Jz);
		for(int i=0; i<RsT.size(); i++)
			fprintf(ofile,"%8.5f\t%8.5f\t%8.5f\t%8.5f\n",
				  RsT[i],vRsT[i],RsBT[i],vRsBT[i]);
		Rs.clear(); vRs.clear();
		double diff=.054*(Rpmax-Rpmin);
		Rpmin+=.1*diff;
		Rpmax-=diff;
	}
}