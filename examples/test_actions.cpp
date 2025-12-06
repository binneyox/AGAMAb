#include "actions_staeckel.h"
#include "actions_newtorus.h"
#include "actions_spherical.h"
#include "potential_factory.h"
#include "potential_multipole.h"
#include "potential_composite.h"
#include "potential_cylspline.h"
#include "potential_analytic.h"
#include "potential_perfect_ellipsoid.h"
#include "galaxymodel_base.h"
#include "orbit.h"
#include "map.h"
//#include <Eigen/Dense>
#include <variant>
#include <random>
#include <vector>
#include <iomanip>
#include <chrono>
#include <iostream>
#include <format> 
#include <fstream>
#include <stdio.h>

std::vector<double> minmaxelement(std::vector<double> A) {
    double max = A[0];
    double min = A[0];
    for (int i = 1;i < A.size();i++) {
        if (A[i] > max) {
            max = A[i];
        }
        else if (A[i] < min) {
            min = A[i];
        }
    }
    std::vector<double> v(2, 0.0);
    v[0] = min, v[1] = max;
    return v;
}
void checkJDerivs(actions::ActionAngles aa, actions::PtrToyMap TM2) {
    double h = 1e-5;
    actions::DerivAct<coord::Cyl> dJ;
    actions::DerivAng<coord::Cyl> dA;
    coord::PosMomCyl Rz=TM2->from_aaT(aa,dJ);
    printf("R:(%f,%f,%f,%f,%f,%f)\n",Rz.R,Rz.z,Rz.phi,Rz.pR,Rz.pz,Rz.pphi);
    TM2->from_aaT(aa, dA);
    aa.Jr -= .5 * h;
    coord::PosMomCyl xv1 = TM2->from_aaT(aa);
    aa.Jr += h;
    coord::PosMomCyl xv2 = TM2->from_aaT(aa);
    printf("dRndJr:%f %f, dzndJr:%f %f, dphidJr:%f %f, dpRndJr:%f %f, dpzndJr:%f %f, dpphidJr:%f %f\n",
        (xv2.R - xv1.R) / h,dJ.dbyJr.R, (xv2.z - xv1.z) / h, dJ.dbyJr.z,(xv2.phi-xv1.phi)/h,dJ.dbyJr.phi,
        (xv2.pR - xv1.pR) / h, dJ.dbyJr.pR, (xv2.pz - xv1.pz) / h,dJ.dbyJr.pz, dJ.dbyJr.pphi,0.0);
    aa.Jr -= .5 * h;
    aa.Jz -= .5 * h;
    xv1 = TM2->from_aaT(aa);
    aa.Jz += h;
    xv2 = TM2->from_aaT(aa);
    printf("dRndJz:%f %f, dzndJz:%f %f, dphidJz:%f %f, dpRndJz:%f %f, dpzndJz:%f %f, dpphidJz:%f %f\n",
        (xv2.R - xv1.R) / h, dJ.dbyJz.R, (xv2.z - xv1.z) / h, dJ.dbyJz.z, (xv2.phi - xv1.phi) / h, dJ.dbyJz.phi,
        (xv2.pR - xv1.pR) / h, dJ.dbyJz.pR, (xv2.pz - xv1.pz) / h,dJ.dbyJz.pz, dJ.dbyJz.pphi, 0.0);
    aa.Jz -= .5 * h;
   // aa.Jphi -= .5 * h;
    xv1 = TM2->from_aaT(aa);
    aa.Jphi += h;
    xv2 = TM2->from_aaT(aa);
    printf("dRndJphi:%f %f, dzndJphi:%f %f, dphidJphi:%f %f, dpRndJphi:%f %f, dpzndJphi:%f %f, dpphidJphi:%f %f\n",
        (xv2.R - xv1.R) / h, dJ.dbyJphi.R, (xv2.z - xv1.z) / h, dJ.dbyJphi.z, (xv2.phi - xv1.phi) / h, dJ.dbyJphi.phi,
        (xv2.pR - xv1.pR) / h, dJ.dbyJphi.pR, (xv2.pz - xv1.pz) / h, dJ.dbyJphi.pz, dJ.dbyJphi.pphi, 1.0);
    aa.Jphi -= h;
    aa.thetar -= .5 * h;
    xv1 = TM2->from_aaT(aa);
    aa.thetar += h;
    xv2 = TM2->from_aaT(aa);
    printf("dRndthetar:%f %f, dzndthetar:%f %f, dphidthetar:%f %f, dpRndthetar:%f %f, dpzndthetar:%f %f, dpphidthetar:%f %f\n",
        (xv2.R - xv1.R) / h, dA.dbythetar.R, (xv2.z - xv1.z) / h, dA.dbythetar.z, (xv2.phi - xv1.phi) / h, dA.dbythetar.phi,
        (xv2.pR - xv1.pR) / h, dA.dbythetar.pR, (xv2.pz - xv1.pz) / h, dA.dbythetar.pz, dA.dbythetar.pphi, 0.);
    aa.thetar -= .5 * h;
    aa.thetaz -= .5 * h;
    xv1 = TM2->from_aaT(aa);
    aa.thetaz += h;
    xv2 = TM2->from_aaT(aa);
    printf("dRndthetaz:%f %f, dzndthetaz:%f %f, dphidthetaz:%f %f, dpRndthetaz:%f %f, dpzndthetaz:%f %f dpphidthetaz:%f %f\n",
        (xv2.R - xv1.R) / h, dA.dbythetaz.R, (xv2.z - xv1.z) / h, dA.dbythetaz.z, (xv2.phi - xv1.phi) / h, dA.dbythetaz.phi,
        (xv2.pR - xv1.pR) / h, dA.dbythetaz.pR, (xv2.pz - xv1.pz) / h, dA.dbythetaz.pz, dA.dbythetaz.pphi, 0.);
    aa.thetaz -= .5 * h;
    aa.thetaphi -= .5 * h;
    xv1 = TM2->from_aaT(aa);
    aa.thetaphi += h;
    xv2 = TM2->from_aaT(aa);
    printf("dRndthetaphi:%f %f, dzndthetaphi:%f %f, dphidthetaphi:%f %f, dpRndthetaphi:%f %f, dpzndthetaphi:%f %f dpphidthetaphi:%f %f\n",
        0., dA.dbythetaphi.R, 0., dA.dbythetaphi.z, (xv2.phi - xv1.phi) / h, dA.dbythetaphi.phi,
        0., dA.dbythetaphi.pR, 0., dA.dbythetaphi.pz, dA.dbythetaphi.pphi, 0.);
    aa.thetaz -= .5 * h;
}

void checkJDerivs2(actions::ActionAngles aa, actions::PtrToyMap TM2) {
    double h = 1e-5;
    actions::DerivAct<coord::Cyl> dJ;
    actions::DerivAng<coord::Cyl> dA;
    coord::PosMomCyl Rz=TM2->from_aaT(aa,dJ);
    printf("R:(%f,%f,%f,%f,%f,%f)\n",Rz.R,Rz.z,Rz.phi,Rz.pR,Rz.pz,Rz.pphi);
    TM2->from_aaT(aa, dA);
    aa.Jr -= .5 * h;
    coord::PosMomCyl xv1 = TM2->from_aaT(aa);
    aa.Jr += h;
    coord::PosMomCyl xv2 = TM2->from_aaT(aa);
    printf("dRndJr:%f %f, dzndJr:%f %f, dphidJr:%f %f, dpRndJr:%f %f, dpzndJr:%f %f, dpphidJr:%f %f\n",
        (xv2.R - xv1.R) / h,dJ.dbyJr.R, (xv2.z - xv1.z) / h, dJ.dbyJr.z,(xv2.phi-xv1.phi)/h,dJ.dbyJr.phi,
        (xv2.pR - xv1.pR) / h, dJ.dbyJr.pR, (xv2.pz - xv1.pz) / h,dJ.dbyJr.pz, dJ.dbyJr.pphi,0.0);
    aa.Jr -= .5 * h;
    aa.Jz -= .5 * h;
    xv1 = TM2->from_aaT(aa);
    aa.Jz += h;
    xv2 = TM2->from_aaT(aa);
    printf("dRndJz:%f %f, dzndJz:%f %f, dphidJz:%f %f, dpRndJz:%f %f, dpzndJz:%f %f, dpphidJz:%f %f\n",
        (xv2.R - xv1.R) / h, dJ.dbyJz.R, (xv2.z - xv1.z) / h, dJ.dbyJz.z, (xv2.phi - xv1.phi) / h, dJ.dbyJz.phi,
        (xv2.pR - xv1.pR) / h, dJ.dbyJz.pR, (xv2.pz - xv1.pz) / h,dJ.dbyJz.pz, dJ.dbyJz.pphi, 0.0);
    aa.Jz -= .5 * h;
   // aa.Jphi -= .5 * h;
    xv1 = TM2->from_aaT(aa);
    aa.Jphi += h;
    xv2 = TM2->from_aaT(aa);
    printf("dRndJphi:%f %f, dzndJphi:%f %f, dphidJphi:%f %f, dpRndJphi:%f %f, dpzndJphi:%f %f, dpphidJphi:%f %f\n",
        (xv2.R - xv1.R) / h, dJ.dbyJphi.R, (xv2.z - xv1.z) / h, dJ.dbyJphi.z, (xv2.phi - xv1.phi) / h, dJ.dbyJphi.phi,
        (xv2.pR - xv1.pR) / h, dJ.dbyJphi.pR, (xv2.pz - xv1.pz) / h, dJ.dbyJphi.pz, dJ.dbyJphi.pphi, 1.0);
    aa.Jphi -= h;
    aa.thetar -= .5 * h;
    xv1 = TM2->from_aaT(aa);
    aa.thetar += h;
    xv2 = TM2->from_aaT(aa);
    printf("dRndthetar:%f %f, dzndthetar:%f %f, dphidthetar:%f %f, dpRndthetar:%f %f, dpzndthetar:%f %f, dpphidthetar:%f %f\n",
        (xv2.R - xv1.R) / h, dA.dbythetar.R, (xv2.z - xv1.z) / h, dA.dbythetar.z, (xv2.phi - xv1.phi) / h, dA.dbythetar.phi,
        (xv2.pR - xv1.pR) / h, dA.dbythetar.pR, (xv2.pz - xv1.pz) / h, dA.dbythetar.pz, dA.dbythetar.pphi, 0.);
    aa.thetar -= .5 * h;
    aa.thetaz -= .5 * h;
    xv1 = TM2->from_aaT(aa);
    aa.thetaz += h;
    xv2 = TM2->from_aaT(aa);
    printf("dRndthetaz:%f %f, dzndthetaz:%f %f, dphidthetaz:%f %f, dpRndthetaz:%f %f, dpzndthetaz:%f %f dpphidthetaz:%f %f\n",
        (xv2.R - xv1.R) / h, dA.dbythetaz.R, (xv2.z - xv1.z) / h, dA.dbythetaz.z, (xv2.phi - xv1.phi) / h, dA.dbythetaz.phi,
        (xv2.pR - xv1.pR) / h, dA.dbythetaz.pR, (xv2.pz - xv1.pz) / h, dA.dbythetaz.pz, dA.dbythetaz.pphi, 0.);
    aa.thetaz -= .5 * h;
    aa.thetaphi -= .5 * h;
    xv1 = TM2->from_aaT(aa);
    aa.thetaphi += h;
    xv2 = TM2->from_aaT(aa);
    printf("dRndthetaphi:%f %f, dzndthetaphi:%f %f, dphidthetaphi:%f %f, dpRndthetaphi:%f %f, dpzndthetaphi:%f %f dpphidthetaphi:%f %f\n",
        0., dA.dbythetaphi.R, 0., dA.dbythetaphi.z, (xv2.phi - xv1.phi) / h, dA.dbythetaphi.phi,
        0., dA.dbythetaphi.pR, 0., dA.dbythetaphi.pz, dA.dbythetaphi.pphi, 0.);
    aa.thetaz -= .5 * h;
}
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
double I3Stack(potential::OblatePerfectEllipsoid pot, coord::PosVelCyl p){
	double E = totalEnergy(pot,p), Lz2 = pow_2(p.R*p.vphi);
	double D2 = pot.coordsys().Delta2, Glambda;
	coord::PosVelProlSph lnu(coord::toPosVel(p, pot.coordsys()));
	double absnu=fabs(lnu.nu);
	pot.evalDeriv(lnu.lambda, &Glambda);
	return lnu.lambda*(E-.5*Lz2/(lnu.lambda-D2)+Glambda
		-.125*pow_2(lnu.lambdadot*(lnu.lambda-absnu))/((lnu.lambda-D2)*pow_2(lnu.lambda)));
}
double I3val(const potential::OblatePerfectEllipsoid& pot, const double E,
	     const double R, const double Lz, double& Phi){
	pot.eval(coord::PosCyl(R,0,0), &Phi);
	double vR=0, vphi=Lz/R, vz=sqrt(2*(E-Phi)-pow_2(vphi)); 
	coord::PosVelCyl Rzp(R,0,0,vR,vz,vphi);
	return I3Stack(pot,Rzp);
}
double lambdadot(potential::OblatePerfectEllipsoid pot, coord::PosVelCyl p){
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
	const potential::OblatePerfectEllipsoid& pot;
	const double E, I3, Lz;
	public:
		RFinder(const potential::OblatePerfectEllipsoid& _pot, 
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
	const potential::OblatePerfectEllipsoid& pot;
	const double R, EmP, I3, Lz;
	public:
		vRFinder(const potential::OblatePerfectEllipsoid& _pot, 
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
void getCurve(const potential::OblatePerfectEllipsoid& pot,int np,
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
//		printf("(%f %f %f %f) ",bot,top,R,vR);
	}
}
	
int main(){
    /*potential::PtrPotential potDiskthin=potential::createPotential(
        utils::KeyValueMap("type=Disk,surfaceDensity=3512,scaleRadius=2.9,scaleHeight=0.3"));
    potential::PtrPotential potDiskthick=potential::createPotential(
        utils::KeyValueMap("type=Disk,surfaceDensity=901,scaleRadius=3,scaleHeight=0.9"));
    potential::PtrPotential potbulge=potential::createPotential(
        utils::KeyValueMap("type=spheroid,alpha=1,gamma=0,beta=1.8,scaleRadius=0.075,q=0.5,outerCutoffRadius=2.1,densityNorm=411164.7"));
    potential::PtrPotential potHalo=potential::createPotential(utils::KeyValueMap("type=spheroid,alpha=1,gamma=1,beta=3, densityNorm=36.385,q=1,scaleRadius=20.2"));
    std::vector<potential::PtrPotential> potvec={potbulge,potDiskthin,potDiskthick,potHalo};
    //potential::PtrPotential pot(new potential::CompositeCyl(potvec));*/
    potential::OblatePerfectEllipsoid PEpot(1,1,0.6);
	potential::PtrPotential pot=potential::createPotential(utils::KeyValueMap("type=PerfectEllipsoid, q=0.6, scaleRadius=1, Mass=1"));
	double tol=1e-6;
	int np=300;
	std::vector<double> Rs,vRs;
    actions::TorusGenerator tg(*pot,tol,"TG.log");
    double Phi0=pot->value(coord::PosCyl(0,0,0));
    double E = .5*Phi0, Lz = 0.05*potential::L_circ(*pot,E);
    int nc=11;
    double Rpmin, Rpmax;
	findPlanarOrbitExtent(*pot,E,Lz,Rpmin,Rpmax);
    double Phi, I3zero = I3val(PEpot,E,Rpmax-1e-10,Lz,Phi);//for zero-v curve
	for(int j=0; j<nc; j++){
		double Rmin=1.1*Rpmin, Rmax=.95*Rpmax, R=Rmax;
		double I3 = I3val(PEpot, E, R, Lz, Phi);
		RFinder RF(PEpot, E, I3, Lz);
		double Rout=R-1e-8, Rin=.9*Rout, outer=RF.value(.95*Rout);
		while(RF.value(Rin)*outer > 0 && Rin>1e-5) Rin *= .9;
		if(Rin>1.2e-5) Rin = math::findRoot(RF,Rin,Rout,tol);
		getCurve(PEpot,np,E,Lz,I3,Rin,Rout,Rs,vRs);
		{
			bool onlyTM=false;//true;//
			double vR=0, vphi=Lz/R, vz=sqrt(2*(E-Phi)-pow_2(vphi)); 
			coord::PosVelCyl Rzp(R,0,0,vR,vz,vphi);
			actions::Frequencies omegas;
			actions::ActionAngles aa(actions::actionAnglesAxisymStaeckel(PEpot, Rzp, &omegas));
			printf("J: (%f %f) omegas %f %f %f %f\n",aa.Jr,aa.Jz,
			       omegas.Omegar,omegas.Omegaz,omegas.Omegaphi,omegas.Omegaz/omegas.Omegaphi);
            actions::Angles t(0.2,0.4,1.7);
			actions::Actions J(aa.Jr,aa.Jz,aa.Jphi);
			actions::Torus T;
			if(onlyTM) T=tg.fitBaseTorus(J);
			else T=tg.fitTorus(J);
			actions::ToyPotType ty=T.ptrTM->getToyMapType();
			switch(ty){
				case actions::ToyPotType::None : printf("No TM\n"); break;
				case actions::ToyPotType::Is   : printf("Is TM\n"); break;
				case actions::ToyPotType::HO   : printf("HO TM\n"); break;
			}
            actions::Angles t1(0,0,0);
            coord::PosVelCyl xv0=coord::toPosVelCyl(T.from_true(t1));
            std::vector<double> R,pR,R1,pR1;
            orbit::makeSoS(xv0,*pot,R,pR,100);
            double Rmin,Rmax,Vmax;
            T.zSoS(R1,pR1,100,Rmin,Rmax,Vmax);
            std::ofstream File2("SoSreal_"+std::to_string(j)+".txt");
            for(int i=0;i<Rs.size();i++){
                File2<<Rs[i]<<" "<<vRs[i]<<"\n";
            }
            File2.close();
            std::ofstream Filen("SoSt_"+std::to_string(j)+".txt");
            for(int i=0;i<R1.size();i++){
                Filen<<R1[i]<<" "<<pR1[i]<<"\n";
            }
            Filen.close();
		}
		Rs.clear(); vRs.clear();
		double diff=.054*(Rpmax-Rpmin);
		Rpmin+=.1*diff;
		Rpmax-=diff;
	}
    return 0;
}