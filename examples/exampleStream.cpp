/* Compute data Figs 8 & 9 of Binney et al (2026). Data written to
 * Stream.dat
 */

#include "math_random.h"
#include "potential_analytic.h"
#include "potential_utils.h"
#include "actions_newtorus.h"
#include "potential_factory.h"
#include "obs.h"

// define internal unit system - arbitrary numbers here! the result should not depend on their choice
const units::InternalUnits intUnits(2.7183*units::Kpc, 3.1416*units::Myr);
const double torad=M_PI/180;

using namespace actions;

math::Matrix<double> compute_dOmegadJ(const Actions JP,const TorusGenerator& TG,
				      std::vector<double>& lambdas){
	double dJ=.001;
	Actions J(JP.Jr+.5*dJ, JP.Jz, JP.Jphi);
	Torus T(TG.fitTorus(J));
	math::Matrix<double> M(3,3);
	M(0,0)=T.freqs.Omegar; M(1,0)=T.freqs.Omegaz; M(2,0)=T.freqs.Omegaphi;
	J.Jr -= dJ;
	T=TG.fitTorus(J);
	M(0,0)-=T.freqs.Omegar; M(1,0)-=T.freqs.Omegaz; M(2,0)-=T.freqs.Omegaphi;
	M(0,0)/=dJ; M(1,0)/=dJ; M(2,0)/=dJ;
	J.Jr += .5*dJ;
	J.Jz += .5*dJ;
	T=TG.fitTorus(J);
	M(0,1)=T.freqs.Omegar; M(1,1)=T.freqs.Omegaz; M(2,1)=T.freqs.Omegaphi;
	J.Jz -= dJ;
	T=TG.fitTorus(J);
	M(0,1)-=T.freqs.Omegar; M(1,1)-=T.freqs.Omegaz; M(2,1)-=T.freqs.Omegaphi;
	M(0,1)/=dJ; M(1,1)/=dJ; M(2,1)/=dJ;
	J.Jz += .5*dJ;
	J.Jphi += .5*dJ;
	T=TG.fitTorus(J);
	M(0,2)=T.freqs.Omegar; M(1,2)=T.freqs.Omegaz; M(2,2)=T.freqs.Omegaphi;
	J.Jphi -= dJ;
	T=TG.fitTorus(J);
	M(0,2)-=T.freqs.Omegar; M(1,2)-=T.freqs.Omegaz; M(2,2)-=T.freqs.Omegaphi;
	M(0,2)/=dJ; M(1,2)/=dJ; M(2,2)/=dJ;
	J.Jphi += .5*dJ;
//Make M symmetric (which it nearly is)
	double a;
	a=.5*(M(0,1)+M(1,0)); M(0,1)=a; M(1,0)=a;
	a=.5*(M(0,2)+M(2,0)); M(0,2)=a; M(2,0)=a;
	a=.5*(M(2,1)+M(1,2)); M(2,1)=a; M(1,2)=a;
	math::SVDecomp SVD(M);
	lambdas=SVD.S();
	math::Matrix<double> U(SVD.U());
	return U;
}

int main(void){
	FILE* ofile;
	fopen_s(&ofile,"Stream.dat","w");
	double h=intUnits.from_Kpc*intUnits.from_kms;
	obs::solarShifter sun(intUnits);
	double Vc=220*intUnits.from_kms; //	if(Vc>0) return 0;
	double Rc=0.1*intUnits.from_Kpc;
	double ZtoX=0.9;
	potential::PtrPotential potL(new potential::Logarithmic(Vc,Rc,1,ZtoX,1000));
	potential::PtrPotential pot=potential::wrapPotential(potL);
	printf("Potential created\n");
	TorusGenerator TG(*pot,.5e-4);
	std::time_t start_t = std::time(NULL);

//Progenitor actions from Bovy 2014
	Actions JP(288.5*h,897.6*h,3173.7*h);
	Torus TP(TG.fitTorus(JP));
	double TrP=2*M_PI/TP.freqs.Omegar;
	coord::PosMomCyl RZperi(TP.from_true(Angles(0,0,0))),
		RZapo(TP.from_true(Angles(M_PI,0,0))),
		RZzmax(TP.from_true(Angles(M_PI,.5*M_PI,0)));
//compute dOmega/dJ
	std::vector<double> lambdas;
	math::Matrix<double> U=compute_dOmegadJ(JP,TG,lambdas);
	double e1[3]={U(0,0), U(1,0), U(2,0)};
	double OmegaparP = e1[0]*TP.freqs.Omegar + e1[1]*TP.freqs.Omegaz
			   + e1[2]*TP.freqs.Omegaphi;
//dispersion from Bovy 2014 cluster dispersion 
	double sigCluster = 0.365*intUnits.from_kms;
	double DJr = sigCluster*(RZapo.R-RZperi.R)/M_PI;
	double DJz = sigCluster*2*RZzmax.z/M_PI;
	double DJphi = sigCluster*RZperi.R;
	double H=220*h;
	Actions JPp = JP/H, dJp(6*DJr/H,6*DJz/H,8*DJphi/H);
	Frequencies OmPp = TP.freqs/intUnits.to_Gyr, dOmp(.2,.2,.2);
	printf("Progenitor freqs (Gyr^{-1}) %f %f %f\n",OmPp.Omegar,OmPp.Omegaz,OmPp.Omegaphi);
//Dispersion of ts about tperi
	double DT = 30*intUnits.from_Myr;
	Angles theta0(0,M_PI,0);
	int Np=1000, Nperi=10;
	double t0=((double)Nperi+.2)*TrP, total=(Nperi)*TrP;
	std::vector<double> Xs(2*Np), Ys(2*Np), Zs(2*Np), ells(2*Np), bs(2*Np), tdrft(2*Np);
	std::vector<double> ss(2*Np), Vlos(2*Np);
	std::vector<double> thetapar(2*Np), Omegar(2*Np),Omegaz(2*Np),Omegaphi(2*Np);
	std::vector<Actions> Js(2*Np);
// set up grid width 3*DJr etc around JP
	for(int krep = -1; krep<2; krep+=2){//leading then trailing stream 
		int nx=5, ny=5, nz=5;
		std::vector<Torus> Ts(nx*ny*nz);
		std::vector<double> xs(nx), ys(ny), zs(nz);
		Actions Joff(krep*.7*DJr, krep*1.5*DJz, krep*5*DJphi);//leading stream (greater |Jphi|)
		printf("Computing Tgrid..");
#pragma omp parallel for
		for(int i=0; i<nx; i++){
			printf("%d ",i);
			xs[i]=JP.Jr+Joff.Jr+(i-nx/2)*DJr;
			for(int j=0; j<ny; j++){
				ys[j]=JP.Jz+Joff.Jz+(j-ny/2)*DJz;
				for(int k=0; k<nz; k++){
					zs[k]=JP.Jphi+Joff.Jphi+(k-nz/2)*DJphi;
					Actions J(xs[i],ys[j],zs[k]);
					Ts[(i*nx+j)*ny+k]=TG.fitTorus(J);
				}
			}
		}
		printf("..done\n");
		TorusGrid3 Tgrd(xs,ys,zs,Ts,TG);
// Sample Gaussian  widths DJr etc around JP
		for(int i=0; i<Np; i+=2){
			double r[8];
			for(int k=0; k<7; k+=2) math::getNormalRandomNumbers(r[k],r[k+1]);
			for(int j=0; j<2; j++){
				int ip=(krep+1)/2*Np+i+j;
				double dJphi=DJphi*r[3*j+2], fac=.25*(3.5+krep*dJphi/DJphi);
				Js[ip]=Actions(JP.Jr+Joff.Jr+fac*DJr*r[3*j],
					JP.Jz+Joff.Jz+fac*DJz*r[3*j+1],
					JP.Jphi+Joff.Jphi+dJphi);
				Torus T(Tgrd.T(Js[ip]));
// Advance each group by n*TrP and compute sky coords;
				int n=total*math::random()/TrP; //index of release peri
				double ts=n*TrP+DT*r[3*j+3];	//release time
				tdrft[ip]=(t0-ts)*intUnits.to_Gyr;
				Omegar[ip]=T.freqs.Omegar/intUnits.to_Gyr;
				Omegaz[ip]=T.freqs.Omegaz/intUnits.to_Gyr;
				double Omegapar = e1[0]*T.freqs.Omegar
					+ e1[1]*T.freqs.Omegaz
					+ e1[2]*T.freqs.Omegaphi;
				Omegaphi[ip]=T.freqs.Omegaphi/intUnits.to_Gyr;
//Now the current appearance
				double r1,r2,r3;
				math::getNormalRandomNumbers(r1,r2);
				math::getNormalRandomNumbers(r1,r3);
				Angles dtheta(r1*M_PI/1000,r2*M_PI/1000,r3*M_PI/1000);
				double dthetapar=e1[0]*dtheta.thetar+e1[1]*dtheta.thetaz
					+e1[2]*dtheta.thetaphi;
				Angles theta(theta0+dtheta+(Angles)(TP.freqs*ts+T.freqs*(t0-ts)));
				thetapar[ip]=fabs(dthetapar + (Omegapar-OmegaparP)*(t0-ts));
				coord::PosVelCyl pv(coord::toPosVelCyl(T.from_true(theta)));
				Xs[ip]=coord::toPosCar(pv).x*intUnits.to_Kpc;
				Ys[ip]=coord::toPosCar(pv).y*intUnits.to_Kpc;
				Zs[ip]=pv.z*intUnits.to_Kpc;
				obs::PosVelSky lbmu(sun.toSky(pv, ss[ip], Vlos[ip]));
				ells[ip]=lbmu.pos.l; bs[ip]=lbmu.pos.b;
			}
		}
	}
	std::time_t end_t = std::time(NULL);
	std::cout << end_t-start_t <<" seconds  to create grids\n";
	int np=200;
	std::vector<double> xp,zp,ellp,bp;
	for(int i=0; i<np; i++){
		double t=t0*(1-.01+.02*(double)i/(double)(np-1));
		Angles theta(theta0+Angles(TP.freqs*t));
		coord::PosVelCyl pv(coord::toPosVelCyl(TP.from_true(theta)));
		xp.push_back(coord::toPosCar(pv).x*intUnits.to_Kpc);
		zp.push_back(pv.z*intUnits.to_Kpc);
		double s;
		obs::PosSky lb(sun.toSky(pv,s));
		ellp.push_back(lb.l); bp.push_back(lb.b);
	}
	std::cout << std::time(NULL) - end_t << " secs to create model\n", 
	fprintf(ofile,"%zd X,Z,b,ell coords of stars\n",Xs.size());
	for(int i=0; i<Xs.size(); i++)
		fprintf(ofile,"%8.4f\t%8.4f\t%8.4f\t%8.4f\n",
			  Xs[i],Zs[i],ells[i],bs[i]);
	
	for(int i=0; i<2*Np; i++){
		Xs[i]=Js[i].Jr/(220*h); Ys[i]=Js[i].Jz/(220*h); Zs[i]=Js[i].Jphi/(220*h);
	}
	fprintf(ofile,"%8.4f\t%8.4f\t%8.4f Jr Jz Jphi of progenitor (units 220 kpc kms)\n",
		  TP.J.Jr/(220*h),TP.J.Jz/(220*h),TP.J.Jphi/(220*h));
	fprintf(ofile,"%zd actions Jr, Jz, Jphi of stars  (units 220 kpc kms)\n",Xs.size());
	for(int i=0; i<Xs.size(); i++)
		fprintf(ofile,"%8.4f\t%8.4f\t%8.4f\n",Xs[i],Ys[i],Zs[i]);
	fprintf(ofile,"%8.4f %8.4f %8.4f progenitor frequencies (units Gyr^{-1})\n",
		  TP.freqs.Omegar/intUnits.to_Gyr,
		  TP.freqs.Omegaz/intUnits.to_Gyr,
		  TP.freqs.Omegaphi/intUnits.to_Gyr);
	fprintf(ofile,"%zd frequencies of stars\n",bs.size());
	for(int i=0; i<bs.size(); i++)
		fprintf(ofile,"%8.4f\t%8.4f\t%8.4f\n",
		       Omegar[i],Omegaz[i],Omegaphi[i]);
}