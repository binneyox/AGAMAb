#include "math_random.h"
#include "potential_analytic.h"
#include "potential_utils.h"
#include "actions_newtorus.h"
#include "potential_factory.h"
#include "GD1.h"
#include "/u/sm/mgo.h"

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
/*
	math::Matrix<double> V(SVD.V());
	printf("Checking diagonalisation\n");
	math::Matrix<double> B(3,3),C(3,3);
	for(int i=0; i<3; i++){
		for(int j=0; j<3; j++){
			B(i,j)=0;
			for(int k=0; k<3; k++) B(i,j)+=M(i,k)*V(k,j);
		}
	}
	for(int i=0; i<3; i++){
		for(int j=0; j<3; j++){
			C(i,j)=0;
			for(int k=0; k<3; k++) C(i,j)+=U(k,i)*B(k,j);
			printf("%f ",C(i,j));
		}
		printf("\n");
	}
	for(int i=0;i<3;i++){
		for(int j=0;j<3;j++) printf("%8.3f ",M(i,j));
		for(int j=0;j<3;j++) printf("%8.3f ",U(i,j));
		for(int j=0;j<3;j++) printf("%8.3f ",V(i,j));
		printf("\n");
	}
*/
	return U;
}

int main(void){
	double h=intUnits.from_Kpc*intUnits.from_kms;
	obs::solarShifter sun(intUnits);

	double Vc=220*intUnits.from_kms; //	if(Vc>0) return 0;
	double Rc=0.1*intUnits.from_Kpc;
	double ZtoX=0.9;
	potential::PtrPotential potL(new potential::Logarithmic(Vc,Rc,1,ZtoX,1000));
	potential::PtrPotential pot=potential::wrapPotential(potL);
	TorusGenerator TG(*pot,.5e-4);
	printf("done\n");
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
	printf("OmPp %f %f %f\n",OmPp.Omegar,OmPp.Omegaz,OmPp.Omegaphi);
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
				if(i%Np==0 && j==0){
					printf("%f %f %f %f %f\n",ells[ip],bs[ip],Xs[ip],Ys[ip],Zs[ip]);
					printf("%f %f %f %f\n",Js[ip].Jr,Js[ip].Jphi,Omegar[ip],Omegaphi[ip]);
				}
			}
		}
	}
	std::time_t end_t = std::time(NULL);
	std::cout << end_t-start_t <<" seconds  to create model\n";
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
		
	mgo::plt pl;
	pl.new_plot(0,1.2,0,5,"|\\gD\\gq\\d\\g||","t\\\\ddrft\\\\u (Gyr)");
	pl.points(33.05,thetapar,tdrft);
	pl.grend();
	
	pl.sq_plot(8,14,9,12,"X/kpc","Z/kpc");
	pl.points(33.05,Xs,Zs);
	pl.setcolour("red"); pl.connect(xp,zp);
	pl.grend();

	for(int i=0; i<2*Np; i++){
		Xs[i]=Js[i].Jr/(220*h); Ys[i]=Js[i].Jz/(220*h); Zs[i]=Js[i].Jphi/(220*h);
	}
	pl.multipane(2, -2);
	pl.pane(0,1);
	{
		pl.new_plot(JPp.Jr-dJp.Jr,JPp.Jr+dJp.Jr,JPp.Jz-dJp.Jz,JPp.Jz+dJp.Jz,
			    "J\\dr (220 kpc s\\u-\\u1)","J\\dz",.8);
		pl.setcolour("red");
		pl.relocate(TP.J.Jr/(220*h),TP.J.Jz/(220*h));
		pl.point(80.3); pl.setcolour("black");
	}
	pl.points(33.05,Xs,Ys);
	pl.pane(0,0);
	{
		pl.new_plot(JPp.Jr-dJp.Jr,JPp.Jr+dJp.Jr,JPp.Jphi-dJp.Jphi,JPp.Jphi+dJp.Jphi,
			    "J\\dr (220 kpc km s\\u-\\u1)","J\\d\\gf",.8);
		pl.setcolour("red");
		pl.relocate(TP.J.Jr/(220*h),TP.J.Jphi/(220*h));
		pl.point(80.3); pl.setcolour("black");
	}
	pl.points(33.05,Xs,Zs);
	pl.pane(1,1);
	{
		pl.new_plot(OmPp.Omegar-dOmp.Omegar,OmPp.Omegar+dOmp.Omegar,
			    OmPp.Omegaz-dOmp.Omegaz,OmPp.Omegaz+dOmp.Omegaz,
			    "\\gW\\dr (Gyr\\u-\\u1)","\\gW\\dz",.8);
		pl.setcolour("red");
		pl.relocate(TP.freqs.Omegar/intUnits.to_Gyr,TP.freqs.Omegaz/intUnits.to_Gyr);
		pl.point(80.3); pl.setcolour("black");
	}
	pl.points(33.05,Omegar,Omegaz);
	pl.pane(1,0);
	{
		pl.new_plot(OmPp.Omegar-dOmp.Omegar,OmPp.Omegar+dOmp.Omegar,
			    OmPp.Omegaphi-dOmp.Omegaphi,OmPp.Omegaphi+dOmp.Omegaphi,
			    "\\gW\\dr (Gyr\\u-\\u1)","\\gW\\d\\gf",.8);
		pl.setcolour("red");
		pl.relocate(TP.freqs.Omegar/intUnits.to_Gyr,TP.freqs.Omegaphi/intUnits.to_Gyr);
		pl.point(80.3); pl.setcolour("black");
	}
	pl.points(33.05,Omegar,Omegaphi);
	pl.grend();
	pl.multipane(1,2);
	pl.new_plot(0,55,22,27.5,"\\sl (deg)","b");
	pl.points(33.05,ells,bs);
	pl.setcolour("red"); pl.connect(ellp,bp);
	pl.grend();
}	