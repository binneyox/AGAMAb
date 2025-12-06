/* Code to compute Fig 5 of Binney et al (2026). Orbit data are written
 * on testAFTG.dat. The code is a little slow because a torus is fitted
 * to the J reovered from the Fudge twice to avoid adding functionality
 * to AGAMAb
 */

#include "actions_base.h"
#include "actions_newtorus.h"
#include "actions_newisochrone.h"
#include "potential_analytic.h"
#include "potential_utils.h"
#include "potential_factory.h"
#include "potential_interpolators.h"
#include "orbit.h"

using namespace actions;

int main(int narg,char** args)
{
	FILE* ofile;
	fopen_s(&ofile,"testAFTG.dat","w");
	double Lz=0.5,Ttime=500,dt = .3e-3*Ttime;
	if(narg==2){
		sscanf(args[1],"%lf",&Lz);		
	}
	potential::PtrPotential pot=potential::createPotential(utils::KeyValueMap(
		"type=spheroid gamma=1 beta=3 outercutoffradius=10 axisratioz=0.5"));
	double Phi0=pot->value(coord::PosCyl(0,0,0));
	printf("Potential created\n");
	ActionFinderAxisymFudge afF(pot);
	TorusGenerator TG(*pot,1e-6);
	ActionFinderTG afT(pot,TG,&afF);
	double E  = -0.3;
	double Lc = potential::L_circ(*pot, E);
	Lz *= Lc;
// Compute ICs from shell orbit
	double Rsh = pot->getDelta(E,Lz/L_circ(*pot,E),1/Phi0);
	double vm=sqrt(2*(E-pot->value(coord::PosCyl(Rsh,0,0)))-pow_2(Lz/Rsh));
	const int numang=9, minang=3;
/* Launch at fixed speed in merid plane at increasing angle to r */ 
	for(int a=minang; a<minang+1; a+=2) {
		double ang = M_PI/2 * (double)a/(double)numang;
		coord::PosVelCyl ic(Rsh, 0, 0, vm*cos(ang),
				    vm*sin(ang), Lz/Rsh);
		std::vector<std::pair<coord::PosVelCyl, double> > traj = orbit::integrateTraj(
			ic, Ttime, dt, *pot);
		double Rmin=1e6,Rmax=0,zmax=0;
		/* Compute orbit from ics */
		std::vector<double> RKRs, RKzs;
		for(int i=0; i<traj.size(); i++){
			RKRs.push_back(traj[i].first.R); RKzs.push_back(traj[i].first.z);
			Rmax=fmax(Rmax,RKRs.back()); zmax=fmax(zmax,RKzs.back());
			Rmin=fmin(Rmin,RKRs.back());
		}
		fprintf(ofile,"%zd (R,z) points on orbit from R-K\n",RKRs.size());
		for(int i=0; i<RKRs.size(); i++) fprintf(ofile,"%8.4f\t%8.4f\n",
			RKRs[i],RKzs[i]);
		std::vector<double> Rs, zs;
		/* Use actionFinder to find torus through ICs */
		Frequencies Freqs;
		ActionAngles aaF(afF.actionAngles(ic,&Freqs));
		Actions JF(aaF);
		Torus TF(TG.fitTorus(JF));//Fudged torus
		printf("fitted Fudged torus energy %f strength %f \n",TF.E,TF.GF.strength());
		Torus T;
		ActionAngles aa(afT.actionAnglesTorus(ic, T));
		printf("Vz %f, Fudged Actions (%f %f) freqs (%f %f)\n",ic.vz,JF.Jr,JF.Jz,Freqs.Omegar,Freqs.Omegaphi);
		printf("Vz %f,   TGed Actions (%f %f) freqs (%f %f)\n",ic.vz,T.J.Jr,T.J.Jz,T.freqs.Omegar,T.freqs.Omegaphi);
		printf("Es %f %f\n",TF.E,T.E);
		//Get time sequence from torus
		std::vector<std::pair<coord::PosVelCyl, double> > Ttraj = T.orbit(aa,dt,Ttime);
		for(int i=0; i<Ttraj.size(); i++){
			Rs.push_back(Ttraj[i].first.R); zs.push_back(Ttraj[i].first.z);
		}
		fprintf(ofile,"%zd (R,z) points of orbit from TG torus\n",Rs.size());
		for(int i=0; i<Rs.size(); i++)
			fprintf(ofile,"%8.4f\t%8.4f\n",Rs[i],zs[i]);
		Ttraj.clear();	Ttraj = TF.orbit(aa,dt,Ttime);
		Rs.clear(); zs.clear();
		for(int i=0; i<Ttraj.size(); i++){//plot Fudged torus
			Rs.push_back(Ttraj[i].first.R); zs.push_back(Ttraj[i].first.z);
		}
		fprintf(ofile,"%zd (R,z) points of orbit from Fudged torus\n",Rs.size());
		for(int i=0; i<Rs.size(); i++)
			fprintf(ofile,"%8.4f\t%8.4f\n",Rs[i],zs[i]);
	}
}
