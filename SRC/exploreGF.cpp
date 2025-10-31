#include "actions_newtorus.h"
#include "actions_newisochrone.h"
#include "potential_analytic.h"
#include "potential_utils.h"
#include "potential_bars.h"
#include "potential_factory.h"
#include "orbit.h"
#include "/u/sm/mgo.h"

void exploreGF(actions::GenFnc& GF, actions::Actions& J, mgo::plt& pl){
	double rmax=0,zmax=0;
	std::pair<int,int> maxIs(GF.maxIndices());
	for(int i=0;i<2*maxIs.first;i++){
		double thetar=M_PI*i/(2*(double)maxIs.first);
		for(int j=0; j<2*maxIs.second; j++){
			double thetaz=M_PI*i/(2*(double)maxIs.second);
			actions::Angles thetaT(thetar,thetaz,0);
			actions::Actions JT(GF.toyJ(J,thetaT));
			rmax=fmax(rmax,fabs(J.Jr-JT.Jr));
			zmax=fmax(zmax,fabs(J.Jz-JT.Jz));
		}
	}
	printf("Biggest fractional J changes: %f %f\n",rmax/J.Jr,zmax/J.Jz);

	pl.new_plot(0,(double)maxIs.first+.5,0,1.,"n\\dr","S\\\\dn+1/\\\\uS\\dn");
	std::vector<double> xs,Ss; double Slast;
	for(int i=1; i<=maxIs.first; i++){
		double S = GF.giveValue(actions::GenFncIndex(i,2,0));
		if(!std::isnan(S) && i>1){
			xs.push_back(i); Ss.push_back(S/Slast);
		}
		Slast=S;
	}
	pl.points(43.2,xs,Ss);
	printf("Ratios S(i+1,2)/S(i,2))\n");
	for(int i=0; i<Ss.size(); i++) printf("%f ",Ss[i]);
	xs.clear(); Ss.clear();
	for(int i=1; i<=maxIs.first; i++){
		double S = GF.giveValue(actions::GenFncIndex(i,0,0));
		if(!std::isnan(S) && i>1){
			xs.push_back(i); Ss.push_back(S/Slast);
		}
		Slast = S;
	}
	pl.setcolour("red");
	pl.points(33.2,xs,Ss);
	pl.grend();
	printf("Ratios S(i+1,0)/S(i,0))\n");
	for(int i=0; i<Ss.size(); i++) printf("%f ",Ss[i]);
}

int main(void){
	FILE* ifile; fopen_s(&ifile,fname,"r");
	actions::GenFncFitSeries GFFS;
	std::vector<double> params;
	GFFS.read(ifile,params);
	actions::GenFncDerivs derivs(GFFS.numTerms());
	actions::GenFnc GF(GFFS.indices,params,derivs);
	mgo::plt pl;
	actions::Actions J(0.75,0.5,2);
	exploreGF(GF J, pl);
}
	