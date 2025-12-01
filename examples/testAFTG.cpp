#include "actions_newtorus.h"
#include "actions_newisochrone.h"
#include "potential_analytic.h"
#include "potential_utils.h"
#include "potential_bars.h"
#include "potential_factory.h"
#include "orbit.h"

int main(void){
	std::vector<double> Rs,zs,vRs,vzs,thetas,pthetas;
	double Rmin=7, Rmax=30, Zmax=15, Vmax=.062;
//	mgo::plt pl;
	potential::PtrPotential pot = potential::createPotential(utils::KeyValueMap(
		"type=spheroid, gamma=1, beta=4, scaleradius=1, q=0.6"));
	actions::TorusGenerator TG(*pot,1e-4);
	actions::Actions J(.3,.4,1);
	actions::Torus T(TG.fitTorus(J));
	coord::PosVelCyl Rzv(coord::toPosVelCyl(T.from_true(actions::Angles(.4,.4,.4))));
	double duration=80*M_PI/T.freqs.Omegar, dt=duration/(double)2000;
	std::vector<std::pair<coord::PosVelCyl,double> >
			traj(orbit::integrateTraj(Rzv,duration,dt,*pot));
	for(int i=0;i<traj.size();i++){
		Rs.push_back(traj[i].first.R); zs.push_back(traj[i].first.z);
	}

	actions::ActionFinderAxisymFudge AF(pot);
	actions::ActionFinderTG AFTG(pot,AF,TG);
	for(int i=0; i<traj.size(); i+=400){
		actions::ActionAngles aa0(AFTG.actionAngles(traj[i].first));
		actions::Actions JSt(AF.actions(traj[i].first));
		printf("%f %f %f %f\n",aa0.Jr,aa0.Jz,JSt.Jr,JSt.Jz);
	}
//	printf("aa0: Jr %f Jz %f Jphi %f thetar %f thetaz %f\n",aa0.Jr,aa0.Jz,aa0.Jphi,aa0.thetar,aa0.thetaz);
}