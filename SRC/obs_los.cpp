#include "obs.h"
#include "math_ODE.h"


namespace obs{

EXP SunLos::SunLos(const PosSky _pos, const solarShifter _sun, dust::dustModel* _dm)//angles in degrees// added for sKpc
: pos(_pos), cosl(cos(pos.l*_sun.torad)), cosb(cos(pos.b*_sun.torad)),
sinl(sin(pos.l*_sun.torad)), sinb(sin(pos.b*_sun.torad)),
sun(_sun), BaseLos(true, 0, _sun.from_Kpc, _dm)
{
	if(_dm) tab_extinct();
}
EXP SunLos::SunLos(const double l,const double b, const solarShifter _sun,
		   dust::dustModel* _dm): //angles in degrees
    pos(l,b), cosl(cos(l*_sun.torad)), sinl(sin(l*_sun.torad)),
    cosb(cos(b*_sun.torad)), sinb(sin(b*_sun.torad)),
    sun(_sun), BaseLos(true, 0, _sun.from_Kpc, _dm)
{
	if(_dm) tab_extinct();
}

EXP SunLos::SunLos(const coord::PosCar xyz, const solarShifter _sun, dust::dustModel* _dm)
: BaseLos(true, 0, _sun.from_Kpc, _dm), sun(_sun)
{
	double xhel = xyz.x-sun.xyz().x, yhel=xyz.y-sun.xyz().y, zhel=xyz.z-sun.xyz().z;//position wrt Sun
	double s = sqrt(xhel*xhel+yhel*yhel+zhel*zhel);
	pos.b = asin(zhel/s); pos.l = atan2(yhel,xhel);
	cosl = cos(pos.l); sinl = sin(pos.l);
	cosb = cos(pos.b); sinb = sin(pos.b);
	pos.b/=sun.torad; pos.l/=sun.torad;
	if(_dm) tab_extinct();
}

EXP coord::PosCyl SunLos::deepest() const{//locate clostest distance r_p of los to GC
	double r_p = sun.xyz().x*cosb*cosl + sun.xyz().y*cosb*sinl
		     + sun.xyz().z*sinb;
	return Rzphi(r_p);
}

//helper class for tabulating extinctions: we integrate dN/ds = rho
//along LOS. Distances in intUnits
class A_syst: public math::IOdeSystem{
	private:
		const obs::BaseLos* Los;
		const double start;
		const dust::BaseDustModel* dm;
	public:
		A_syst(const obs::BaseLos* _los, double _start, dust::BaseDustModel* _dm)
				: Los(_los), start(_start), dm(_dm){}
		virtual void eval(const double s,const double A[], double dAds[]) const{
			dAds[0] = dm->dens(Los->Rzphi(start+s));
		}
		virtual unsigned int size() const {return 1;}
};

EXP void BaseLos::tab_extinct(void){
	std::vector<double> dists, Ns;
	const int ns=3;
	A_syst syst(this, start(), dm);
	math::OdeSolverDOP853 solver(syst);
	double s=start(), rho=fmax(0,dm->dens(Rzphi(s))), N=0, rhoPeak=0;
	int iz=0;
	solver.init(&N);
	dists.push_back(s); Ns.push_back(N);
	while(iz<200 && rho >= 0.01*rhoPeak){
		solver.doStep();
		double s_last=s, ds=(start()+solver.getTime()-s_last)/(double)ns;
		for(int i=0; i<ns; i++){// store ns points over step
			s+=ds;
			dists.push_back(s); Ns.push_back(solver.getSol(s-start(),0));
		}
		rho=dm->dens(Rzphi(s)); rhoPeak= rho>rhoPeak? rho : rhoPeak;
		iz++;
	}
	extinct = math::CubicSpline(dists,Ns);
}

EXP extLos::extLos(units::InternalUnits _intUnits, double _xs, double _ys, double _incl, double _s0,
		   dust::BaseDustModel* _dm)
:  intUnits(_intUnits), xs(_xs), ys(_ys), incl(_incl), 
BaseLos(false, _s0, _intUnits.from_Kpc, _dm),
cosi(cos(_incl*M_PI/180.)), sini(sin(_incl*M_PI/180.))
{
	if(_dm) tab_extinct();
}

} //namespace obs