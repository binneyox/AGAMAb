/** \file    obs_los.h
    \brief   Defines lines of sight
    \author  James Binney
    \date    2022
*/
#pragma once

namespace obs {

class EXP BaseLos{
	private:
		math::CubicSpline extinct;
		dust::BaseDustModel* dm;
	public:
		bool semi_inf;
		const double s0, from_Kpc;
		//BaseLos(){};
		virtual ~BaseLos(){};
		BaseLos(const bool _semi_inf, const double _s0, const double _from_Kpc,
			  dust::BaseDustModel* _dm=NULL):
		    semi_inf(_semi_inf), s0(_s0), from_Kpc(_from_Kpc), dm(_dm) {
		}
		virtual coord::PosCyl deepest() const{
			printf("Erroneous call to deepest()");
			return coord::PosCyl(0,0,0);
		};
		virtual coord::PosCar xyz(const double s) const{
			printf("Error: xyz call in BaseLos\n");
			return coord::PosCar(0,0,0);
		}
		coord::PosCyl Rzphi(const double s) const{// s in internal units
			return coord::PosCyl(coord::toPosCyl(xyz(s)));			
		}
		virtual double s(const coord::PosCar& xyz) const{// s in internal units
			printf("Error: s call in BaseLos\n");
			return 0;
		}
		double s(const coord::PosCyl& Rz) const{
			return s(coord::toPosCar(Rz));
		}
		virtual std::pair<double,double> sVlos(const coord::PosVelCyl& xv) const{// s Vlos int units
			return std::make_pair(0,0);
		}
		virtual double sMod_H(const coord::PosVelCyl) const=0;//m-M including A
		virtual double reduc_H(const coord::PosVelCyl) const=0;//m-M from A only
		double A_V(const double s) const{
			if(!dm) return 0;
			else{
				double A,dA;
				extinct.evalDeriv(s,&A,&dA);
				return A;
			}
		}
		double A_B(const double s) const{
			return 1.324 * A_V(s);
		}
		double A_R(const double s) const{
			return 0.748 * A_V(s);
		}
		double A_H(const double s) const{
			return 0.175 * A_V(s);
		}
		double A_K(const double s) const{
			return 0.112 * A_V(s);
		}
		virtual double start() const{
			printf("Erroneous call to start\n");
			return 0;
		}
		virtual void tab_extinct(void);
};

/*
 *class for lines of sight through external galaxies
 *(xs,ys) coords (in kpc) wrt apparent maj/min axies
 *incl inclination (deg), s0 distance to galaxy (kpc)
 */
class EXP extLos: public obs::BaseLos {
	public: 
		units::InternalUnits intUnits;
		const double xs, ys, incl, cosi, sini;//xs ys in intUnits
		extLos(units::InternalUnits _intUnits,
		       double _xs, double _ys, double _incl, double s0,
		       dust::BaseDustModel* _dm=NULL);
		virtual coord::PosCyl deepest() const{// Point of closest approach to GC
			return Rzphi(-ys*cosi/sini);
		}
		// s in internal units is measured from crossing of xz plane
		virtual coord::PosCar xyz(const double s) const{
			double x = xs;
			double y = s*sini;
			double z = s*cosi - ys/sini;
			return coord::PosCar(x, y, z);
		}
		// s in intUnits is measured from crossing of xz plane
		virtual double s(const coord::PosCar& _p) const{
			if(_p.x!=xs || fabs(_p.y-(_p.z*sini - ys)/cosi)>.0001)
				printf("Error p not on LOS in extLos.s\n"); 
			double s = pow_2(_p.y) + pow_2(_p.z-ys/sini);
			s = s>0?  sqrt(s) : 0;
			return _p.y>0? s : -s;
		}
		virtual double s(const coord::PosCyl& p) const{
			return s(coord::toPosCar(p));
		}
		virtual std::pair<double,double> sVlos(const coord::PosVelCyl& p) const{
			coord::PosVelCar xv(coord::toPosVelCar(p));
			return std::make_pair(s(xv), (xv.vy*sini + xv.vz*cosi));
		}
		virtual double start() const{//starting s (intUnits) value for Av integration 
			return sini*cosi!=0? 2*ys/(sini*cosi) : -10*fmax(xs,ys);
//			fabs(sini)>0? (1-ys/intUnits.from_Kpc/sini)/cosi : -30;
		}
		virtual double sMod_H(const coord::PosVelCyl pv) const{
			double s1 = s(pv);
			return 5*log10(100*(s0+s1)*intUnits.to_Kpc) + A_H(s1);
		}
		virtual double reduc_H(const coord::PosVelCyl pv) const{
			return pow(10,-0.4*A_H(s(pv)));
		}
};

//class line of sight from Sun, specified by (l,b) and location of Sun.
//Extinction in BVRHK down it will be computed if a dustModel isspecified;
//otherwise zero extinction
class EXP SunLos: public obs::BaseLos {
	public: 
		PosSky pos;// l, b in degrees
		double cosl, cosb, sinl, sinb;
		const obs::solarShifter sun;
		SunLos(
		       const PosSky _pos, const obs::solarShifter _sun,
		       dust::dustModel* _dm=NULL);
		SunLos(const double l,const double b, const obs::solarShifter _sun,
		       dust::dustModel* _dm=NULL);
		SunLos(const coord::PosCar xyz, const obs::solarShifter _sun,
		       dust::dustModel* _dm=NULL);
		// Point of closest approach to GC
		virtual coord::PosCyl deepest() const;
		// s in internal units
		virtual coord::PosCar xyz(const double s) const{
			double x = sun.xyz().x + s*cosb*cosl;
			double y = sun.xyz().y + s*cosb*sinl;
			double z = sun.xyz().z + s*sinb;
			coord::PosCar Car(x, y, z);
			return Car;
		}
		// s in internal units (this isn't to do with a los)
		virtual double s(const coord::PosCar& p) const{
			double dx = p.x-sun.xyz().x, dy = p.y-sun.xyz().y;
			double dR2 = dx*dx+dy*dy, dz = p.z-sun.xyz().z;
			double s2 = dR2+dz*dz;
			double s= s2>0? sqrt(s2) : 0;
			if(fabs(pos.l*sun.torad - atan2(dy,dx))>1e-5 ||
			   fabs(dz-s*sinb)>1e-5)
				printf("Error p not on LOS in sun.s()");
			return s;
		}
		virtual std::pair<double,double> sVlos(const coord::PosVelCyl& p) const{
			coord::PosVelCar xv(coord::toPosVelCar(p));
			double Vlos=(xv.x*cosl+xv.y*sinl)*cosb+xv.z*sinb;
			return std::make_pair(s(xv), Vlos);
		}
		virtual double sMod_H(const coord::PosVelCyl pv) const{
			double s1=s(coord::toPosCar(pv));
			return 5*log10(100*s1/sun.from_Kpc) + A_H(s1/sun.from_Kpc);
		}
		virtual double reduc_H(const coord::PosVelCyl pv) const{
			double s1=s(coord::toPosCar(pv));
			return pow(10,-0.4*A_H(s1/sun.from_Kpc));
		}
		virtual double start() const{
			return 0;
		}
};

}// namespace obs