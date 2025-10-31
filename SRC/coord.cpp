#include "coord.h"
#include <stdexcept>


namespace coord{

ProlSph::ProlSph(double Delta) : 
    Delta2(Delta*Delta)
{
	if(Delta<=0)
		throw std::invalid_argument("Invalid parameters for Prolate Spheroidal coordinate system");
}

PosT<ProlMod>::PosT(double _rho, double _tau, double _phi, const ProlMod& coordsys) :
    rho(_rho), tau(_tau), phi(_phi), chi(sqrt(pow_2(_rho) + pow_2(coordsys.D))) {}


//--------  angular momentum functions --------//

template<> EXP double Ltotal(const PosVelCar& p) {
    return sqrt(pow_2(p.y*p.vz-p.z*p.vy) + pow_2(p.z*p.vx-p.x*p.vz) + pow_2(p.x*p.vy-p.y*p.vx));
}
template<> EXP double Ltotal(const PosVelCyl& p) {
	return sqrt((pow_2(p.R) + pow_2(p.z)) * pow_2(p.vphi) + pow_2(p.R*p.vz-p.z*p.vR));
}
template<> EXP double Ltotal(const PosMomCyl& p) {
	return p.R>0? sqrt((1 + pow_2(p.z)/pow_2(p.R)) * pow_2(p.pphi) + pow_2(p.R*p.pz-p.z*p.pR))
			: fabs(p.R*p.pz-p.z*p.pR);
}
template<> EXP double Ltotal(const PosVelSph& p) {
	return sqrt(pow_2(p.vtheta) + pow_2(p.vphi)) * p.r;
}
template<> EXP double Ltotal(const PosMomSph& p) {
	return sqrt(pow_2(p.ptheta) + pow_2(p.pphi));
}
template<> EXP double Lz(const PosVelCar& p) { volatile double a = p.x * p.vy, b = p.y * p.vx; return a-b; }
template<> EXP double Lz(const PosVelCyl& p) { return p.R * p.vphi; }
template<> EXP double Lz(const PosVelSph& p) { return p.r * sin(p.theta) * p.vphi; }

// multiply two numbers, replacing {anything including INFINITY} * 0 with 0
inline double mul(double x, double y) { return y==0 ? 0 : x*y; }

//--------  position conversion functions ---------//

template<>
EXP PosCar toPos(const PosCyl& p) {
    double sinphi, cosphi;
    math::sincos(p.phi, sinphi, cosphi);
    return PosCar(mul(p.R, cosphi), mul(p.R, sinphi), p.z);
}
template<>
EXP PosCar toPos(const PosSph& p) {
    double sintheta, costheta, sinphi, cosphi;
    math::sincos(p.theta, sintheta, costheta);
    math::sincos(p.phi, sinphi, cosphi);
    return PosCar(mul(p.r, sintheta*cosphi), mul(p.r, sintheta*sinphi), mul(p.r, costheta));
}
template<>
EXP PosCyl toPos(const PosCar& p) {
    return PosCyl(sqrt(pow_2(p.x) + pow_2(p.y)), p.z, atan2(p.y, p.x));
}
template<>
EXP PosCyl toPos(const PosSph& p) {
    double sintheta, costheta;
    math::sincos(p.theta, sintheta, costheta);
    return PosCyl(mul(p.r, sintheta), mul(p.r, costheta), p.phi);
}
template<>
EXP PosSph toPos(const PosCar& p) {
    return PosSph(sqrt(pow_2(p.x)+pow_2(p.y)+pow_2(p.z)), 
        atan2(sqrt(pow_2(p.x) + pow_2(p.y)), p.z), atan2(p.y, p.x));
}
template<>
EXP PosSph toPos(const PosCyl& p) {
    return PosSph(sqrt(pow_2(p.R) + pow_2(p.z)), atan2(p.R, p.z), p.phi);
}
template<>
EXP PosCyl toPos(const PosProlSph& p) {
	if(fabs(p.nu)>p.coordsys.Delta2 || p.lambda<p.coordsys.Delta2)
		throw std::invalid_argument("Incorrect ProlSph coordinates");
	const double R = sqrt( (p.lambda-p.coordsys.Delta2) * (1 - fabs(p.nu) / p.coordsys.Delta2) );
	const double z = sqrt( p.lambda * fabs(p.nu) / p.coordsys.Delta2) * (p.nu>=0 ? 1 : -1);
	return PosCyl(R, z, p.phi);
}
template<>
EXP PosSph toPos(const PosSphMod& p) {
	double tau2=pow_2(p.tau);
	return PosSph(p.r,(1-tau2)/(1+tau2),p.phi);
}


// declare an instantiation which will be defined later
template<>
EXP PosProlSph toPosDeriv(const PosCyl& from, const ProlSph& cs,
			  PosDerivT<Cyl, ProlSph>* derivs, PosDeriv2T<Cyl, ProlSph>* derivs2);
template<>
EXP PosProlSph toPos(const PosCyl& from, const ProlSph& cs) {
	return toPosDeriv<Cyl,ProlSph>(from, cs, NULL, NULL);
}

// declare an instantiation which will be defined later
template<>
EXP PosCyl toPosDerivR(const PosSphMod& p, PosDerivT<Cyl, SphMod>* derivs,
		  PosDeriv2T<Cyl, SphMod>* derivs2);
template<>
EXP PosCyl toPos(const PosSphMod& p) {
	return PosCyl(toPosDerivR<SphMod, Cyl>(p, NULL, NULL));
/*	double tau2=pow_2(p.tau);
	double sintheta=(1-tau2)/(1+tau2);
	double rt = sintheta<1? sqrt(1-pow_2(sintheta)) : 0;
	double costheta = p.tau*(1+sintheta);
	return PosCyl(p.r*sintheta,p.r*costheta,p.phi);*/
}
EXP UVSph interpUVSph(const double x, const UVSph& C0, const UVSph& C1){
	return UVSph(x*C0.Delta+(1-x)*C1.Delta);
}
template<>
EXP PosCyl toPos(const PosUVSph& p) {
	if(p.u<0)
		throw std::invalid_argument("Incorrect UVSph coordinates");
	const double R = p.coordsys.Delta * sinh(p.u) * sin(p.v);
	const double z = p.coordsys.Delta * cosh(p.u) * cos(p.v);
	return PosCyl(R, z, p.phi);
}
// declare an instantiation which will be defined later
template<>
EXP PosUVSph toPosDeriv(const PosCyl& from, const UVSph& cs,
			  PosDerivT<Cyl, UVSph>* derivs, PosDeriv2T<Cyl, UVSph>* derivs2);
template<>
EXP PosUVSph toPos(const PosCyl& from, const UVSph& cs) {
	return toPosDeriv<Cyl,UVSph>(from, cs, NULL, NULL);
}

//-------- position conversion with derivatives --------//

template<>
EXP PosCyl toPosDeriv(const PosCar& p, PosDerivT<Car, Cyl>* deriv, PosDeriv2T<Car, Cyl>* deriv2) 
{
    const double R2=pow_2(p.x)+pow_2(p.y), R=sqrt(R2);
    if(R==0) {
        // degenerate case, but provide something meaningful nevertheless,
        // assuming that these numbers will be multiplied by 0 anyway
        if(deriv!=NULL)
            deriv->dRdx=deriv->dRdy=deriv->dphidx=deriv->dphidy=1.;
        if(deriv2!=NULL)
            deriv2->d2Rdx2=deriv2->d2Rdy2=deriv2->d2Rdxdy=
            deriv2->d2phidx2=deriv2->d2phidy2=deriv2->d2phidxdy=1.;
        return PosCyl(0, p.z, 0);
    }
    const double cosphi=p.x/R, sinphi=p.y/R;
    if(deriv!=NULL) {
        deriv->dRdx=cosphi;
        deriv->dRdy=sinphi;
        deriv->dphidx=-sinphi/R;
        deriv->dphidy=cosphi/R;
    }
    if(deriv2!=NULL) {
        deriv2->d2Rdx2 =pow_2(sinphi)/R;
        deriv2->d2Rdy2 =pow_2(cosphi)/R;
        deriv2->d2Rdxdy=-sinphi*cosphi/R;
        deriv2->d2phidx2 =2*sinphi*cosphi/R2;
        deriv2->d2phidy2 =-deriv2->d2phidx2;
        deriv2->d2phidxdy=(pow_2(sinphi)-pow_2(cosphi))/R2;
    }
    return PosCyl(R, p.z, atan2(p.y, p.x));
}

template<>
EXP PosCar toPosDeriv(const PosCyl& p, PosDerivT<Cyl, Car>* deriv, PosDeriv2T<Cyl, Car>* deriv2) {
    double sinphi, cosphi;
    math::sincos(p.phi, sinphi, cosphi);
    const double x=mul(p.R, cosphi), y=mul(p.R, sinphi);
    if(deriv!=NULL) {
        deriv->dxdR=cosphi;
        deriv->dydR=sinphi;
        deriv->dxdphi=-y;
        deriv->dydphi= x;
    }
    if(deriv2!=NULL) {
        deriv2->d2xdRdphi=-sinphi;
        deriv2->d2ydRdphi=cosphi;
        deriv2->d2xdphi2=-x;
        deriv2->d2ydphi2=-y;
    }
    return PosCar(x, y, p.z);
}

template<>
EXP PosCar toPosDeriv(const PosSph& p, PosDerivT<Sph, Car>* deriv, PosDeriv2T<Sph, Car>* deriv2) {
    double sintheta, costheta, sinphi, cosphi;
    math::sincos(p.theta, sintheta, costheta);
    math::sincos(p.phi, sinphi, cosphi);
    const double R=mul(p.r, sintheta), x=mul(R, cosphi), y=mul(R, sinphi), z=mul(p.r, costheta);
    if(deriv!=NULL) {
        deriv->dxdr=sintheta*cosphi;
        deriv->dydr=sintheta*sinphi;
        deriv->dzdr=costheta;
        deriv->dxdtheta=z*cosphi;
        deriv->dydtheta=z*sinphi;
        deriv->dzdtheta=-R;
        deriv->dxdphi=-y;
        deriv->dydphi= x;
    }
    if(deriv2!=NULL) {
        deriv2->d2xdrdtheta=costheta*cosphi;
        deriv2->d2ydrdtheta=costheta*sinphi;
        deriv2->d2zdrdtheta=-sintheta;
        deriv2->d2xdrdphi=-sintheta*sinphi;
        deriv2->d2ydrdphi= sintheta*cosphi;
        deriv2->d2xdtheta2=-x;
        deriv2->d2ydtheta2=-y;
        deriv2->d2zdtheta2=-z;
        deriv2->d2xdthetadphi=-z*sinphi;
        deriv2->d2ydthetadphi= z*cosphi;
        deriv2->d2xdphi2=-x;
        deriv2->d2ydphi2=-y;
    }
    return PosCar(x, y, z);
}

template<>
EXP PosSph toPosDeriv(const PosCar& p, PosDerivT<Car, Sph>* deriv, PosDeriv2T<Car, Sph>* deriv2) {
	const double x2=pow_2(p.x), y2=pow_2(p.y), z2=pow_2(p.z);
	const double R2=x2+y2, R=sqrt(R2);
	const double r2=R2+z2, r=sqrt(r2), invr=1/r;
//    if(R==0)
//        throw std::runtime_error("PosDeriv Car=>Sph: R=0, degenerate case!");
	if(deriv!=NULL) {
		deriv->drdx=p.x*invr;
		deriv->drdy=p.y*invr;
		deriv->drdz=p.z*invr;
		const double temp=p.z/(R*r2);
		deriv->dthetadx=p.x*temp;
		deriv->dthetady=p.y*temp;
		deriv->dthetadz=-R/r2;
		deriv->dphidx=-p.y/R2;
		deriv->dphidy=p.x/R2;
	}
	if(deriv2!=NULL) {
		const double invr3=invr/r2;
		deriv2->d2rdx2=(r2-x2)*invr3;
		deriv2->d2rdy2=(r2-y2)*invr3;
		deriv2->d2rdz2=R2*invr3;
		deriv2->d2rdxdy=-p.x*p.y*invr3;
		deriv2->d2rdxdz=-p.x*p.z*invr3;
		deriv2->d2rdydz=-p.y*p.z*invr3;
		const double invr4=1/(r2*r2);
		const double temp=p.z*invr4/(R*R2);
		deriv2->d2thetadx2=(r2*y2-2*R2*x2)*temp;
		deriv2->d2thetady2=(r2*x2-2*R2*y2)*temp;
		deriv2->d2thetadz2=2*R*p.z*invr4;
		deriv2->d2thetadxdy=-p.x*p.y*(r2+2*R2)*temp;
		const double temp2=(R2-z2)*invr4/R;
		deriv2->d2thetadxdz=p.x*temp2;
		deriv2->d2thetadydz=p.y*temp2;
		deriv2->d2phidx2=2*p.x*p.y/pow_2(R2);
		deriv2->d2phidy2=-deriv2->d2phidx2;
		deriv2->d2phidxdy=(y2-x2)/pow_2(R2);
	}
	return PosSph(r, atan2(R, p.z), atan2(p.y, p.x));
}

template<>
EXP PosSph toPosDeriv(const PosCyl& p, PosDerivT<Cyl, Sph>* deriv, PosDeriv2T<Cyl, Sph>* deriv2) {
	const double r = sqrt(pow_2(p.R) + pow_2(p.z));
//    if(r==0)
//        throw std::runtime_error("PosDeriv Cyl=>Sph: r=0, degenerate case!");
	const double rinv= 1./r;
	const double costheta=p.z*rinv, sintheta=p.R*rinv;
	if(deriv!=NULL) {
		deriv->drdR=sintheta;
		deriv->drdz=costheta;
		deriv->dthetadR=costheta*rinv;
		deriv->dthetadz=-sintheta*rinv;
	}
	if(deriv2!=NULL) {
		deriv2->d2rdR2=pow_2(costheta)*rinv;
		deriv2->d2rdz2=pow_2(sintheta)*rinv;
		deriv2->d2rdRdz=-costheta*sintheta*rinv;
		deriv2->d2thetadR2=-2*costheta*sintheta*pow_2(rinv);
		deriv2->d2thetadz2=-deriv2->d2thetadR2;
		deriv2->d2thetadRdz=(pow_2(sintheta)-pow_2(costheta))*pow_2(rinv);
	}
	return PosSph(r, atan2(p.R, p.z), p.phi);
}

template<>
EXP PosCyl toPosDeriv(const PosSph& p, PosDerivT<Sph, Cyl>* deriv, PosDeriv2T<Sph, Cyl>* deriv2) {
    double sintheta, costheta;
    math::sincos(p.theta, sintheta, costheta);
    const double R=mul(p.r, sintheta), z=mul(p.r, costheta);
    if(deriv!=NULL) {
        deriv->dRdr=sintheta;
        deriv->dRdtheta=z;
        deriv->dzdr=costheta;
        deriv->dzdtheta=-R;
    }
    if(deriv2!=NULL) {
        deriv2->d2Rdrdtheta=costheta;
        deriv2->d2Rdtheta2=-p.r*sintheta;
        deriv2->d2zdrdtheta=-sintheta;
        deriv2->d2zdtheta2=-p.r*costheta;
    }
    return PosCyl(R, z, p.phi);
}

template<>
EXP PosCyl toPosDeriv(const PosProlSph& p, PosDerivT<ProlSph, Cyl>* deriv, PosDeriv2T<ProlSph, Cyl>* deriv2)
{
	const double absnu = fabs(p.nu);
	const double sign = p.nu>=0 ? 1 : -1;
	const double lminusd = p.lambda-p.coordsys.Delta2;
	const double nminusd = absnu-p.coordsys.Delta2;  // note: |nu|<=Delta^2
	if(nminusd>0 || lminusd<0)
		throw std::invalid_argument("Incorrect ProlSph coordinates");
	const double R = sqrt( lminusd * (1 - absnu / p.coordsys.Delta2) );
	const double z = sqrt( p.lambda * absnu /  p.coordsys.Delta2 ) * (p.nu>=0 ? 1 : -1);
	if(deriv!=NULL) {
		deriv->dRdlambda = 0.5*R/lminusd;
		deriv->dRdnu     = 0.5*R/nminusd * sign;
		deriv->dzdlambda = 0.5*z/p.lambda;
		deriv->dzdnu     = 0.5*z/p.nu;
	}
	if(deriv2!=NULL) {
		deriv2->d2Rdlambda2   = -0.25*R / pow_2(lminusd);
		deriv2->d2Rdnu2       = -0.25*R / pow_2(nminusd);
		deriv2->d2Rdlambdadnu = -0.25*R / (lminusd * nminusd * sign);
		deriv2->d2zdlambda2   = -0.25*z / pow_2(p.lambda);
		deriv2->d2zdnu2       = -0.25*z / pow_2(p.nu);
		deriv2->d2zdlambdadnu = -0.25*z / (p.lambda * p.nu);
	}
	return PosCyl(R, z, p.phi);
}

template<>
EXP PosProlSph toPosDeriv(const PosCyl& from, const ProlSph& cs,
			  PosDerivT<Cyl, ProlSph>* deriv, PosDeriv2T<Cyl, ProlSph>* deriv2)
{
    // lambda and nu are roots "t" of equation  R^2/(t-Delta^2) + z^2/t = 1
	double R2     = pow_2(from.R), z2 = pow_2(from.z);
	double signz  = from.z>=0 ? 1 : -1;   // nu will have the same sign as z
	double sum    = R2+z2+cs.Delta2;
	double dif    = R2+z2-cs.Delta2;
	double sqD    = sqrt(pow_2(dif) + 4*R2*cs.Delta2);   // determinant is always non-negative
	if(z2==0) sqD = sum;
	if(R2==0) sqD = fabs(dif);
	double lmd, dmn;  // lambda-Delta^2, Delta^2-|nu| - separately from lambda and nu, to avoid roundoffs
	if(dif >= 0) {
		lmd       = 0.5 * (sqD + dif);
		dmn       = R2>0 ? cs.Delta2 * R2 / lmd : 0;
	} else {
		dmn       = 0.5 * (sqD - dif);
		lmd       = cs.Delta2 * R2 / dmn;
	}
	double lambda = cs.Delta2 + lmd;
	double absnu  = 2 * cs.Delta2 / (sum + sqD) * z2;
	if(absnu*2 > cs.Delta2)             // compare |nu| and Delta^2-|nu|
		absnu     = cs.Delta2 - dmn;    // avoid roundoff errors when Delta^2-|nu| is small
	else
		dmn       = cs.Delta2 - absnu;  // same in the opposite case, when |nu| is small
	if(deriv!=NULL || deriv2!=NULL) {
		if(sqD==0)
			throw std::runtime_error("Error in coordinate conversion Cyl=>ProlSph: "
				"the special case lambda = nu = Delta^2 is not implemented");
		if(deriv!=NULL) {  // accurate expressions valid for arbitrary large/small values (no cancellations)
			deriv->dlambdadR = from.R * 2*lambda / sqD;
			deriv->dlambdadz = from.z * 2*lmd    / sqD;
			deriv->dnudR     = from.R * 2*-absnu / sqD * signz;
			deriv->dnudz     = from.z * 2*dmn    / sqD * signz;
		}
		if(deriv2!=NULL) {  // here no attempts were made to avoid cancellation errors
			double common = 8 * cs.Delta2 * R2 * z2 / pow_3(sqD);
			deriv2->d2lambdadR2 = 1 + sum/sqD - common;
			deriv2->d2lambdadz2 = 1 + dif/sqD + common;
			deriv2->d2nudR2     =(1 - sum/sqD + common) * signz;
			deriv2->d2nudz2     =(1 - dif/sqD - common) * signz;
			deriv2->d2lambdadRdz= 2 * from.R * from.z * (1 - sum * dif / pow_2(sqD)) / sqD;
			deriv2->d2nudRdz    = -deriv2->d2lambdadRdz * signz;
		}
	}
	return PosProlSph(lambda, absnu*signz, from.phi, cs);
}

template<>
EXP PosCyl toPosDeriv(const PosProlMod& p, PosDerivT<ProlMod, Cyl>* deriv, PosDeriv2T<ProlMod, Cyl>* deriv2)
{
	double sinv = (1 - pow_2(p.tau)) / (1 + pow_2(p.tau));
	double cosv = 2 * p.tau / (1 + pow_2(p.tau));
	if(deriv) {
		deriv->dRdrho = sinv;
		deriv->dRdtau =-p.rho * cosv * (1+sinv);
		deriv->dzdrho = p.rho * cosv / p.chi;
		deriv->dzdtau = p.chi * sinv * (1+sinv);
	}
	if(deriv2) {
		deriv2->d2Rdrho2    = 0;
		deriv2->d2Rdtau2    = p.rho * pow_2(1+sinv) * (1-2*sinv);
		deriv2->d2Rdrhodtau = -cosv * (1+sinv);
		deriv2->d2zdrho2    = (pow_2(p.chi) - pow_2(p.rho)) / pow_3(p.chi) * cosv;
		deriv2->d2zdtau2    = -p.chi * cosv * (1+sinv) * (1+2*sinv);
		deriv2->d2zdrhodtau = p.rho / p.chi * sinv * (1+sinv);
	}
	return PosCyl( p.rho * sinv, p.chi * cosv, p.phi);
}

template<>
EXP PosCyl toPosDeriv(const PosUVSph& p, PosDerivT<UVSph, Cyl>* deriv, PosDeriv2T<UVSph, Cyl>* deriv2)
{
	const double sn = sin(p.v), cs = cos(p.v);
	const double sh = sinh(p.u), ch = cosh(p.u);
	const double R = p.coordsys.Delta * sh * sn;
	const double z = p.coordsys.Delta * ch * cs;
	if(deriv!=NULL) {
		deriv->dRdu = p.coordsys.Delta * ch * sn;
		deriv->dRdv = p.coordsys.Delta * sh * cs;
		deriv->dzdu = p.coordsys.Delta * sh * cs;
		deriv->dzdv =-p.coordsys.Delta * ch * sn;
	}
	if(deriv2!=NULL) {
		deriv2->d2Rdu2  = R;
		deriv2->d2Rdv2  =-R;
		deriv2->d2Rdudv = z;
		deriv2->d2zdu2  = z;
		deriv2->d2zdv2  =-z;
		deriv2->d2zdudv =-R;
	}
	return PosCyl(R, z, p.phi);
}
template<>
EXP PosUVSph toPosDeriv(const PosCyl& from, const UVSph& cs,
			PosDerivT<Cyl, UVSph>* deriv, PosDeriv2T<Cyl, UVSph>* deriv2)
{
	double R2 = pow_2(from.R), z2 = pow_2(from.z);
	double R2_z2 = (R2+z2)/cs.Delta2;
	double shu2 = .5*(R2_z2-1 + sqrt(pow_2(1-R2_z2)+4*R2/cs.Delta2));
	double shu = sqrt(shu2);
	double chu2 = 1+shu2, chu=sqrt(chu2), ch2u=chu2+shu2;
	double cosv = from.z/(cs.Delta*chu), cosv2=cosv*cosv;
	double sinv = sqrt(1-cosv2), cos2v = cosv2 - pow_2(sinv);;
	double u = asinh(shu);
	double v = acos(cosv);
	if(deriv!=NULL){
		deriv->dudR = cs.Delta * chu * sinv / (cs.Delta2 * ch2u - R2 - z2);
		deriv->dudz = cs.Delta * shu * cosv / (cs.Delta2 * ch2u - R2 - z2);
		deriv->dvdR = cs.Delta * shu * cosv / (R2+z2-cs.Delta2*cos2v);
		deriv->dvdz =-cs.Delta * chu * sinv / (R2+z2-cs.Delta2*cos2v);
	}
	if(deriv2!=NULL) {
		//Missing code
	}
	return PosUVSph(u, v, from.phi, cs);
}

// fragment of code shared between toPosDeriv(Cyl=>ProlMod) and toPosVel(ProlMod=>Cyl)
inline void derivCyl2ProlMod(double rho, double chi, double sinv, double cosv,
    PosDerivT<Cyl, ProlMod>& deriv)
{
    double invdet= 1 / (pow_2(rho*cosv) + pow_2(chi*sinv)); // 1 / (rho^2 + (D*sinv)^2)
    deriv.drhodR = invdet * chi * sinv * chi;
    deriv.drhodz = invdet * rho * cosv * chi;
    deriv.dtaudR =-invdet * rho * cosv / (1 + sinv);
    deriv.dtaudz = invdet * chi * sinv / (1 + sinv);
}

template<>
EXP PosProlMod toPosDeriv(const PosCyl& p, const ProlMod& cs,
    PosDerivT<Cyl, ProlMod>* deriv, PosDeriv2T<Cyl, ProlMod>* deriv2)
{
    double r2  = pow_2(p.R) + pow_2(p.z);
    double sum = 0.5 * (r2 + pow_2(cs.D));
    double dif = 0.5 * (r2 - pow_2(cs.D));
    double det = p.z==0 || cs.D==0 ? sum : sqrt(pow_2(dif) + pow_2(p.R*cs.D));
    double chi = sqrt(det + sum);
    double cosv= p.z / chi;
    double rho, sinv;  // accurate treatment to avoid cancellation
    if(dif >= 0) {
        rho  = sqrt(det + dif);
        sinv = p.R!=0 ? p.R / rho : 0;
    } else {
        sinv = sqrt(det - dif) / cs.D;
        rho  = p.R / sinv;
    }
    if(deriv)
        derivCyl2ProlMod(rho, chi, sinv, cosv, *deriv);
    if(deriv2) {
        assert(!"Deriv2 Cyl=>ProlMod not implemented");
    }
    return PosProlMod(rho, cosv / (1 + sinv), p.phi, chi);
}
template<>
EXP PosCyl toPosDerivR(const PosSphMod& p, PosDerivT<Cyl, SphMod>* derivs,
		       PosDeriv2T<Cyl, SphMod>* derivs2){
	double tau2=pow_2(p.tau);
	double sintheta = (1-tau2)/(1+tau2);
	double costheta = p.tau*(1+sintheta);
	PosCyl Rz(p.r*sintheta, p.r*costheta, p.phi);
	if(derivs){ 
		derivs->drdR = sintheta;
		derivs->dtaudR = -p.tau/p.r;
		derivs->drdz = costheta;
		derivs->dtaudz = sintheta/(1+sintheta)/p.r;
	}
	return Rz;
}

template<>
EXP PosCyl toPosDerivR(const PosSph& p, PosDerivT<Cyl, Sph>* derivs,
		       PosDeriv2T<Cyl, Sph>* derivs2){
	double costheta, sintheta;
	math::sincos(p.theta, sintheta, costheta);
	PosCyl Rz(p.r*sintheta, p.r*costheta, p.phi);
	if(derivs){ 
		derivs->drdR = sintheta;
		derivs->dthetadR = costheta/p.r;
		derivs->drdz = costheta;
		derivs->dthetadz = -sintheta/p.r;
	}
	return Rz;
}

template<>
EXP PosSph toPosDerivR(const PosCyl& p, PosDerivT<Sph, Cyl>* derivs,
		       PosDeriv2T<Sph, Cyl>* derivs2){
	double r=sqrt(p.R*p.R+p.z*p.z), costheta=p.z/r;
	double theta=acos(costheta), sintheta=sin(theta);
	PosSph rtheta(r, theta, p.phi);
	if(derivs){ 
		derivs->dRdr = sintheta;
		derivs->dRdtheta = r*costheta;
		derivs->dzdr = costheta;
		derivs->dzdtheta = -r*sintheta;
	}
	return rtheta;
}

template<>
EXP PosSphMod toPosDerivR(const PosCyl& p, PosDerivT<SphMod, Cyl>* derivs,
			 PosDeriv2T<SphMod, Cyl>* derivs2){
	double theta=atan2(p.R,p.z);
	double cs = cos(theta), sn = sin(theta);
	double tau = cs/(1+sn), r = sqrt(pow_2(p.R)+pow_2(p.z)), dthdtau = -(1+sn);
	PosSphMod rtau(r,tau,p.phi);
	if(derivs){ 
		derivs->dRdr = sn;
		derivs->dzdr = cs;
		derivs->dRdtau = p.z*dthdtau;
		derivs->dzdtau = -p.R*dthdtau;
	}
	return rtau;
}


//--------  position+velocity conversion functions  ---------//
EXP void UVmomenta(const PosUVSph& uv,const VelCyl& Rz,double& pu,double& pv){
	double chu = cosh(uv.u), shu = sinh(uv.u), csv = cos(uv.v), snv = sin(uv.v);
	pu = uv.coordsys.Delta * (chu*snv*Rz.vR + shu*csv*Rz.vz);
	pv = uv.coordsys.Delta * (shu*csv*Rz.vR - chu*snv*Rz.vz);
}
EXP void UVmomenta(const PosVelUVSph& uvt,double& pu,double& pv){
	double fac = .5*uvt.coordsys.Delta2*(pow_2(cosh(uvt.u)) - pow_2(sin(uvt.v)));
	pu = fac*uvt.udot; pv = fac*uvt.vdot;
}
EXP void Rzmomenta(const PosUVSph& uv,const double pu,const double pv,VelCyl& Rz){
	double chu = cosh(uv.u), shu = sinh(uv.u), csv = cos(uv.v), snv = sin(uv.v);
	Rz.vR = (pu*chu*snv + pv*shu*csv)/(uv.coordsys.Delta*(snv*snv + shu*shu));
	Rz.vz = (pu*shu*csv - pv*chu*snv)/(uv.coordsys.Delta*(chu*chu - csv*csv));
}

template<>
EXP PosVelCar toPosVel(const PosVelCyl& p) {
    double sinphi, cosphi;
    math::sincos(p.phi, sinphi, cosphi);
    const double vx=p.vR*cosphi-p.vphi*sinphi;
    const double vy=p.vR*sinphi+p.vphi*cosphi;
    return PosVelCar(p.R*cosphi, p.R*sinphi, p.z, vx, vy, p.vz);
}

template<>
EXP PosVelCar toPosVel(const PosVelSph& p) {
    double sintheta, costheta, sinphi, cosphi;
    math::sincos(p.theta, sintheta, costheta);
    math::sincos(p.phi, sinphi, cosphi);
    const double R=p.r*sintheta, vmer=p.vr*sintheta + p.vtheta*costheta;
    const double vx=vmer*cosphi - p.vphi*sinphi;
    const double vy=vmer*sinphi + p.vphi*cosphi;
    const double vz=p.vr*costheta - p.vtheta*sintheta;
    return PosVelCar(R*cosphi, R*sinphi, p.r*costheta, vx, vy, vz); 
}

template<>
EXP PosVelCyl toPosVel(const PosVelCar& p) {
	const double R=sqrt(pow_2(p.x)+pow_2(p.y));
	if(R==0)  // determine phi from vy/vx rather than y/x
		return PosVelCyl(R, p.z, atan2(p.vy, p.vx), sqrt(pow_2(p.vx)+pow_2(p.vy)), p.vz, 0);
	const double cosphi=p.x/R, sinphi=p.y/R;
	const double vR  = p.vx*cosphi+p.vy*sinphi;
	const double vphi=-p.vx*sinphi+p.vy*cosphi;
	return PosVelCyl(R, p.z, atan2(p.y, p.x), vR, p.vz, vphi);
}

template<>
EXP PosVelCyl toPosVel(const PosVelSph& p) {
	double sintheta, costheta;
	math::sincos(p.theta, sintheta, costheta);
	const double R=p.r*sintheta, z=p.r*costheta;
	const double vR=p.vr*sintheta+p.vtheta*costheta;
	const double vz=p.vr*costheta-p.vtheta*sintheta;
	return PosVelCyl(R, z, p.phi, vR, vz, p.vphi);
}

template<>
EXP PosVelCyl toPosVel(const PosMomSph& p) {
	double sintheta, costheta;
	math::sincos(p.theta, sintheta, costheta);
	const double R=p.r*sintheta, z=p.r*costheta;
	const double vR=p.pr*sintheta+p.ptheta/p.r*costheta;
	const double vz=p.pr*costheta-p.ptheta/p.r*sintheta;
	return PosVelCyl(R, z, p.phi, vR, vz, p.pphi/(p.r*sintheta));
}

template<>
EXP PosMomSph toPosMom(const PosMomCyl& p) {
	PosDerivT<Sph, Cyl> derivs;	
	const PosSph rtheta(toPosDerivR<Cyl, Sph> (p, &derivs));
	return PosMomSph(rtheta, MomSph(
					derivs.dRdr * p.pR + derivs.dzdr * p.pz,
					derivs.dRdtheta * p.pR + derivs.dzdtheta * p.pz,
					p.pphi));
}

template<>
EXP PosMomSphMod toPosMom(const PosMomCyl& p) {
	const double r   = sqrt(pow_2(p.R) + pow_2(p.z));
	const double tau = p.z / (p.R + r);
	const double pr  = (p.R * p.pR + p.z * p.pz) / r;
	const double ptau= (p.R * p.pz - p.z * p.pR) * (1 + p.R/r);
	return PosMomSphMod(r, tau, p.phi, pr, ptau, p.pphi);
}

template<>
EXP PosVelSph toPosVel(const PosVelCar& p) {
    const double R2=pow_2(p.x)+pow_2(p.y), R=sqrt(R2);
    const double r2=R2+pow_2(p.z), r=sqrt(r2), invr=1/r;
    if(R==0) {
        const double vR=sqrt(pow_2(p.vx)+pow_2(p.vy));
        const double phi=atan2(p.vy, p.vx);
        if(p.z==0) 
            return PosVelSph(0, atan2(vR, p.vz), phi, sqrt(vR*vR+p.vz*p.vz), 0, 0);
        return PosVelSph(r, p.z>0?0:M_PI, phi, p.vz*(p.z>0?1:-1), vR*(p.z>0?1:-1), 0);
    }
    const double temp   = p.x*p.vx+p.y*p.vy;
    const double vr     = (temp+p.z*p.vz)*invr;
    const double vtheta = (temp*p.z/R-p.vz*R)*invr;
    const double vphi   = (p.x*p.vy-p.y*p.vx)/R;
    return PosVelSph(r, atan2(R, p.z), atan2(p.y, p.x), vr, vtheta, vphi);
}

template<>
EXP PosVelSph toPosVel(const PosVelCyl& p) {
    const double r=sqrt(pow_2(p.R)+pow_2(p.z));
    if(r==0) {
        return PosVelSph(0, atan2(p.vR, p.vz), p.phi, sqrt(p.vR*p.vR+p.vz*p.vz), 0, 0);
    }
    const double rinv= 1./r;
    const double costheta=p.z*rinv, sintheta=p.R*rinv;
    const double vr=p.vR*sintheta+p.vz*costheta;
    const double vtheta=p.vR*costheta-p.vz*sintheta;
    return PosVelSph(r, atan2(p.R, p.z), p.phi, vr, vtheta, p.vphi);
}

template<>
EXP PosVelProlSph toPosVel(const PosVelCyl& from, const ProlSph& cs) {
	PosDerivT<Cyl, ProlSph> derivs;
	const PosProlSph pprol = toPosDeriv<Cyl, ProlSph> (from, cs, &derivs);
	double lambdadot = derivs.dlambdadR*from.vR + derivs.dlambdadz*from.vz;
	double nudot     = derivs.dnudR    *from.vR + derivs.dnudz    *from.vz;
	double phidot    = from.vphi!=0 ? from.vphi/from.R : 0;
	return PosVelProlSph(pprol, lambdadot, nudot, phidot);
}

template<>
EXP PosMomProlMod toPosMom(const PosMomCyl& p, const ProlMod& cs) {
	const PosProlMod pprol = toPosDeriv<Cyl, ProlMod>(p, cs, NULL);
	PosDerivT<ProlMod, Cyl> derivs;  // need derivs of _inverse_ transformation
	toPosDeriv<ProlMod, Cyl>(pprol, &derivs);
	double prho = derivs.dRdrho * p.pR + derivs.dzdrho * p.pz;
	double ptau = derivs.dRdtau * p.pR + derivs.dzdtau * p.pz;
	double pphi = p.pphi;
	return PosMomProlMod(pprol, MomProlMod(prho, ptau, pphi));
}

template<>
EXP PosMomUVSph toPosMom(const PosMomCyl& p, const UVSph& cs) {
	const PosUVSph pUVSph = toPosDeriv<Cyl, UVSph>(p, cs, NULL);
	PosDerivT<UVSph, Cyl> derivs;  // need derivs of _inverse_ transformation
	toPosDeriv<UVSph, Cyl>(pUVSph, &derivs);
	double pu = derivs.dRdu * p.pR + derivs.dzdu * p.pz;
	double pv = derivs.dRdv * p.pR + derivs.dzdv * p.pz;
	double pphi = p.pphi;
	return PosMomUVSph(pUVSph, MomUVSph(pu, pv, pphi));
}

template<>
EXP PosMomCyl toPosMom(const PosMomUVSph& pUVSph) {
	PosDerivT<UVSph, Cyl> derivs;  // need derivs of transformation
	const PosCyl p = toPosDeriv<UVSph, Cyl>(pUVSph, &derivs);
	double det = derivs.dRdu * derivs.dzdv - derivs.dRdv * derivs.dzdu;
	double pR = (derivs.dzdv * pUVSph.pu - derivs.dzdu * pUVSph.pv)/det;
	double pz = (derivs.dRdu * pUVSph.pv - derivs.dRdv * pUVSph.pu)/det;
	return PosMomCyl(p, MomCyl(pR, pz, pUVSph.pphi));
}

template<>
EXP PosMomCyl toPosMom(const PosMomSphMod& p) {
	PosDerivT<Cyl, SphMod> derivs;  // need derivs of _inverse_ transformation
	const PosCyl Rz(toPosDerivR<SphMod, Cyl> (p, &derivs));
	return PosMomCyl(Rz, MomCyl(
				    derivs.drdR * p.pr + derivs.dtaudR * p.ptau,
				    derivs.drdz * p.pr + derivs.dtaudz * p.ptau,
				    p.pphi));
}

template<>
EXP PosMomCyl toPosMom(const PosMomSph& p) {
	PosDerivT<Cyl, Sph> derivs;  // need derivs of _inverse_ transformation
	const PosCyl Rz(toPosDerivR<Sph, Cyl> (p, &derivs));
	return PosMomCyl(Rz, MomCyl(
				    derivs.drdR * p.pr + derivs.dthetadR * p.ptheta,
				    derivs.drdz * p.pr + derivs.dthetadz * p.ptheta,
				    p.pphi));
}

template<>
EXP PosMomCyl toPosMom(const PosMomProlMod& p) {
	double sinv = (1 - pow_2(p.tau)) / (1 + pow_2(p.tau));
	double cosv = 2 * p.tau / (1 + pow_2(p.tau));
	PosDerivT<Cyl, ProlMod> derivs;  // need derivs of _inverse_ transformation
	derivCyl2ProlMod(p.rho, p.chi, sinv, cosv, derivs);
	return PosMomCyl(p.rho * sinv, p.chi * cosv, p.phi,
			 derivs.drhodR * p.prho + derivs.dtaudR * p.ptau,
			 derivs.drhodz * p.prho + derivs.dtaudz * p.ptau,
			 p.pphi);
}

template<>
EXP PosMomCar toPosMom(const PosMomCyl& p) {
	double csp=cos(p.phi), snp=sin(p.phi), vphi=p.pphi/p.R;
	return PosMomCar(p.R*csp, p.R*snp, p.z,
			 p.pR*csp-vphi*snp,p.pR*snp+vphi*csp,p.pz);
}

template<>
EXP PosMomCar toPosMom(const PosMomSph& p) {
	PosMomCyl pC = toPosMom<Sph,Cyl>(p);
	return toPosMom<Cyl, Car>(pC);
}

template<>
EXP PosVelUVSph toPosVel(const PosVelCyl& from, const UVSph& cs) {
	PosDerivT<Cyl, UVSph> derivs;
	const PosUVSph puv = toPosDeriv<Cyl, UVSph> (from, cs, &derivs);
	double udot = derivs.dudR * from.vR + derivs.dudz * from.vz;
	double vdot     = derivs.dvdR * from.vR + derivs.dvdz * from.vz;
	double phidot    = from.vphi!=0 ? from.vphi/from.R : 0;
	return PosVelUVSph(puv, udot, vdot, phidot);
}

template<>
EXP PosVelCyl toPosVel(const PosVelUVSph& from) {
	PosDerivT<UVSph, Cyl> derivs;
	const PosCyl Rz = toPosDeriv<UVSph, Cyl> (from, &derivs);
	double vR = derivs.dRdu * from.udot + derivs.dRdv * from.vdot;
	double vz = derivs.dzdu * from.udot + derivs.dzdv * from.vdot;
	double vphi = from.phidot!=0 ? from.phidot*Rz.R : 0;
	const VelCyl vRz(vR,vz,vphi);
	return PosVelCyl(Rz, vRz);
}

template<>
EXP PosVelCyl toPosVel(const PosMomCyl& p) {
	return PosVelCyl(p.R,p.z,p.phi,p.pR,p.pz,p.pphi/p.R);
}
/** trivial conversions */
template<> inline PosMomCar toPosMom<Car,Car>(const PosMomCar& p) { return p;}
template<> inline PosMomCyl toPosMom<Cyl,Cyl>(const PosMomCyl& p) { return p;}
template<> inline PosMomSph toPosMom<Sph,Sph>(const PosMomSph& p) { return p;}
/* Vel to Mom conversions */
template<>
EXP PosMomCar toPosMom(const PosVelCar& p) {
	return PosMomCar(p.x, p.y, p.z, p.vx, p.vy, p.vz);
}

template<>
EXP PosMomCyl toPosMom(const PosVelCyl& p) {
	return PosMomCyl(p.R, p.z, p.phi, p.vR, p.vz, p.vphi*p.R);
}

template<>
EXP PosMomSph toPosMom(const PosVelSph& p) {
	return PosMomSph(p.r, p.theta, p.phi, p.vr, p.vtheta*p.r, p.vphi*p.r*sin(p.theta));
}


//-------- implementations of functions that convert gradients --------//
// note: the code below is machine-generated

template<>
EXP GradCar toGrad(const GradCyl& src, const PosDerivT<Car, Cyl>& deriv) {
    GradCar dest;
    dest.dx = src.dR*deriv.dRdx + src.dphi*deriv.dphidx;
    dest.dy = src.dR*deriv.dRdy + src.dphi*deriv.dphidy;
    dest.dz = src.dz;
    return dest;
}

template<>
EXP GradCar toGrad(const GradSph& src, const PosDerivT<Car, Sph>& deriv) {
    GradCar dest;
    dest.dx = src.dr*deriv.drdx + src.dtheta*deriv.dthetadx + src.dphi*deriv.dphidx;
    dest.dy = src.dr*deriv.drdy + src.dtheta*deriv.dthetady + src.dphi*deriv.dphidy;
    dest.dz = src.dr*deriv.drdz + src.dtheta*deriv.dthetadz;
    return dest;
}

template<>
EXP GradCyl toGrad(const GradCar& src, const PosDerivT<Cyl, Car>& deriv) {
    GradCyl dest;
    dest.dR = src.dx*deriv.dxdR + src.dy*deriv.dydR;
    dest.dz = src.dz;
    dest.dphi = src.dx*deriv.dxdphi + src.dy*deriv.dydphi;
    return dest;
}

template<>
EXP GradCyl toGrad(const GradSph& src, const PosDerivT<Cyl, Sph>& deriv) {
    GradCyl dest;
    dest.dR = src.dr*deriv.drdR + src.dtheta*deriv.dthetadR;
    dest.dz = src.dr*deriv.drdz + src.dtheta*deriv.dthetadz;
    dest.dphi = src.dphi;
    return dest;
}

template<>
EXP GradSph toGrad(const GradCar& src, const PosDerivT<Sph, Car>& deriv) {
    GradSph dest;
    dest.dr = src.dx*deriv.dxdr + src.dy*deriv.dydr + src.dz*deriv.dzdr;
    dest.dtheta = src.dx*deriv.dxdtheta + src.dy*deriv.dydtheta + src.dz*deriv.dzdtheta;
    dest.dphi = src.dx*deriv.dxdphi + src.dy*deriv.dydphi;
    return dest;
}

template<>
EXP GradSph toGrad(const GradCyl& src, const PosDerivT<Sph, Cyl>& deriv) {
    GradSph dest;
    dest.dr = src.dR*deriv.dRdr + src.dz*deriv.dzdr;
    dest.dtheta = src.dR*deriv.dRdtheta + src.dz*deriv.dzdtheta;
    dest.dphi = src.dphi;
    return dest;
}

template<>
EXP GradCyl toGrad(const GradProlSph& src, const PosDerivT<Cyl, ProlSph>& deriv) {
    GradCyl dest;
    dest.dR   = src.dlambda*deriv.dlambdadR + src.dnu*deriv.dnudR;
    dest.dz   = src.dlambda*deriv.dlambdadz + src.dnu*deriv.dnudz;
    dest.dphi = src.dphi;
    return dest;
}

template<>
EXP GradProlSph toGrad(const GradCyl& src, const PosDerivT<ProlSph, Cyl>& deriv) {
    GradProlSph dest;
    dest.dlambda = src.dR*deriv.dRdlambda + src.dz*deriv.dzdlambda;
    dest.dnu     = src.dR*deriv.dRdnu     + src.dz*deriv.dzdnu;
    dest.dphi    = src.dphi;
    return dest;
}

template<>
EXP GradCyl toGrad(const GradProlMod& src, const PosDerivT<Cyl, ProlMod>& deriv) {
    GradCyl dest;
    dest.dR   = src.drho*deriv.drhodR + src.dtau*deriv.dtaudR;
    dest.dz   = src.drho*deriv.drhodz + src.dtau*deriv.dtaudz;
    dest.dphi = src.dphi;
    return dest;
}

template<>
EXP GradProlMod toGrad(const GradCyl& src, const PosDerivT<ProlMod, Cyl>& deriv) {
    GradProlMod dest;
    dest.drho = src.dR*deriv.dRdrho + src.dz*deriv.dzdrho;
    dest.dtau = src.dR*deriv.dRdtau + src.dz*deriv.dzdtau;
    dest.dphi = src.dphi;
    return dest;
}

//-------- implementations of functions that convert hessians --------//
// note: the code below is machine-generated and is not intended to be human-readable

template<>
EXP HessCar toHess(const GradCyl& srcGrad, const HessCyl& srcHess,
    const PosDerivT<Car, Cyl>& deriv, const PosDeriv2T<Car, Cyl>& deriv2) {
    HessCar dest;
    dest.dx2 = 
        (srcHess.dR2*deriv.dRdx + srcHess.dRdphi*deriv.dphidx)*deriv.dRdx + 
        (srcHess.dRdphi*deriv.dRdx + srcHess.dphi2*deriv.dphidx)*deriv.dphidx + 
        srcGrad.dR*deriv2.d2Rdx2 + srcGrad.dphi*deriv2.d2phidx2;
    dest.dxdy = 
        (srcHess.dR2*deriv.dRdy + srcHess.dRdphi*deriv.dphidy)*deriv.dRdx + 
        (srcHess.dRdphi*deriv.dRdy + srcHess.dphi2*deriv.dphidy)*deriv.dphidx + 
        srcGrad.dR*deriv2.d2Rdxdy + srcGrad.dphi*deriv2.d2phidxdy;
    dest.dxdz = 
        srcHess.dRdz*deriv.dRdx + 
        srcHess.dzdphi*deriv.dphidx;
    dest.dy2 = 
        (srcHess.dR2*deriv.dRdy + srcHess.dRdphi*deriv.dphidy)*deriv.dRdy + 
        (srcHess.dRdphi*deriv.dRdy + srcHess.dphi2*deriv.dphidy)*deriv.dphidy + 
        srcGrad.dR*deriv2.d2Rdy2 + srcGrad.dphi*deriv2.d2phidy2;
    dest.dydz = 
        srcHess.dRdz*deriv.dRdy + 
        srcHess.dzdphi*deriv.dphidy;
    dest.dz2 = srcHess.dz2;
    return dest;
}

template<>
EXP HessCar toHess(const GradSph& srcGrad, const HessSph& srcHess,
    const PosDerivT<Car, Sph>& deriv, const PosDeriv2T<Car, Sph>& deriv2) {
    HessCar dest;
    dest.dx2 = 
        (srcHess.dr2*deriv.drdx + srcHess.drdtheta*deriv.dthetadx + srcHess.drdphi*deriv.dphidx)*deriv.drdx + 
        (srcHess.drdtheta*deriv.drdx + srcHess.dtheta2*deriv.dthetadx + srcHess.dthetadphi*deriv.dphidx)*deriv.dthetadx + 
        (srcHess.drdphi*deriv.drdx + srcHess.dthetadphi*deriv.dthetadx + srcHess.dphi2*deriv.dphidx)*deriv.dphidx + 
        srcGrad.dr*deriv2.d2rdx2 + srcGrad.dtheta*deriv2.d2thetadx2 + srcGrad.dphi*deriv2.d2phidx2;
    dest.dxdy = 
        (srcHess.dr2*deriv.drdy + srcHess.drdtheta*deriv.dthetady + srcHess.drdphi*deriv.dphidy)*deriv.drdx + 
        (srcHess.drdtheta*deriv.drdy + srcHess.dtheta2*deriv.dthetady + srcHess.dthetadphi*deriv.dphidy)*deriv.dthetadx + 
        (srcHess.drdphi*deriv.drdy + srcHess.dthetadphi*deriv.dthetady + srcHess.dphi2*deriv.dphidy)*deriv.dphidx + 
        srcGrad.dr*deriv2.d2rdxdy + srcGrad.dtheta*deriv2.d2thetadxdy + srcGrad.dphi*deriv2.d2phidxdy;
    dest.dxdz = 
        (srcHess.dr2*deriv.drdz + srcHess.drdtheta*deriv.dthetadz)*deriv.drdx + 
        (srcHess.drdtheta*deriv.drdz + srcHess.dtheta2*deriv.dthetadz)*deriv.dthetadx + 
        (srcHess.drdphi*deriv.drdz + srcHess.dthetadphi*deriv.dthetadz)*deriv.dphidx + 
        srcGrad.dr*deriv2.d2rdxdz + srcGrad.dtheta*deriv2.d2thetadxdz;
    dest.dy2 = 
        (srcHess.dr2*deriv.drdy + srcHess.drdtheta*deriv.dthetady + srcHess.drdphi*deriv.dphidy)*deriv.drdy + 
        (srcHess.drdtheta*deriv.drdy + srcHess.dtheta2*deriv.dthetady + srcHess.dthetadphi*deriv.dphidy)*deriv.dthetady + 
        (srcHess.drdphi*deriv.drdy + srcHess.dthetadphi*deriv.dthetady + srcHess.dphi2*deriv.dphidy)*deriv.dphidy + 
        srcGrad.dr*deriv2.d2rdy2 + srcGrad.dtheta*deriv2.d2thetady2 + srcGrad.dphi*deriv2.d2phidy2;
    dest.dydz = 
        (srcHess.dr2*deriv.drdz + srcHess.drdtheta*deriv.dthetadz)*deriv.drdy + 
        (srcHess.drdtheta*deriv.drdz + srcHess.dtheta2*deriv.dthetadz)*deriv.dthetady + 
        (srcHess.drdphi*deriv.drdz + srcHess.dthetadphi*deriv.dthetadz)*deriv.dphidy + 
        srcGrad.dr*deriv2.d2rdydz + srcGrad.dtheta*deriv2.d2thetadydz;
    dest.dz2 = 
        (srcHess.dr2*deriv.drdz + srcHess.drdtheta*deriv.dthetadz)*deriv.drdz + 
        (srcHess.drdtheta*deriv.drdz + srcHess.dtheta2*deriv.dthetadz)*deriv.dthetadz + 
        srcGrad.dr*deriv2.d2rdz2 + srcGrad.dtheta*deriv2.d2thetadz2;
    return dest;
}

template<>
EXP HessCyl toHess(const GradCar& srcGrad, const HessCar& srcHess,
    const PosDerivT<Cyl, Car>& deriv, const PosDeriv2T<Cyl, Car>& deriv2) {
    HessCyl dest;
    dest.dR2 = 
        (srcHess.dx2*deriv.dxdR + srcHess.dxdy*deriv.dydR)*deriv.dxdR + 
        (srcHess.dxdy*deriv.dxdR + srcHess.dy2*deriv.dydR)*deriv.dydR;
    dest.dRdz = 
        srcHess.dxdz*deriv.dxdR + 
        srcHess.dydz*deriv.dydR;
    dest.dRdphi = 
        (srcHess.dx2*deriv.dxdphi + srcHess.dxdy*deriv.dydphi)*deriv.dxdR + 
        (srcHess.dxdy*deriv.dxdphi + srcHess.dy2*deriv.dydphi)*deriv.dydR + 
        srcGrad.dx*deriv2.d2xdRdphi + srcGrad.dy*deriv2.d2ydRdphi;
    dest.dz2 = srcHess.dz2;
    dest.dzdphi = 
        (srcHess.dxdz*deriv.dxdphi + srcHess.dydz*deriv.dydphi);
    dest.dphi2 = 
        (srcHess.dx2*deriv.dxdphi + srcHess.dxdy*deriv.dydphi)*deriv.dxdphi + 
        (srcHess.dxdy*deriv.dxdphi + srcHess.dy2*deriv.dydphi)*deriv.dydphi + 
        srcGrad.dx*deriv2.d2xdphi2 + srcGrad.dy*deriv2.d2ydphi2;
    return dest;
}

template<>
EXP HessCyl toHess(const GradSph& srcGrad, const HessSph& srcHess,
    const PosDerivT<Cyl, Sph>& deriv, const PosDeriv2T<Cyl, Sph>& deriv2) {
    HessCyl dest;
    dest.dR2 = 
        (srcHess.dr2*deriv.drdR + srcHess.drdtheta*deriv.dthetadR)*deriv.drdR + 
        (srcHess.drdtheta*deriv.drdR + srcHess.dtheta2*deriv.dthetadR)*deriv.dthetadR + 
        srcGrad.dr*deriv2.d2rdR2 + srcGrad.dtheta*deriv2.d2thetadR2;
    dest.dRdz = 
        (srcHess.dr2*deriv.drdz + srcHess.drdtheta*deriv.dthetadz)*deriv.drdR + 
        (srcHess.drdtheta*deriv.drdz + srcHess.dtheta2*deriv.dthetadz)*deriv.dthetadR + 
        srcGrad.dr*deriv2.d2rdRdz + srcGrad.dtheta*deriv2.d2thetadRdz;
    dest.dRdphi = 
        srcHess.drdphi*deriv.drdR + 
        srcHess.dthetadphi*deriv.dthetadR;
    dest.dz2 = 
        (srcHess.dr2*deriv.drdz + srcHess.drdtheta*deriv.dthetadz)*deriv.drdz + 
        (srcHess.drdtheta*deriv.drdz + srcHess.dtheta2*deriv.dthetadz)*deriv.dthetadz + 
        srcGrad.dr*deriv2.d2rdz2 + srcGrad.dtheta*deriv2.d2thetadz2;
    dest.dzdphi = 
        srcHess.drdphi*deriv.drdz + 
        srcHess.dthetadphi*deriv.dthetadz;
    dest.dphi2 = srcHess.dphi2;
    return dest;
}

template<>
EXP HessSph toHess(const GradCar& srcGrad, const HessCar& srcHess,
    const PosDerivT<Sph, Car>& deriv, const PosDeriv2T<Sph, Car>& deriv2) {
    HessSph dest;
    dest.dr2 = 
        (srcHess.dx2*deriv.dxdr + srcHess.dxdy*deriv.dydr + srcHess.dxdz*deriv.dzdr)*deriv.dxdr + 
        (srcHess.dxdy*deriv.dxdr + srcHess.dy2*deriv.dydr + srcHess.dydz*deriv.dzdr)*deriv.dydr + 
        (srcHess.dxdz*deriv.dxdr + srcHess.dydz*deriv.dydr + srcHess.dz2*deriv.dzdr)*deriv.dzdr;
    dest.drdtheta = 
        (srcHess.dx2*deriv.dxdtheta + srcHess.dxdy*deriv.dydtheta + srcHess.dxdz*deriv.dzdtheta)*deriv.dxdr + 
        (srcHess.dxdy*deriv.dxdtheta + srcHess.dy2*deriv.dydtheta + srcHess.dydz*deriv.dzdtheta)*deriv.dydr + 
        (srcHess.dxdz*deriv.dxdtheta + srcHess.dydz*deriv.dydtheta + srcHess.dz2*deriv.dzdtheta)*deriv.dzdr + 
        srcGrad.dx*deriv2.d2xdrdtheta + srcGrad.dy*deriv2.d2ydrdtheta + srcGrad.dz*deriv2.d2zdrdtheta;
    dest.drdphi = 
        (srcHess.dx2*deriv.dxdphi + srcHess.dxdy*deriv.dydphi)*deriv.dxdr + 
        (srcHess.dxdy*deriv.dxdphi + srcHess.dy2*deriv.dydphi)*deriv.dydr + 
        (srcHess.dxdz*deriv.dxdphi + srcHess.dydz*deriv.dydphi)*deriv.dzdr + 
        srcGrad.dx*deriv2.d2xdrdphi + srcGrad.dy*deriv2.d2ydrdphi;
    dest.dtheta2 = 
        (srcHess.dx2*deriv.dxdtheta + srcHess.dxdy*deriv.dydtheta + srcHess.dxdz*deriv.dzdtheta)*deriv.dxdtheta + 
        (srcHess.dxdy*deriv.dxdtheta + srcHess.dy2*deriv.dydtheta + srcHess.dydz*deriv.dzdtheta)*deriv.dydtheta + 
        (srcHess.dxdz*deriv.dxdtheta + srcHess.dydz*deriv.dydtheta + srcHess.dz2*deriv.dzdtheta)*deriv.dzdtheta + 
        srcGrad.dx*deriv2.d2xdtheta2 + srcGrad.dy*deriv2.d2ydtheta2 + srcGrad.dz*deriv2.d2zdtheta2;
    dest.dthetadphi = 
        (srcHess.dx2*deriv.dxdphi + srcHess.dxdy*deriv.dydphi)*deriv.dxdtheta + 
        (srcHess.dxdy*deriv.dxdphi + srcHess.dy2*deriv.dydphi)*deriv.dydtheta + 
        (srcHess.dxdz*deriv.dxdphi + srcHess.dydz*deriv.dydphi)*deriv.dzdtheta + 
        srcGrad.dx*deriv2.d2xdthetadphi + srcGrad.dy*deriv2.d2ydthetadphi;
    dest.dphi2 = 
        (srcHess.dx2*deriv.dxdphi + srcHess.dxdy*deriv.dydphi)*deriv.dxdphi + 
        (srcHess.dxdy*deriv.dxdphi + srcHess.dy2*deriv.dydphi)*deriv.dydphi + 
        srcGrad.dx*deriv2.d2xdphi2 + srcGrad.dy*deriv2.d2ydphi2;
    return dest;
}

template<>
EXP HessSph toHess(const GradCyl& srcGrad, const HessCyl& srcHess,
    const PosDerivT<Sph, Cyl>& deriv, const PosDeriv2T<Sph, Cyl>& deriv2) {
    HessSph dest;
    dest.dr2 = 
        (srcHess.dR2*deriv.dRdr + srcHess.dRdz*deriv.dzdr)*deriv.dRdr + 
        (srcHess.dRdz*deriv.dRdr + srcHess.dz2*deriv.dzdr)*deriv.dzdr;
    dest.drdtheta = 
        (srcHess.dR2*deriv.dRdtheta + srcHess.dRdz*deriv.dzdtheta)*deriv.dRdr + 
        (srcHess.dRdz*deriv.dRdtheta + srcHess.dz2*deriv.dzdtheta)*deriv.dzdr + 
        srcGrad.dR*deriv2.d2Rdrdtheta + srcGrad.dz*deriv2.d2zdrdtheta;
    dest.drdphi = 
        srcHess.dRdphi*deriv.dRdr + 
        srcHess.dzdphi*deriv.dzdr;
    dest.dtheta2 = 
        (srcHess.dR2*deriv.dRdtheta + srcHess.dRdz*deriv.dzdtheta)*deriv.dRdtheta + 
        (srcHess.dRdz*deriv.dRdtheta + srcHess.dz2*deriv.dzdtheta)*deriv.dzdtheta + 
        srcGrad.dR*deriv2.d2Rdtheta2 + srcGrad.dz*deriv2.d2zdtheta2;
    dest.dthetadphi = 
        srcHess.dRdphi*deriv.dRdtheta + 
        srcHess.dzdphi*deriv.dzdtheta;
    dest.dphi2 = srcHess.dphi2;
    return dest;
}

template<>
EXP HessCyl toHess(const GradProlSph& srcGrad, const HessProlSph& srcHess,
    const PosDerivT<Cyl, ProlSph>& deriv, const PosDeriv2T<Cyl, ProlSph>& deriv2) {
    HessCyl dest;
    dest.dR2 = 
        (srcHess.dlambda2*deriv.dlambdadR + srcHess.dlambdadnu*deriv.dnudR)*deriv.dlambdadR + 
        (srcHess.dlambdadnu*deriv.dlambdadR + srcHess.dnu2*deriv.dnudR)*deriv.dnudR + 
        srcGrad.dlambda*deriv2.d2lambdadR2 + srcGrad.dnu*deriv2.d2nudR2;
    dest.dRdz = 
        (srcHess.dlambda2*deriv.dlambdadz + srcHess.dlambdadnu*deriv.dnudz)*deriv.dlambdadR + 
        (srcHess.dlambdadnu*deriv.dlambdadz + srcHess.dnu2*deriv.dnudz)*deriv.dnudR + 
        srcGrad.dlambda*deriv2.d2lambdadRdz + srcGrad.dnu*deriv2.d2nudRdz;
    dest.dz2 = 
        (srcHess.dlambda2*deriv.dlambdadz + srcHess.dlambdadnu*deriv.dnudz)*deriv.dlambdadz + 
        (srcHess.dlambdadnu*deriv.dlambdadz + srcHess.dnu2*deriv.dnudz)*deriv.dnudz + 
        srcGrad.dlambda*deriv2.d2lambdadz2 + srcGrad.dnu*deriv2.d2nudz2;
    dest.dRdphi = dest.dzdphi = dest.dphi2 = 0;  // assuming no dependence on phi
    return dest;
}

//------ conversion of derivatives of f(r) into gradients/hessians in different coord.sys. ------//
template<>
EXP void evalAndConvertSph(const math::IFunction& F,
    const PosCar& pos, double* value, GradCar* deriv, HessCar* deriv2)
{
    assert(F.numDerivs()>=2);
    const double r=sqrt(pow_2(pos.x)+pow_2(pos.y)+pow_2(pos.z));
    if(deriv==NULL && deriv2==NULL) {
        F.evalDeriv(r, value);
        return;
    }
    double der, der2;
    F.evalDeriv(r, value, &der, deriv2!=NULL ? &der2 : 0);
    double x_over_r=pos.x/r, y_over_r=pos.y/r, z_over_r=pos.z/r;
    if(r==0) {
        x_over_r=y_over_r=z_over_r=0;
    }
    if(deriv) {
        deriv->dx = x_over_r*der;
        deriv->dy = y_over_r*der;
        deriv->dz = z_over_r*der;
    }
    if(deriv2) {
        double der_over_r=der/r, dd=der2-der_over_r;
        if(r==0) { 
            dd=0;
            if(der==0) der_over_r=der2;
        }
        deriv2->dx2 = pow_2(x_over_r)*dd + der_over_r;
        deriv2->dy2 = pow_2(y_over_r)*dd + der_over_r;
        deriv2->dz2 = pow_2(z_over_r)*dd + der_over_r;
        deriv2->dxdy= x_over_r*y_over_r*dd;
        deriv2->dydz= y_over_r*z_over_r*dd;
        deriv2->dxdz= x_over_r*z_over_r*dd;
    }
}

template<>
EXP void evalAndConvertSph(const math::IFunction& F,
    const PosCyl& pos, double* value, GradCyl* deriv, HessCyl* deriv2)
{
    assert(F.numDerivs()>=2);
    const double r=sqrt(pow_2(pos.R)+pow_2(pos.z));
    if(deriv==NULL && deriv2==NULL) {
        F.evalDeriv(r, value);
        return;
    }
    double der, der2;
    F.evalDeriv(r, value, &der, deriv2!=NULL ? &der2 : 0);
    double R_over_r=pos.R/r, z_over_r=pos.z/r;
    if(r==0) {
        R_over_r=z_over_r=0;
    }
    if(deriv) {
        deriv->dR = R_over_r*der;
        deriv->dz = z_over_r*der;
        deriv->dphi = 0;
    }
    if(deriv2) {
        double der_over_r=der/r, dd=der2-der_over_r;
        if(r==0) {
            dd=0;
            if(der==0) der_over_r=der2;
        }
        deriv2->dR2 = pow_2(R_over_r)*dd + der_over_r;
        deriv2->dz2 = pow_2(z_over_r)*dd + der_over_r;
        deriv2->dRdz= R_over_r*z_over_r*dd;
        deriv2->dRdphi=deriv2->dzdphi=deriv2->dphi2=0;
    }
}

template<>
EXP void evalAndConvertSph(const math::IFunction& F,
    const PosSph& pos, double* value, GradSph* deriv, HessSph* deriv2)
{
    assert(F.numDerivs()>=2);
    double der, der2;
    F.evalDeriv(pos.r, value, deriv!=NULL ? &der : NULL, deriv2!=NULL ? &der2 : NULL);
    if(deriv) {
        deriv->dr = der;
        deriv->dtheta = deriv->dphi = 0;
    }
    if(deriv2) {
        deriv2->dr2 = der2;
        deriv2->dtheta2 = deriv2->dphi2 = deriv2->drdtheta = deriv2->drdphi = deriv2->dthetadphi = 0;
    }
}


//------ 3x3 matrix representing a [passive] rotation specified by Euler angles ------//

EXP void makeRotationMatrix(double alpha, double beta, double gamma, double mat[9])
{
	double sa, ca, sb, cb, sc, cc;
	math::sincos(alpha, sa, ca);
	math::sincos(beta,  sb, cb);
	math::sincos(gamma, sc, cc);
	mat[0] =  ca * cc - sa * cb * sc;
	mat[1] =  sa * cc + ca * cb * sc;
	mat[2] =  sb * sc;
	mat[3] = -ca * sc - sa * cb * cc;
	mat[4] = -sa * sc + ca * cb * cc;
	mat[5] =  sb * cc;
	mat[6] =  sa * sb;
	mat[7] = -ca * sb;
	mat[8] =  cb;
}

EXP std::vector<double> makeRotationMatrix(double alpha, double beta, double gamma)
{
	std::vector<double> mat(9);
	double sa, ca, sb, cb, sc, cc;
	math::sincos(alpha, sa, ca);
	math::sincos(beta,  sb, cb);
	math::sincos(gamma, sc, cc);
	mat[0] =  ca*cc-sc*sb*sa;
	mat[1] = -sa*cc-sc*cb*ca;
	mat[2] =  sc*sb;
	mat[3] =  ca*sc+cc*cb*sa;
	mat[4] = -sa*sc+cc*cb*ca;
	mat[5] = -sb*cc;
	mat[6] =  sa*sb;
	mat[7] =  sb*ca;
	mat[8] =  cb;
//	printf("%f\t%f\t%f\n%f\t%f\t%f\n%f\t%f\t%f\n",mat[0],mat[1],mat[2],
//	       mat[3],mat[4],mat[5],mat[6],mat[7],mat[8]);
	return mat;
}

}  // namespace coord
