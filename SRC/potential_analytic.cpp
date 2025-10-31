#include "potential_analytic.h"
#include <cmath>

namespace potential{

void Plummer::evalDeriv(double r,
    double* potential, double* deriv, double* deriv2) const
{
	double invrsq = mass?  1. / (pow_2(r) + pow_2(scaleRadius)) : 0;  // if mass=0, output 0
	double pot = -mass * sqrt(invrsq);
    if(potential)
        *potential = pot;
    if(deriv)
        *deriv = -pot * r * invrsq;
    if(deriv2)
        *deriv2 = pot * (2 * pow_2(r * invrsq) - pow_2(scaleRadius * invrsq));
}
double Plummer::Phi(double r) const{
	double invrsq = mass?  1. / (pow_2(r) + pow_2(scaleRadius)) : 0;  // if mass=0, output 0
	return -mass * sqrt(invrsq);
}
double Plummer::dPhidr(double r) const{
	double invrsq = mass?  1. / (pow_2(r) + pow_2(scaleRadius)) : 0;  // if mass=0, output 0
	return mass * r * pow(invrsq,1.5);
}
double Plummer::densitySph(const coord::PosSph &pos) const
{
	double invrsq = 1. / (pow_2(pos.r) + pow_2(scaleRadius));
	return 0.75/M_PI * mass * pow_2(scaleRadius * invrsq) * sqrt(invrsq);
}

double Plummer::enclosedMass(double r) const
{
    if(scaleRadius==0)
        return mass;
    return mass / pow_3(sqrt(pow_2(scaleRadius/r) + 1));
}

void Isochrone::evalDeriv(double r,
    double* potential, double* deriv, double* deriv2) const
{
    double rb  = sqrt(pow_2(r) + pow_2(scaleRadius));
    double brb = scaleRadius + rb;
    double pot = -mass / brb;
    if(potential)
        *potential = pot;
    if(deriv)
        *deriv = -pot * r / (rb * brb);
    if(deriv2)
        *deriv2 = pot * (2*pow_2(r / (rb * brb)) - pow_2(scaleRadius / (rb * brb)) * (1 + scaleRadius / rb));
}
double Isochrone::density(double r) const{
	double a=sqrt(pow_2(scaleRadius) + r*r);
	return mass/(4*M_PI)*(3*(scaleRadius + a)*a*a-r*r*(scaleRadius + 3*a))/
			(pow_3((scaleRadius + a)*a));
}
double Isochrone::Phi(double r) const{
	double a=sqrt(pow_2(scaleRadius) + r*r);
	return -mass/(scaleRadius+a);
}
double Isochrone::dPhidr(double r) const{
	double a=sqrt(pow_2(scaleRadius) + r*r);
	return mass*r/(a*pow_2(scaleRadius + r));
}
double Isochrone::Mass(double r) const{
	double a=sqrt(pow_2(scaleRadius)+r*r);
	return mass*pow_3(r)/(pow_2(scaleRadius + a)*a);
}			

void Hernquist::evalDeriv(double r,
			  double* potential, double* deriv, double* deriv2) const{
	if(potential)
 		*potential = Phi(r);
	if(deriv)
		*deriv = dPhidr(r);
	if(deriv2)
		*deriv2 = -2*mass/pow_3(scaleRadius + r);
}
double Hernquist::density(double r) const{
	double x=r/scaleRadius;
	return rho0/(x*pow_3(1+x));
}
double Hernquist::Phi(double r) const{
	return -mass/(scaleRadius + r);
}
double Hernquist::dPhidr(double r) const{
	return mass/pow_2(scaleRadius + r);
}
double Hernquist::Mass(double r) const{
	double x=r/scaleRadius;
	return 2*mass*pow_2(x/(1+x));
}			

void NFW::evalDeriv(double r,
		    double* potential, double* deriv, double* deriv2) const
{
	double rrel = r / scaleRadius;
	double ln_over_r = r==INFINITY ? 0 :
			   rrel > 1e-3 ? log(1 + rrel) / r :
	// accurate asymptotic expansion at r->0
			   (1 - 0.5 * rrel * (1 - 2./3 * rrel * (1 - 3./4 * rrel))) / scaleRadius;
	if(potential)
		*potential = -mass * ln_over_r;
	if(deriv)
		*deriv = mass * (r==0 ? 0.5 / pow_2(scaleRadius) :
				 (ln_over_r - 1/(r+scaleRadius)) / r );
	if(deriv2)
		*deriv2 = -mass * (r==0 ? 2./3 / pow_3(scaleRadius) :
				   (2*ln_over_r - (2*scaleRadius + 3*r) / pow_2(scaleRadius+r) ) / pow_2(r) );
}

void purePower::evalDeriv(double r,
		    double* potential, double* deriv, double* deriv2) const
{
	if(potential)
		if(index==0) *potential = log(r)-100;
		else *potential = index<0? -pow(r,index) : pow(r,index)-100;// Phi always increases outwards
	if(deriv)
		if(index==0) *deriv = 1/r;
		else *deriv = index<0? -index*pow(r,index-1) : index*pow(r,index-1);
	if(deriv2)
		if(index==0) *deriv2 = - 1/(r*r);
		else *deriv2 = index<0? -index*(index-1)*pow(r,index-2) : index*(index-1)*pow(r,index-2);
}

void MiyamotoNagai::evalCyl(const coord::PosCyl &pos,
    double* potential, coord::GradCyl* deriv, coord::HessCyl* deriv2) const
{
    double zb    = sqrt(pow_2(pos.z) + pow_2(scaleRadiusB));
    double azb2  = pow_2(scaleRadiusA + zb);
    double den2  = 1. / (pow_2(pos.R) + azb2);
    double denom = sqrt(den2);
    double Rsc   = pos.R * denom;
    double zsc   = pos.z * denom;
    if(potential)
        *potential = -mass * denom;
    if(deriv) {
        deriv->dR  = mass * den2 * Rsc;
        deriv->dz  = mass * den2 * zsc * (1 + scaleRadiusA/zb);
        deriv->dphi= 0;
    }
    if(deriv2) {
        double mden3 = mass * denom * den2;
        deriv2->dR2  = mden3 * (azb2*den2 - 2*pow_2(Rsc));
        deriv2->dz2  = mden3 * ( (pow_2(Rsc) - 2*azb2*den2) * pow_2(pos.z/zb) +
            pow_2(scaleRadiusB) * (scaleRadiusA/zb + 1) * (pow_2(Rsc) + azb2*den2) / pow_2(zb) );
        deriv2->dRdz = mden3 * -3 * Rsc * zsc * (scaleRadiusA/zb + 1);
        deriv2->dRdphi = deriv2->dzdphi = deriv2->dphi2 = 0;
    }
}

void Logarithmic::evalCar(const coord::PosCar &pos,
			  double* potential, coord::GradCar* deriv, coord::HessCar* deriv2) const
{
	double m2 = 1 + (pow_2(pos.x) + pow_2(pos.y)/p2 + pow_2(pos.z)/q2)/Rc2;
	double L = log(m2), oLpLm = 1/(1+L/Lm), oLpLmsq = pow_2(oLpLm);
	if(potential)
		*potential = -0.5 * sigma2 * Lm * oLpLm;
	if(deriv) {
		deriv->dx = pos.x * sigma2/(Rc2*m2)    * oLpLmsq;
		deriv->dy = pos.y * sigma2/(Rc2*m2*p2) * oLpLmsq;
		deriv->dz = pos.z * sigma2/(Rc2*m2*q2) * oLpLmsq;
	}
	if(deriv2) {
		double dLdx = 2*pos.x/(m2*Rc2), dLdy = 2*pos.y/(p2*m2*Rc2), dLdz = 2*pos.z/(q2*m2*Rc2);
		deriv2->dx2 = sigma2 * oLpLmsq * ((1/(m2*Rc2)    - 2 * pow_2(pos.x / (m2*Rc2)))
					- dLdx * dLdx /Lm * oLpLm);
		deriv2->dy2 = sigma2 * oLpLmsq * ((1/(Rc2*m2*p2) - 2 * pow_2(pos.y / (Rc2*m2*p2)))
					- dLdy * dLdy /Lm *oLpLm);
		deriv2->dz2 = sigma2 * oLpLmsq * ((1/(Rc2*m2*q2) - 2 * pow_2(pos.z / (Rc2*m2*q2)))
					- dLdz * dLdz /Lm *oLpLm);
		deriv2->dxdy=-sigma2 * oLpLmsq * (pos.x * pos.y * 2 / (pow_2(Rc2*m2) * p2)
					+ dLdx * dLdy /Lm * oLpLm);
		deriv2->dydz=-sigma2 * oLpLmsq * (pos.y * pos.z * 2 / (pow_2(Rc2*m2) * p2 * q2)
					+ dLdy * dLdz /Lm * oLpLm);
		deriv2->dxdz=-sigma2 * oLpLmsq * (pos.z * pos.x * 2 / (pow_2(Rc2*m2) * q2)
					+ dLdx * dLdz /Lm * oLpLm);
	}}

void logRe::evalCyl(const coord::PosCyl &pos,
			  double* potential, coord::GradCyl* deriv, coord::HessCyl* deriv2) const
{
	double x2=pow_2(pos.R), z2=pow_2(pos.z);
	double m2 = coreRadius2 + x2 + z2/q2, r=sqrt(x2+z2);
	if(ifRe) m2-=.5*r/Re*(x2-z2);
	if(potential)
		*potential = 0.5 * sigma2 * log(m2);
	if(deriv) {
		deriv->dR = pos.R * sigma2/m2;
		deriv->dz = pos.z * sigma2/m2/q2;
		deriv->dphi = 0;
		if(ifRe){
			deriv->dR-=(.5*pos.R*(x2-z2)/r+r*pos.R)*sigma2/(2*Re*m2);
			deriv->dz-=(.5*pos.z*(x2-z2)/r-r*pos.z)*sigma2/(2*Re*m2);
		}
	}
	if(deriv2) {
		printf("Hessian not implemented\n");
	}
}

void Harmonic::evalCar(const coord::PosCar &pos,
		       double* potential, coord::GradCar* deriv, coord::HessCar* deriv2) const
{
	if(potential)
		*potential = 0.5*Omega2 * (pow_2(pos.x) + pow_2(pos.y)/p2 + pow_2(pos.z)/q2);
	if(deriv) {
		deriv->dx = pos.x*Omega2;
		deriv->dy = pos.y*Omega2/p2;
		deriv->dz = pos.z*Omega2/q2;
	}
	if(deriv2) {
		deriv2->dx2 = Omega2;
		deriv2->dy2 = Omega2/p2;
		deriv2->dz2 = Omega2/q2;
		deriv2->dxdy=deriv2->dydz=deriv2->dxdz=0;
	}
}

void HenonHeiles::evalCyl(const coord::PosCyl &pos,
		       double* potential, coord::GradCyl* deriv, coord::HessCyl* deriv2) const
{
	if(potential)
		*potential = 0.5* (pow_2(pos.R) + pow_2(pos.z))
			     + pow_2(pos.z)*pos.R - pow_3(pos.R)/3;
	if(deriv) {
		deriv->dR = pos.R + pow_2(pos.z) - pow_2(pos.R);
		deriv->dz = pos.z + 2*pos.z*pos.R;
		deriv->dphi = 0;
	}
	if(deriv2) {
		deriv2->dR2 = 1 - 2*pos.R;
		deriv2->dz2 = 1 + 2*pos.R;
		deriv2->dphi2 = 0;
		deriv2->dRdz = 2*pos.z;
		deriv2->dzdphi=deriv2->dRdphi=0;
	}
}

void Sormani::evalCyl(const coord::PosCyl &pos,
		       double* potential, coord::GradCyl* deriv, coord::HessCyl* deriv2) const
{
	double R2 = pow_2(pos.R), z2 = pow_2(pos.z)/q2;
	double m2 = Rb2 + R2 + z2;
	double Phi2 = K*R2*pow(m2,-2.5);
	double cos2phi = cos(2*pos.phi);
	if(potential)
		*potential = -Phi2 * cos2phi;
	if(deriv) {
		deriv->dR = -K*pos.R*pow(m2,-2.5)*(2-5*R2/m2) * cos2phi;
		deriv->dz = 5*Phi2/m2*pos.z/q2 * cos2phi;
		deriv->dphi = Phi2 * 2*sin(2*pos.phi);
	}
	if(deriv2) {
		deriv2->dR2 = -K*pow(m2,-2.5)*(2-25*R2/m2+35*pow_2(R2/m2)) * cos2phi;
		deriv2->dz2 = K*R2/q2*pow(m2,-3.5)*(5-35/q2*pow_2(pos.z)/m2) * cos2phi;
		deriv2->dphi2 = Phi2 * 4*cos2phi;
		deriv2->dRdz = K*pos.R*pos.z/q2*pow(m2,-3.5)*(10-35*R2/m2) * cos2phi;
		deriv2->dRdphi = K*pos.R*pow(m2,-2.5)*(2-5*R2/m2) * 2*sin(2*pos.phi);
		deriv2->dzdphi = -5*Phi2/m2*pos.z/q2 * 2*sin(2*pos.phi);
	}
}

}  // namespace potential
