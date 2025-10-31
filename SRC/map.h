/* Stuff to do with toy maps, ie combinatons of an HJ map and a point
 * map
 */
#pragma once
#include "actions_base.h"
#include "actions_harmonicoscillator.h"
#include "actions_newisochrone.h"
#include "coord.h"
#include "math_core.h"
#include "math_linalg.h"

namespace actions{
/* Class of maps from spherical polar coords to cylindrical polar
 * coords using confocal ellipsoidal coords. Used with isochrone HJ map
 */
class EXP PTIso: public BasePointTransform{
	double t2v(const double theta, double* dvdt = NULL, double* d2vdt2 = NULL, 
		   double *dvdP = NULL, double* ddvdtdP = NULL) const;
	double v2t(const double, double* = NULL) const;
	public:
		coord::UVSph cs;
		std::vector<double> paramsF;
		std::vector<double> paramsFr;
		math::ScalingInfTh sc;
		int N, Nr;
		PTIso():
		    paramsF({}),paramsFr({}),N(0),Nr(0){}
		PTIso(double _D):
		    cs(_D),paramsF({}),paramsFr({}),N(0),Nr(0){}
		PTIso(coord::UVSph _cs):
		    cs(_cs),paramsF({}),paramsFr({}),N(0),Nr(0){}
		PTIso(double _D, const std::vector<double>& _p) :
		    cs(_D), paramsF(_p), paramsFr({}), N(_p.size()),Nr(0) {}
		PTIso(coord::UVSph _cs, const std::vector<double>& _p) :
		    cs(_cs), paramsF(_p), paramsFr({}), N(_p.size()),Nr(0) {}
		PTIso(double _D, math::ScalingInfTh _sc, std::vector<double>& _p, std::vector<double> _pr) :
		    cs(_D), paramsF(_p), paramsFr(_pr), sc(_sc), N(_p.size()),Nr(_pr.size()) {}
		PTIso(coord::UVSph _cs, math::ScalingInfTh _sc, std::vector<double>& _p, std::vector<double> _pr) :
		    cs(_cs), paramsF(_p), paramsFr(_pr), sc(_sc), N(_p.size()), Nr(_pr.size()) {}
		virtual coord::PosMomCyl map(const coord::PosMomSph &point, coord::PosMomCyl* dRzdP,
					     coord::PosMomCyl *dRzdFr, coord::PosMomCyl *dRzdFz) const;
		virtual coord::PosMomCyl map(const coord::PosMomSph &point) const{
			return map(point, NULL, NULL, NULL);
		}
		virtual coord::PosMomCyl map(const coord::PosMomCyl &point) const{
			return map(coord::toPosMomSph(point));
		}
		virtual coord::PosMomCyl map(const coord::PosMomSph &point, math::Matrix<double>& dRzdrt,
					     coord::PosMomCyl* dRzdP, coord::PosMomCyl *dRzdFr,coord::PosMomCyl *dRzdFz) const;
		virtual coord::PosMomCyl map(const coord::PosMomSph &point, math::Matrix<double>& dRzdrt) const{
			return map(point,dRzdrt,NULL,NULL,NULL);
		}
		virtual coord::PosMomCyl map(const coord::PosMomCyl &point, math::Matrix<double>& dRzdrt) const{
			return map(coord::toPosMomSph(point));
		}
		virtual coord::PosMomSph revmapSph(const coord::PosMomCyl &point) const;
		virtual coord::PosMomCyl revmap(const coord::PosMomCyl &point) const{
			return coord::toPosMomCyl(revmapSph(point));
		}
		virtual const char* name() const{return "Isochrone";}
		virtual void getParams(double* params=NULL, double* Fourr=NULL, double* Fourz=NULL) const{
			if(params){
				if(N==0&&Nr==0){
					params[0]=cs.Delta;
				}else{
					params[0]=cs.Delta;
					params[1]=sc.x0;
				}
			}
			if(Fourr){
				for(int i=0;i<Nr;i++)Fourr[i]=paramsFr[i];
			}
			if(Fourz){
				for(int i=0;i<N;i++)Fourz[i]=paramsF[i];
			}
		}
		virtual int numParams() const{
			if(N==0&&Nr==0)return 1;
			return 2;
		}
		virtual int FourierSizer() const{return Nr;}
		virtual int FourierSizez() const{return N;}
};
/* Class of maps between cylindrical polars for use with HO HJ map 
 */
class EXP PTHarm : public BasePointTransform{
	double zntov(const double zn, double* dvdt = NULL, double* d2vdt2 = NULL, 
		     double *dvdz0=NULL, double* ddvdzdz0=NULL, double *dvdFC = NULL, double *ddvdtdFC = NULL) const;
	double vtozn(const double v, double* dzdv = NULL) const;
	public:
		coord::UVSph cs;
		std::vector<double> paramsF;
		std::vector<double> paramsFr;
		math::ScalingInfTh sc;
		math::ScalingInfTh scz;
		int N, Nr;
		PTHarm():
		    paramsF({}),paramsFr({}),N(0),Nr(0){}
		PTHarm(double _D, math::ScalingInfTh _sc, math::ScalingInfTh _scz, std::vector<double>& _p, std::vector<double> _pr) :
		    cs(_D), paramsF(_p), paramsFr(_pr), sc(_sc), scz(_scz), N(_p.size()),Nr(_pr.size()) {
		}
		PTHarm(coord::UVSph _cs, math::ScalingInfTh _sc, math::ScalingInfTh _scz, std::vector<double>& _p, std::vector<double> _pr) :
		    cs(_cs), paramsF(_p), paramsFr(_pr), sc(_sc), scz(_scz), N(_p.size()), Nr(_pr.size()) {
		}
		virtual coord::PosMomCyl map(const coord::PosMomCyl &point, coord::PosMomCyl* dRzdP,
					     coord::PosMomCyl *dRzdFr, coord::PosMomCyl *dRzdFz) const;
		virtual coord::PosMomCyl map(const coord::PosMomCyl &point) const{
			return map(point,NULL,NULL,NULL);
		}
		virtual coord::PosMomCyl map(const coord::PosMomSph &point) const{
			return map(coord::toPosMomCyl(point));
		}
		virtual coord::PosMomCyl map(const coord::PosMomCyl &point, math::Matrix<double>& dRzdrt,
					     coord::PosMomCyl* dRzdP, coord::PosMomCyl *dRzdFr,coord::PosMomCyl *dRzdFz) const;
		virtual coord::PosMomCyl map(const coord::PosMomCyl &point, math::Matrix<double>& dRzdrt) const{
			return map(point, dRzdrt, NULL, NULL, NULL);
		}
		virtual coord::PosMomCyl map(const coord::PosMomSph &point, math::Matrix<double>& dRzdrt) const{
			return map(coord::toPosMomCyl(point),dRzdrt);
		}
		virtual coord::PosMomCyl revmap(const coord::PosMomCyl &point) const;
		virtual coord::PosMomSph revmapSph(const coord::PosMomCyl &point) const{
			return coord::toPosMomSph(revmap(point));
		}
		virtual void getParams(double* params=NULL, double* Fourr=NULL, double* Fourz=NULL) const{
			if(params){
				params[0]=cs.Delta;
				params[1]=sc.x0;
				params[2]=scz.x0;
			}
			if(Fourr){
				for(int i=0;i<Nr;i++)Fourr[i]=paramsFr[i];
			}
			if(Fourz){
				for(int i=0;i<N;i++)Fourz[i]=paramsF[i];
			}
		}
		virtual const char* name() const{return "Harmonic Oscillator";}
		virtual std::vector<double> FourSerr() const{
			return paramsFr;
		}
		virtual std::vector<double> FourSerz() const{
			return paramsF;
		}
		virtual int numParams() const{return 3;}
		virtual int FourierSizer() const{return Nr;}
		virtual int FourierSizez() const{return N;}
};

    EXP PTIso interpPTIso(double x, const PTIso&, const PTIso&);
    EXP PTHarm interpPTHarm(double x, const PTHarm&, const PTHarm&);
    EXP PtrPointTransform interpPointTrans(double x, const PtrPointTransform, const PtrPointTransform);
    
/* Class of toy maps that combines point trans with isochrone HJ map
 */
    class EXP ToyMapIso:public BaseToyMap{
	    private:
		    PtrPointTransform PtrPT;
	    public:
		    Isochrone Is;
		    PTIso PT;
		    ToyMapIso(){}
		    ToyMapIso(Isochrone _Is, PTIso _PT) : Is(_Is), PT(_PT){
			//PtrPT=std::make_shared<const BasePointTransform>(_PT);
		    };
		    virtual Actions pq2J(const coord::PosMomCyl Rzp) const {
			    return Is.pq2J(PT.revmapSph(Rzp));
		    }
		    virtual ActionAngles pq2aa(const coord::PosMomCyl& Rz) const {
			    coord::PosMomSph rp(PT.revmapSph(Rz));
			    if (Is.H(rp) < 0) return Is.pq2aa(rp);
			    printf("ToyMap::pq2aa: Iso H = %g >=0!",Is.H(rp));
			    exit(0);
		    }
		    virtual coord::PosMomCyl from_aaT(const ActionAngles& aaT) const {
			    return PT.map(Is.aa2pq(aaT));
		    }
		    virtual coord::PosMomCyl from_aaT(const Actions& J, const Angles& theta) const {
			    return from_aaT(ActionAngles(J, theta));
		    }
		    virtual coord::PosMomCyl from_aaT(const ActionAngles& aaT, coord::PosMomCyl* dRzdPpot,
			    coord::PosMomCyl* dRzdPPT, coord::PosMomCyl* dRzdFr, coord::PosMomCyl* dRzdFz) const;
		    virtual coord::PosMomCyl from_aaT(const Actions& J, Angles& theta) const {
			    return from_aaT(ActionAngles(J, theta));
		    }
		    virtual coord::PosMomCyl from_aaT(const ActionAngles& aaT, DerivAct<coord::Cyl>& dRzdJ) const;
		    virtual coord::PosMomCyl from_aaT(const Actions& J, const Angles& theta,
			    DerivAct<coord::Cyl>& dRzdJ) const {
			    return from_aaT(ActionAngles(J, theta), dRzdJ);
		    }
		    virtual coord::PosMomCyl from_aaT(const ActionAngles& aaT, DerivAng<coord::Cyl>& dRzdT) const;
		    virtual coord::PosMomCyl from_aaT(const Actions& J, const Angles& theta,
			    DerivAng<coord::Cyl>& dRzdt) const {
			    return from_aaT(ActionAngles(J, theta), dRzdt);
		    }
		    virtual PtrPointTransform getPointTrans() const{return PtrPT;}
		    virtual ToyPotType getToyMapType() const{return ToyPotType::Is;}
		    virtual int PTParamSize() const{return PT.N+PT.Nr+2;}
		    virtual int potParamSize() const{return 2;}
		    virtual const char* name() const { return myName(); }
		    static const char* myName() { static const char* text = "Isochrone"; return text; }
		    void getParams(double* Js=NULL,double* b=NULL) const{
			    if(Js)*Js=Is.Js;
			    if(b)*b=Is.b;
		    };
    };

    /* Class of toy maps that combine a point transformation with the
     * HO HJ map
     */
    class EXP ToyMapHarm:public BaseToyMap{
	private:
		PtrPointTransform PtrPT;
	public:
		HarmonicOscilattor HOs;
		PTHarm PT;
		ToyMapHarm() {}
		ToyMapHarm(HarmonicOscilattor _HOs, PTHarm _PT) : HOs(_HOs),PT(_PT){
			//PtrPT=std::make_shared<const BasePointTransform>(_PT);
		}
		Actions pq2J(const coord::PosMomCyl Rzp) const {
			return HOs.pq2J(PT.revmap(Rzp));
		}
		ActionAngles pq2aa(const coord::PosMomCyl& Rz) const {
			coord::PosMomCyl rp(PT.revmap(Rz));
			return HOs.pq2aa(rp);
		}
		coord::PosMomCyl from_aaT(const ActionAngles& aaT) const {
			return PT.map(HOs.aa2pq(aaT));
		}
		coord::PosMomCyl from_aaT(const Actions& J, const Angles& theta) const {
			return from_aaT(ActionAngles(J, theta));
		}
		coord::PosMomCyl from_aaT(const ActionAngles& aaT, coord::PosMomCyl* dRzdPpot,
			coord::PosMomCyl* dRzdPPT, coord::PosMomCyl* dRzdFr, coord::PosMomCyl* dRzdFz) const;
		coord::PosMomCyl from_aaT(const ActionAngles& aaT, DerivAct<coord::Cyl>& dRzdJ) const;
		coord::PosMomCyl from_aaT(const Actions& J, const Angles& theta,
			DerivAct<coord::Cyl>& dRzdJ) const {
			return from_aaT(ActionAngles(J, theta), dRzdJ);
		}
		coord::PosMomCyl from_aaT(const ActionAngles& aaT, DerivAng<coord::Cyl>& dRzdT) const;
		coord::PosMomCyl from_aaT(const Actions& J, const Angles& theta,
			DerivAng<coord::Cyl>& dRzdt) const {
			return from_aaT(ActionAngles(J, theta), dRzdt);
		}
		virtual PtrPointTransform getPointTrans() const{return PtrPT;}
		virtual ToyPotType getToyMapType() const{return ToyPotType::HO;}
		virtual int PTParamSize() const{return PT.N+PT.Nr+3;}
		virtual int potParamSize() const{return 2;}
		virtual const char* name() const { return myName(); }
		static const char* myName() {
			static const char* text = "Harmonic Oscillator"; return text;
		}
		void getParams(double* omegar=NULL,double* omegaz=NULL) const{
			if(omegar)*omegar=HOs.omegaR;
			if(omegaz)*omegaz=HOs.omegaR;
		}
	};
    EXP PtrToyMap interpPtrToyMap(double, const PtrToyMap&, const PtrToyMap&);
    /*
	  Map of equatorial plane of (r,theta,phi) to system in (x,y) plane
	  extended along y
	*/
	class PsiCoord {
	private:
		const double a, b;
	public:
		PsiCoord(const double _a, const double _b) : a(_a), b(_b) {}
		double Psi(const double phi, double& dpsidphi) const {
			dpsidphi = 1 + 2 * a * cos(2 * phi) + 4 * b * cos(4 * phi);
			return phi + a * sin(2 * phi) + b * sin(4 * phi);
		}
		double Phi(const double psi, double& dpsidphi) const {
			double phi = psi, D = Psi(phi, dpsidphi) - psi;
			while (fabs(D) > 1e-5) {
				phi -= D / dpsidphi;
				D = Psi(phi, dpsidphi) - psi;
			}
			return phi;
		}
	};
}