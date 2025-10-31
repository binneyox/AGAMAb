/* Stuff to do with toy maps, ie combinatons of an HJ map and a point
 * map
 */
#include "map.h"
namespace actions{
    namespace{
        //function used to invert v(theta) to get theta(v)
        class vthet : public math::IFunction {
        private:
            std::vector<double> p;
            double v;
        public:
            vthet(std::vector<double> _p, double _v) :p(_p),v(_v) {};
            virtual unsigned int numDerivs()const { return 1; }
            virtual void evalDeriv(const double x,
                double* value, double* deriv, double* deriv2) const {
                double f1 = x - v;
                double f2 = 1.;
                for (int i = 0;i < p.size();i++) {
                    if (value)f1 += sin(2 * (i + 1) * x) * p[i];
                    if (deriv)f2 += 2 * (i + 1) * cos(2 * (i + 1) * x) * p[i];
                }
                if (value)*value = f1;
                if (deriv)*deriv = f2;
            }
        };
        class xfind : public math::IFunction {
        private:
            std::vector<double> p;
            double x;
        public:
            xfind(std::vector<double> _p, double _x) :p(_p), x(_x) {};
            virtual unsigned int numDerivs()const { return 1; }
            virtual void evalDeriv(const double x1,
                double* value, double* deriv, double* deriv2) const {
                if (x1 == 1.) {
                    if (value)*value=1.;
                    if (deriv)*deriv = -1.;
                    return;
                }
                double f1 = atanh(2 * x1 - 1) - x;
                double f2 = 2 / (1 - pow_2(2 * x1 - 1));
                for (int i = 0;i < p.size();i++) {
                    if (value)f1 += sin(2 * (i + 1) * M_PI * x1) * p[i];
                    if (deriv)f2 += 2 * (i + 1) * M_PI * cos(2 * (i + 1) * M_PI * x1) * p[i];
                }
                if (value)*value = f1;
                if (deriv)*deriv = f2;
            }
        };
        class zfind : public math::IFunction {
        private:
            std::vector<double> p;
            double v0;
        public:
            zfind(std::vector<double> _p, double _v0) :p(_p),v0(_v0)  {};
            virtual unsigned int numDerivs()const { return 1; }
            virtual void evalDeriv(const double x1,
                double* value, double* deriv, double* deriv2) const {
                double dvdsz = -1.;
                double deltav = M_PI - x1-v0;
                for (int i = 0; i < p.size(); i++) {
                    deltav += p[i] * sin(2 * (i + 1) * x1);
                    dvdsz += 2 * (i + 1) * p[i] * cos(2 * (i + 1) * x1);
                }
                if (value)*value = deltav;
                if (deriv)*deriv = dvdsz;
            }
        };
        double FourSer(double x,std::vector<double> values, double *dSdx=NULL, double *d2Sdx2=NULL,
         double *dSdP=NULL, double *ddSdxdP=NULL){
            double S=0,dS=0,d2S=0;
            for(int i=0;i<values.size();i++){
                S+=values[i] * sin(2 * (i + 1) * x);
                if(dSdx)dS += 2 * (i + 1) * values[i] * cos(2 * (i + 1) * x);
			    if(d2Sdx2)d2S += -pow_2(2 * (i + 1)) * values[i] * sin(2 * (i + 1) * x);
                if(dSdP)dSdP[i]=sin(2*(i+1)*x);
                if(ddSdxdP)ddSdxdP[i]=2*(i+1)*cos(2*(i+1)*x);
            }
            if(dSdx)*dSdx=dS;
            if(d2Sdx2)*d2Sdx2=d2S;
            return S;
        }
        double r2rn(const double r, const std::vector<double> values, const math::ScalingInfTh &sc, 
            double* drndr=NULL, double* d2rndr2=NULL,double *drndx0=NULL,double *ddrndrdx0=NULL, 
            double *drndFC=NULL, double *ddrndrdFC=NULL){
            double dsdr, d2sdr2;
            double sr = M_PI*math::scale(sc, r, &dsdr, &d2sdr2);
            dsdr*=M_PI; d2sdr2*=M_PI;
            double drFds = 0, d2rFds2 = 0;
            bool needDeriv=drndr||d2rndr2||drndx0||ddrndrdx0;
            double rn=r+sc.x0*FourSer(sr,values,needDeriv?&drFds:NULL,d2rndr2||ddrndrdx0?&d2rFds2:NULL,
                drndFC,ddrndrdFC);
            int sgn=math::sign(rn);
            drFds*=sc.x0;
            d2rFds2*=sc.x0;
            if(ddrndrdFC){
                for(int i=0;i<values.size();i++){
                    drndFC[i]*=sc.x0*sgn;
                    ddrndrdFC[i]*=dsdr*sgn*sc.x0;
                }
            }
            if (drndr)*drndr = (1 + drFds * dsdr)*sgn;
            if (d2rndr2) {
                double d2rnFdr2 = pow_2(dsdr) * d2rFds2 + drFds * d2sdr2;
                *d2rndr2 = d2rnFdr2*sgn;
            }
            if(drndx0||ddrndrdx0){
                double dsdr0=dsdr*(-r/sc.x0);
                double ddsrdrdr0=d2sdr2*(-r/sc.x0)-dsdr/sc.x0;
                if(drndx0)*drndx0=sgn*drFds*dsdr0;
                if(ddrndrdx0)*ddrndrdx0=(d2rFds2*dsdr0*dsdr+drFds*ddsrdrdr0)*sgn;
            }
            return fabs(rn);
        }
        double rn2r(const double rn,const std::vector<double> values,
            const math::ScalingInfTh sc, double* drdrn=NULL){
            if(values.size()==0){
                if(drdrn)*drdrn=1;
                return rn;
            }
            double v = rn/sc.x0;
            xfind x1(values, v);
            double t = math::findRoot(x1, 0.5, 1., 1e-9);
            double r1 = math::unscale(sc, t);
            if (drdrn) {
                double drndr;
                double rn2=r2rn(r1, values, sc, &drndr);
                printf("rn:%f %f\n",rn,rn2);
                *drdrn = 1 / drndr;
            }
            return r1;
        }
    }//internal ns
    EXP PTIso interpPTIso(double x, const PTIso& PT0, const PTIso& PT1){
        double xp = 1 - x;
		coord::UVSph cs1(x * PT0.cs.Delta + xp * PT1.cs.Delta);
		double R1 = x * PT0.sc.x0 + xp * PT1.sc.x0;
		std::vector<double> p, pr;
		if (PT0.N == PT1.N)
			for (int i = 0; i < PT0.N; i++)
				p.push_back(x * PT0.paramsF[i] + xp * PT1.paramsF[i]);
		else if (PT0.N > PT1.N) {
			for (int i = 0; i < PT1.N; i++)
				p.push_back(x * PT0.paramsF[i] + xp * PT1.paramsF[i]);
			for (int i = PT1.N; i < PT0.N; i++)
				p.push_back(x * PT0.paramsF[i]);
		}
		else {
			for (int i = 0; i < PT0.N; i++)
				p.push_back(x * PT0.paramsF[i] + xp * PT1.paramsF[i]);
			for (int i = PT0.N; i < PT1.N; i++)
				p.push_back(xp * PT1.paramsF[i]);
		}
		if (PT0.Nr == PT1.Nr)
			for (int i = 0; i < PT0.Nr; i++)
				pr.push_back(x * PT0.paramsFr[i] + xp * PT1.paramsFr[i]);
		else if (PT0.Nr > PT1.Nr) {
			for (int i = 0; i < PT1.Nr; i++)
				pr.push_back(x * PT0.paramsFr[i] + xp * PT1.paramsFr[i]);
			for (int i = PT1.Nr; i < PT0.Nr; i++)
				pr.push_back(x * PT0.paramsFr[i]);
		}
		else {
			for (int i = 0; i < PT0.Nr; i++)
				pr.push_back(x * PT0.paramsFr[i] + xp * PT1.paramsFr[i]);
			for (int i = PT0.Nr; i < PT1.Nr; i++)
				pr.push_back(xp * PT1.paramsFr[i]);
		}
        math::ScalingInfTh sc(R1);
		return PTIso(cs1, sc, p, pr);
    }
    EXP PTHarm interpPTHarm(double x, const PTHarm& PT0, const PTHarm& PT1){
        double xp = 1 - x;
		coord::UVSph cs1(x * PT0.cs.Delta + xp * PT1.cs.Delta);
		double R1 = x * PT0.sc.x0 + xp * PT1.sc.x0;
		double z1 = x * PT0.sc.x0 + xp * PT1.sc.x0;
		std::vector<double> p, pr;
		if (PT0.N == PT1.N)
			for (int i = 0; i < PT0.N; i++)
				p.push_back(x * PT0.paramsF[i] + xp * PT1.paramsF[i]);
		else if (PT0.N > PT1.N) {
			for (int i = 0; i < PT1.N; i++)
				p.push_back(x * PT0.paramsF[i] + xp * PT1.paramsF[i]);
			for (int i = PT1.N; i < PT0.N; i++)
				p.push_back(x * PT0.paramsF[i]);
		}
		else {
			for (int i = 0; i < PT0.N; i++)
				p.push_back(x * PT0.paramsF[i] + xp * PT1.paramsF[i]);
			for (int i = PT0.N; i < PT1.N; i++)
				p.push_back(xp * PT1.paramsF[i]);
		}
		if (PT0.Nr == PT1.Nr)
			for (int i = 0; i < PT0.Nr; i++)
				pr.push_back(x * PT0.paramsFr[i] + xp * PT1.paramsFr[i]);
		else if (PT0.Nr > PT1.Nr) {
			for (int i = 0; i < PT1.Nr; i++)
				pr.push_back(x * PT0.paramsFr[i] + xp * PT1.paramsFr[i]);
			for (int i = PT1.Nr; i < PT0.Nr; i++)
				pr.push_back(x * PT0.paramsFr[i]);
		}
		else {
			for (int i = 0; i < PT0.Nr; i++)
				pr.push_back(x * PT0.paramsFr[i] + xp * PT1.paramsFr[i]);
			for (int i = PT0.Nr; i < PT1.Nr; i++)
				pr.push_back(xp * PT1.paramsFr[i]);
		}
        math::ScalingInfTh sc(R1),scz(z1);
		return PTHarm(cs1, sc, scz, p, pr);
    }
    EXP PtrPointTransform interpPointTrans(double x, const PtrPointTransform PT0, const PtrPointTransform PT1) {
        if(PT0->name()!=PT1->name()){
            printf("Error:Cannot interpolate between a harmonic oscillator and isochrone point map.\n");
            exit(0);
        }
        static const char* IsText= "Isochrone";
        static const char* HOText = "Harmonic Oscillator";
        if(PT0->name()==IsText){
            double par0[2];
            int Nr0=PT0->FourierSizer(), Nz0=PT0->FourierSizez();
            std::vector<double> p0(Nz0),pr0(Nr0);
            PT0->getParams(par0,&pr0[0],&p0[0]);
            math::ScalingInfTh sc0(par0[1]);
            PTIso PTIs0(par0[0],sc0,p0,pr0);
            double par1[3];
            int Nr1=PT1->FourierSizer(), Nz1=PT1->FourierSizez();
            std::vector<double> p1(Nz1),pr1(Nr1);
            PT0->getParams(&par1[0],&pr1[0],&p1[0]);
            math::ScalingInfTh sc1(par1[1]);
            PTIso PTIs1(par1[0],sc1,p1,pr1);
            PTIso PTIsInterp=interpPTIso(x,PTIs0, PTIs1);
            return PtrPointTransform(
                new PTIso(PTIsInterp.cs,PTIsInterp.sc,PTIsInterp.paramsF,PTIsInterp.paramsFr));
        }
        else if(PT0->name()==HOText){
            double par0[3];
            int Nr0=PT0->FourierSizer(), Nz0=PT0->FourierSizez();
            std::vector<double> p0(Nz0),pr0(Nr0);
            PT0->getParams(par0,&pr0[0],&p0[0]);
            math::ScalingInfTh sc0(par0[1]), scz0(par0[2]);
            PTHarm PTH0(par0[0],sc0,scz0,p0,pr0);
            double par1[3];
            int Nr1=PT1->FourierSizer(), Nz1=PT1->FourierSizez();
            std::vector<double> p1(Nz1),pr1(Nr1);
            PT0->getParams(par1,&pr1[0],&p1[0]);
            math::ScalingInfTh sc1(par1[1]), scz1(par1[2]);
            PTHarm PTH1(par1[0],sc1,scz1,p1,pr1);
            PTHarm PTHarmInterp=interpPTHarm(x,PTH0, PTH1);
            return PtrPointTransform(
                new PTHarm(PTHarmInterp.cs,PTHarmInterp.sc,PTHarmInterp.scz,PTHarmInterp.paramsF,PTHarmInterp.paramsFr));
        }
        return PtrPointTransform();
	}
	double PTIso::t2v(const double theta, double* dvdt, double* d2vdt2, double *dvdP, double *ddvdtdP) const {
		double v = theta;
        double deriv=0,deriv2=0;
        v+=FourSer(theta,paramsF,dvdt?&deriv:NULL,d2vdt2?&deriv2:NULL,dvdP?dvdP:NULL,ddvdtdP?ddvdtdP:NULL);
        deriv+=1;
		if (dvdt)*dvdt = deriv;
		if (d2vdt2)*d2vdt2 = deriv2;
		return v;
	}
    double PTIso::v2t(const double v, double* dtdv) const {
		if (N==0) {
			if (dtdv) *dtdv = 1; return v;
		}
		vthet V(paramsF, v);
		double t = math::findRoot(V, 0., M_PI, 1e-9);
		if (dtdv) {
			double dvdt;
			t2v(t, &dvdt);
			*dtdv = 1. / dvdt;
		}
		return t;
	}
    double PTHarm::zntov(const double zn, double* dvdt, double* d2vdt2, double *dvdz0, 
        double *ddvdzdz0, double *dvdFC, double *ddvdtdFC) const {
        double dszdz, d2szdz2;
		double sz = M_PI * math::scale(scz, zn, &dszdz, &d2szdz2);
		dszdz *= M_PI;
		d2szdz2 *= M_PI;
        double dvdsz=0,d2vdsz2=0;
        bool needDeriv=(dvdt||d2vdt2||dvdz0||ddvdzdz0);
        double v=M_PI-sz+FourSer(sz,paramsF,needDeriv?&dvdsz:NULL,d2vdt2||ddvdzdz0?&d2vdsz2:NULL,
        (dvdFC)?dvdFC:NULL,(ddvdtdFC)?ddvdtdFC:NULL);
        if(ddvdtdFC){
            for(int i=0;i<N;i++) ddvdtdFC[i]*=dszdz;
        }
        dvdsz-=1;
        double dvdz=dvdsz*dszdz;
        double d2vdzn2=d2vdsz2 * pow_2(dszdz) + dvdsz * d2szdz2;
        if (dvdt)*dvdt = dvdz;
		if (d2vdt2)*d2vdt2 = d2vdzn2;
        if(dvdz0||ddvdzdz0){ 
            double dszdz0=dszdz*(-zn/scz.x0);
            double ddszdzdz0=d2szdz2*(-zn/scz.x0)-dszdz/scz.x0;
            if(dvdz0)*dvdz0=dvdsz*dszdz0;
            if(ddvdzdz0) *ddvdzdz0=d2vdsz2*dszdz0*dszdz+dvdsz*ddszdzdz0;
        }
		return v;
    }
	double PTHarm::vtozn(const double v, double* dzdv) const {
		zfind zr(paramsF, v);
		double s = math::findRoot(zr, 0, M_PI, 1e-9);
		double zn = math::unscale(scz, s/M_PI);
		if (dzdv) {
			double drdR;
			zntov(zn, &drdR);
			*dzdv = 1. / drdR;
		}
		return zn;
	}
	coord::PosMomSph PTIso::revmapSph(const coord::PosMomCyl &Rz) const {
		double R2 = pow_2(Rz.R), z2 = pow_2(Rz.z);
		double B = R2 + z2 - cs.Delta2;
		double r2 = .5 * (B + sqrt(B * B + 4 * R2 * cs.Delta2));
		double r0 = sqrt(r2);
		double drdrn;
		double r = rn2r(r0, paramsFr, sc, &drdrn);
		double rt = sqrt(r2 + cs.Delta2);
		double v = acos(Rz.z / rt);double dtdv;
		double snt, cst; math::sincos(v, snt, cst);
		double pr = snt * Rz.pR + r0 / rt * cst * Rz.pz;
		double ptheta = r0 * cst * Rz.pR - rt * snt * Rz.pz;
		double psi = v2t(v, &dtdv);
		return coord::PosMomSph(r, psi, Rz.phi, pr / drdrn, ptheta / dtdv, Rz.pphi);
	}
	coord::PosMomCyl PTHarm::revmap(const coord::PosMomCyl &Rz) const {
		double R2 = pow_2(Rz.R), z2 = pow_2(Rz.z);
		double B = R2 + z2 - pow_2(cs.Delta);
		double r2 = .5 * (B + sqrt(B * B + 4 * R2 * pow_2(cs.Delta)));
		double r0 = sqrt(r2);
		double drdrn;
		double r = rn2r(r0, paramsFr, sc, &drdrn);
		double rt = sqrt(r2 + pow_2(cs.Delta));
		double v = acos(Rz.z / rt), dzdv;
		double snt, cst; math::sincos(v, snt, cst);
		double pr = snt * Rz.pR + r0 / rt * cst * Rz.pz;
		double ptheta = r0 * cst * Rz.pR - rt * snt * Rz.pz;
		double z = vtozn(v, &dzdv);
		return coord::PosMomCyl(r, z, Rz.phi, pr / drdrn, ptheta / dzdv, Rz.pphi);
	}
	coord::PosMomCyl PTIso::map(const coord::PosMomSph &rp, coord::PosMomCyl* dRzdP,
            coord::PosMomCyl *dRzdFr,coord::PosMomCyl *dRzdFz) const {
		double drndr=0, drndx0=0, ddrndrdx0=0;
        std::vector<double> drndP(Nr),ddrndrdP(Nr),dvdP(N),ddvdtdP(N);
		double r = r2rn(rp.r, paramsFr, sc, &drndr, NULL, dRzdP?&drndx0:NULL, dRzdP?&ddrndrdx0:NULL,
            dRzdFr?&drndP[0]:NULL,dRzdFr?&ddrndrdP[0]:NULL);
		double sq = pow_2(r) + cs.Delta2, rt = sqrt(sq);
		double dvdt;
		double v = t2v(rp.theta, &dvdt, NULL, &dvdP[0], &ddvdtdP[0]);
		double snt, cst; math::sincos(v, snt, cst);
		double R = r * snt, z = rt * cst;
		double bot1 = pow_2(r) + cs.Delta2 * pow_2(snt); //(sq * pow_2(snt) + pow_2(rp.r * cst));
        double bot2 = (pow_2(r * cst) / rt + rt * pow_2(snt));
		double pR = (sq * snt * rp.pr / drndr + r * cst * rp.ptheta / dvdt) / bot1;
		double pz = (r * cst * rp.pr / drndr - snt * rp.ptheta / dvdt) / bot2;
        double dbot1dv = 2 * cs.Delta2 * cst * snt;
        double dbot2dv = dbot1dv / rt;
        double dbot1dr = 2 * r;
		double dbot2dr = -r / pow_3(rt) * bot1 + 2 * r / rt;
		if(dRzdP){
			if(N==0&&Nr==0){
				dRzdP[0].R = dRzdP[0].phi = dRzdP[0].pphi = 0; 
				dRzdP[0].z = cs.Delta / rt * cst;
				dRzdP[0].pR = 2 * cs.Delta * snt * rp.pr/ (drndr*bot1)
					- pR / bot1 * 2 * cs.Delta * pow_2(snt);
				dRzdP[0].pz = -pz / bot2 * cs.Delta / rt * (-pow_2(r * cst) / sq
					+ pow_2(snt));
				dRzdP[0].phi=dRzdP[0].pphi;
			}else{
				dRzdP[0].R = dRzdP[0].phi = dRzdP[0].pphi = 0; 
				dRzdP[0].z = cs.Delta / rt * cst;
				dRzdP[0].pR = 2 * cs.Delta * snt * rp.pr/ (drndr*bot1)
					- pR / bot1 * 2 * cs.Delta * pow_2(snt);
				dRzdP[0].pz = -pz / bot2 * cs.Delta / rt * (-pow_2(r * cst) / sq
					+ pow_2(snt));
				dRzdP[0].phi=dRzdP[0].pphi=dRzdP[1].phi=dRzdP[1].pphi=0;
				dRzdP[1].R=snt*drndx0;
				dRzdP[1].z=r/rt*cst*drndx0;
				dRzdP[1].phi=dRzdP[1].pphi=0;
				dRzdP[1].pR=((2*r*snt*rp.pr/drndr + cst * rp.ptheta/dvdt) / bot1 - pR/bot1*dbot1dr)*drndx0
				-sq * snt * rp.pr/(bot1*pow_2(drndr))*ddrndrdx0;
				dRzdP[1].pz=(cst * rp.pr /(drndr*bot2)-pz/bot2*dbot2dr)*drndx0-r * cst * rp.pr/(bot2*pow_2(drndr))*ddrndrdx0;
			}
        }
        if(dRzdFr)
            for(int i=0;i<Nr;i++){
                dRzdFr[i].R=snt*drndP[i];
                dRzdFr[i].z=r/rt*cst*drndP[i];
                dRzdFr[i].phi=dRzdFr[i].pphi=0;
                dRzdFr[i].pR=((2*r*snt*rp.pr/drndr+cst * rp.ptheta/dvdt) / bot1-pR/bot1*dbot1dr)*drndP[i]
                    -sq * snt * rp.pr/(bot1*pow_2(drndr))*ddrndrdP[i];
                dRzdFr[i].pz=(cst * rp.pr / (drndr*bot2)-pz/bot2*dbot2dr)*drndP[i]
                    -r * cst * rp.pr/(bot2*pow_2(drndr))*ddrndrdP[i];
				dRzdFr[i].phi=dRzdFr[i].pphi=0;
            }
        if(dRzdFz)
            for(int i=0;i<N;i++){
                dRzdFz[i].R=r*cst*dvdP[i];
                dRzdFz[i].z=-rt*snt*dvdP[i];
                dRzdFz[i].phi=dRzdFz[i].pphi=0;
                dRzdFz[i].pR=((sq * cst * rp.pr/drndr - r * snt * rp.ptheta/dvdt) / bot1-pR/bot1*dbot1dv)*dvdP[i]
                    -r*cst*rp.ptheta/(bot1*pow_2(dvdt))*ddvdtdP[i];
                dRzdFz[i].pz=((-r*snt*rp.pr/drndr-cst*rp.ptheta/dvdt) / bot2-pz/bot2*dbot2dv)*dvdP[i]
                    +snt*rp.ptheta/(bot2*pow_2(dvdt))*ddvdtdP[i];
				dRzdFz[i].phi=dRzdFz[i].pphi=0;
            }
		return coord::PosMomCyl(R, z, rp.phi, pR, pz, rp.pphi);
	}

	coord::PosMomCyl PTIso::map(const coord::PosMomSph &rp, math::Matrix<double>& dRzdrt,
		coord::PosMomCyl* dRzdP, coord::PosMomCyl* dRzdFr, coord::PosMomCyl* dRzdFz) const { 
		double drndr, d2rndr2, drndx0=0, ddrndrdx0=0;;
        std::vector<double> drndP(Nr),ddrndrdP(Nr),dvdP(N),ddvdtdP(N);
		double r = r2rn(rp.r, paramsFr, sc, &drndr, &d2rndr2, dRzdP?&drndx0:NULL, dRzdP?&ddrndrdx0:NULL,
            dRzdFr?&drndP[0]:NULL,dRzdFr?&ddrndrdP[0]:NULL);
		double sq = pow_2(r) + cs.Delta2, rt = sqrt(sq);
		double dvdt, d2vdt2;
		double v = t2v(rp.theta, &dvdt, &d2vdt2,dRzdFz?&dvdP[0]:NULL,dRzdFz?&ddvdtdP[0]:NULL);
		double snt, cst; math::sincos(v, snt, cst);
		double R = r * snt, z = rt * cst;
		double bot1 = pow_2(r) + cs.Delta2 * pow_2(snt); //(sq * pow_2(snt) + pow_2(rp.r * cst));
		double dbot1dr = 2 * r;
		double dbot1dtheta = 2 * cs.Delta2 * cst * snt;
		double bot2 = (pow_2(r * cst) / rt + rt * pow_2(snt));
		double dbot2dr = -r / pow_3(rt) * bot1 + 2 * r / rt;
		double dbot2dtheta = dbot1dtheta / rt;
		double pR = (sq * snt * rp.pr / drndr + r * cst * rp.ptheta / dvdt) / bot1;
		double pz = (r * cst * rp.pr / drndr - snt * rp.ptheta / dvdt) / bot2;
		dRzdrt(0, 0) = drndr * snt; dRzdrt(0, 1) = r * cst * dvdt; dRzdrt(0, 2) = 0;//dR/dr dR/dpsi
		dRzdrt(1, 0) = drndr * r / rt * cst;   dRzdrt(1, 1) = -rt * snt * dvdt; dRzdrt(1, 2) = 0;//dz/dr dz/dpsi
		dRzdrt(2, 0) = ((2 * r * snt * rp.pr / drndr + cst * rp.ptheta / dvdt) / bot1
			- pR / bot1 * dbot1dr) * drndr - sq * snt * rp.pr / (pow_2(drndr) * bot1) * d2rndr2; //dpR/dr 
		dRzdrt(2, 1) = ((sq * cst * rp.pr / drndr - r * snt * rp.ptheta / dvdt) / bot1
			- pR / bot1 * dbot1dtheta) * dvdt - r * cst / bot1 * rp.ptheta / pow_2(dvdt) * d2vdt2;//dpR/dpsi
		dRzdrt(3, 0) = (cst * rp.pr / (bot2 * drndr) - pz / bot2 * dbot2dr) * drndr 
        - r * cst * rp.pr / (bot2 * pow_2(drndr)) * d2rndr2;// dpz/dr
		dRzdrt(3, 1) = ((-r * snt * rp.pr / drndr - cst * rp.ptheta / dvdt) / bot2
			- pz / bot2 * dbot2dtheta) * dvdt + snt / bot2 * rp.ptheta / pow_2(dvdt) * d2vdt2;// dpz/dpsi
		dRzdrt(0, 2) = 0; dRzdrt(0, 3) = 0;//dR/dpr ppsi
		dRzdrt(1, 2) = 0; dRzdrt(1, 3) = 0;//dz/dpr ppsi
		dRzdrt(2, 2) = sq * snt / (bot1 * drndr); dRzdrt(2, 3) = r * cst / (bot1 * dvdt);//dpR/dpr ppsi
		dRzdrt(3, 2) = r * cst / (bot2 * drndr); dRzdrt(3, 3) = -snt / (bot2 * dvdt);//dpz/dpr ppsi
        
		if(dRzdP){
			if(N==0&&Nr==0){
				dRzdP[0].R = dRzdP[0].phi = dRzdP[0].pphi = 0; 
				dRzdP[0].z = cs.Delta / rt * cst;
				dRzdP[0].pR = 2 * cs.Delta * snt * rp.pr/ (drndr*bot1)
					- pR / bot1 * 2 * cs.Delta * pow_2(snt);
				dRzdP[0].pz = -pz / bot2 * cs.Delta / rt * (-pow_2(r * cst) / sq
					+ pow_2(snt));
				dRzdP[0].phi=dRzdP[0].pphi;
			}else{
				dRzdP[0].R = dRzdP[0].phi = dRzdP[0].pphi = 0; 
				dRzdP[0].z = cs.Delta / rt * cst;
				dRzdP[0].pR = 2 * cs.Delta * snt * rp.pr/ (drndr*bot1)
					- pR / bot1 * 2 * cs.Delta * pow_2(snt);
				dRzdP[0].pz = -pz / bot2 * cs.Delta / rt * (-pow_2(r * cst) / sq
					+ pow_2(snt));
				dRzdP[0].phi=dRzdP[0].pphi=dRzdP[1].phi=dRzdP[1].pphi=0;
				dRzdP[1].R=snt*drndx0;
				dRzdP[1].z=r/rt*cst*drndx0;
				dRzdP[1].phi=dRzdP[1].pphi=0;
				dRzdP[1].pR=((2*r*snt*rp.pr/drndr + cst * rp.ptheta/dvdt) / bot1 - pR/bot1*dbot1dr)*drndx0
				-sq * snt * rp.pr/(bot1*pow_2(drndr))*ddrndrdx0;
				dRzdP[1].pz=(cst * rp.pr /(drndr*bot2)-pz/bot2*dbot2dr)*drndx0-r * cst * rp.pr/(bot2*pow_2(drndr))*ddrndrdx0;
			}
        }
        if(dRzdFr)
            for(int i=0;i<Nr;i++){
                dRzdFr[i].R=snt*drndP[i];
                dRzdFr[i].z=r/rt*cst*drndP[i];
                dRzdFr[i].phi=dRzdFr[i].pphi=0;
                dRzdFr[i].pR=((2*r*snt*rp.pr/drndr+cst * rp.ptheta/dvdt) / bot1-pR/bot1*dbot1dr)*drndP[i]
                    -sq * snt * rp.pr/(bot1*pow_2(drndr))*ddrndrdP[i];
                dRzdFr[i].pz=(cst * rp.pr / (drndr*bot2)-pz/bot2*dbot2dr)*drndP[i]
                    -r * cst * rp.pr/(bot2*pow_2(drndr))*ddrndrdP[i];
				dRzdFr[i].phi=dRzdFr[i].pphi=0;
            }
        if(dRzdFz)
            for(int i=0;i<N;i++){
                dRzdFz[i].R=r*cst*dvdP[i];
                dRzdFz[i].z=-rt*snt*dvdP[i];
                dRzdFz[i].phi=dRzdFz[i].pphi=0;
                dRzdFz[i].pR=((sq * cst * rp.pr/drndr - r * snt * rp.ptheta/dvdt) / bot1-pR/bot1*dbot1dtheta)*dvdP[i]
                    -r*cst*rp.ptheta/(bot1*pow_2(dvdt))*ddvdtdP[i];
                dRzdFz[i].pz=((-r*snt*rp.pr/drndr-cst*rp.ptheta/dvdt) / bot2-pz/bot2*dbot2dtheta)*dvdP[i]
                    +snt*rp.ptheta/(bot2*pow_2(dvdt))*ddvdtdP[i];
				dRzdFz[i].phi=dRzdFz[i].pphi=0;
            }
		return coord::PosMomCyl(R, z, rp.phi, pR, pz, rp.pphi);
	}
	coord::PosMomCyl PTHarm::map(const coord::PosMomCyl &point, coord::PosMomCyl* dRzdP,
        coord::PosMomCyl *dRzdFr,coord::PosMomCyl *dRzdFz) const {
		double dvdz=0, dvdz0=0, ddvdzdz0=0;
        std::vector<double> dvdP(N),d2vdzdP(N),drndP(Nr),ddrndrdP(Nr);
		double v = zntov(point.z, &dvdz, NULL, dRzdP?&dvdz0:NULL, dRzdP?&ddvdzdz0:NULL, 
            dRzdFz?&dvdP[0]:NULL,dRzdFz?&d2vdzdP[0]:NULL);
		double snt, cst; math::sincos(v, snt, cst);
		double drdx, drndx0, ddrndrdx0;
		double r = r2rn(point.R, paramsFr, sc, &drdx, NULL, dRzdP?&drndx0:NULL, 
            (dRzdP)?&ddrndrdx0:NULL, (dRzdFr)?&drndP[0]:NULL, (dRzdFr)?&ddrndrdP[0]:NULL);
		double pr = point.pR / drdx;
		double sq = pow_2(r) + cs.Delta2, rt = sqrt(sq);
		double R = r * snt, z = rt * cst;
		double bot1 = (sq * pow_2(snt) + pow_2(r * cst));
		double bot2 = (pow_2(r * cst) / rt + rt * pow_2(snt));
		double pv = point.pz / dvdz;
		double pR = (sq * snt * pr + r * cst * pv) / bot1;
		double pz = (r * cst * pr - snt * pv) / bot2;
        double dbot1dv = 2 * cs.Delta2 * cst * snt;
        double dbot2dv = dbot1dv / rt;
        double dbot1dr = 2 * r;
		double dbot2dr = -r / pow_3(rt) * bot1 + 2 * r / rt;
        if(dRzdP){
			dRzdP[0].R =dRzdP[0].phi=dRzdP[0].pphi=0; 
            dRzdP[0].z = cs.Delta / rt * cst;
			dRzdP[0].pR = 2 * cs.Delta * snt * pr/ bot1
				- pR / bot1 * 2 * cs.Delta * pow_2(snt);
			dRzdP[0].pz = -pz / bot2 * cs.Delta / rt * (-pow_2(r * cst) / sq
				+ pow_2(snt));
			dRzdP[0].phi=dRzdP[0].pphi=0;
            dRzdP[1].R=snt*drndx0;
            dRzdP[1].z=r/rt*cst*drndx0;
            dRzdP[1].phi=dRzdP[1].pphi=0;
            dRzdP[1].pR=((2*r*snt*pr+cst * pv) / bot1-pR/bot1*dbot1dr)*drndx0-sq * snt * pr/(bot1*drdx)*ddrndrdx0;
            dRzdP[1].pz=(cst * pr / bot2-pz/bot2*dbot2dr)*drndx0-r * cst * pr/(bot2*drdx)*ddrndrdx0;
            dRzdP[2].R=r*cst*dvdz0;
            dRzdP[2].z=-rt*snt*dvdz0;
            dRzdP[2].phi= dRzdP[2].pphi=0;
            dRzdP[2].pR=((sq * cst * pr - r * snt * pv) / bot1-pR/bot1*dbot1dv)*dvdz0-r*cst*pv/(bot1*dvdz)*ddvdzdz0;
            dRzdP[2].pz=((-r*snt*pr-cst*pv) / bot2-pz/bot2*dbot2dv)*dvdz0+snt*pv/(bot2*dvdz)*ddvdzdz0;
        }
        if(dRzdFz)
            for(int i=0;i<N;i++){
                dRzdFz[i].R=r*cst*dvdP[i];
                dRzdFz[i].z=-rt*snt*dvdP[i];
                dRzdFz[i].phi=dRzdFz[i].pphi=0;
                dRzdFz[i].pR=((sq * cst * pr - r * snt * pv) / bot1-pR/bot1*dbot1dv)*dvdP[i]-r*cst*pv/(bot1*dvdz)*d2vdzdP[i];
                dRzdFz[i].pz=((-r*snt*pr-cst*pv) / bot2-pz/bot2*dbot2dv)*dvdP[i]+snt*pv/(bot2*dvdz)*d2vdzdP[i];
            }
        if(dRzdFr)
            for(int i=0;i<Nr;i++){
                dRzdFr[i].R=snt*drndP[i];
                dRzdFr[i].z=r/rt*cst*drndP[i];
                dRzdFr[i].phi=dRzdFr[i].pphi=0;
                dRzdFr[i].pR=((2*r*snt*pr+cst * pv) / bot1-pR/bot1*dbot1dr)*drndP[i]-sq * snt * pr/(bot1*drdx)*ddrndrdP[i];
                dRzdFr[i].pz=(cst * pr / bot2-pz/bot2*dbot2dr)*drndP[i]-r * cst * pr/(bot2*drdx)*ddrndrdP[i];
            }
		return coord::PosMomCyl(R, z, point.phi, pR, pz, point.pphi);
	}
	coord::PosMomCyl PTHarm::map(const coord::PosMomCyl &rp, math::Matrix<double>& dRzdrt, coord::PosMomCyl* dRzdP,
        coord::PosMomCyl *dRzdFr,coord::PosMomCyl *dRzdFz) const {
		double dvdz, d2vdzn2, dvdz0=0, ddvdzdz0=0;
        std::vector<double> dvdP(N),d2vdzdP(N),drndP(Nr),ddrndrdP(Nr);
		double v = zntov(rp.z, &dvdz, &d2vdzn2, dRzdP?&dvdz0:NULL, dRzdP?&ddvdzdz0:NULL, 
            dRzdFz?&dvdP[0]:NULL,dRzdFz?&d2vdzdP[0]:NULL);
		double snt, cst; math::sincos(v, snt, cst);
		double drdx, d2rdx2, drndx0, ddrndrdx0;
		double r = r2rn(rp.R, paramsFr, sc, &drdx, &d2rdx2, dRzdP?&drndx0:NULL, 
            dRzdP?&ddrndrdx0:NULL, dRzdFr?&drndP[0]:NULL, dRzdFr?&ddrndrdP[0]:NULL);
            printf("r:%f %f %f %d\n",r,rp.R,math::scale(sc,rp.R),paramsFr.size());
        for(int i=0;i<paramsFr.size();i++)printf("pnew:%f\n",paramsFr[i]);
		double pr = rp.pR / drdx;
		double sq = pow_2(r) + cs.Delta2, rt = sqrt(sq);
		double R = r * snt, z = rt * cst;
		double bot1 = (sq * pow_2(snt) + pow_2(r * cst));
		double bot2 = (pow_2(r * cst) / rt + rt * pow_2(snt));
		double pv = rp.pz / dvdz;
		double pR = (sq * snt * pr + r * cst * pv) / bot1;
		double pz = (r * cst * pr - snt * pv) / bot2;
		double dbot1dtheta = 2 * cs.Delta2 * cst * snt;
		double dbot1dr = 2 * r;
		double dbot2dtheta = dbot1dtheta / rt;
		double dbot2dr = -r / pow_3(rt) * bot1 + 2 * r / rt;
		dRzdrt(0, 0) = drdx * snt; dRzdrt(0, 1) = r * cst * dvdz; dRzdrt(0, 2) = 0;//dR/dr dR/dpsi
		dRzdrt(1, 0) = drdx * r / rt * cst;   dRzdrt(1, 1) = -rt * snt * dvdz; dRzdrt(1, 2) = 0;//dz/dR1 dz/dz1
		dRzdrt(2, 0) = ((2 * r * snt * rp.pR / drdx + cst * rp.pz / dvdz) / bot1
			- pR / bot1 * dbot1dr) * drdx - sq * snt * rp.pR / (pow_2(drdx) * bot1) * d2rdx2; //dpR/dR1 
		dRzdrt(2, 1) = ((sq * cst * rp.pR / drdx - r * snt * rp.pz / dvdz) / bot1
			- pR / bot1 * dbot1dtheta) * dvdz - r * cst / bot1 * rp.pz / pow_2(dvdz) * d2vdzn2;//dpR/dz1
		dRzdrt(3, 0) = (cst * rp.pR / (bot2 * drdx) - pz / bot2 * dbot2dr) * drdx - r * cst * rp.pR / (bot2 * pow_2(drdx)) * d2rdx2;// dpz/dr
		dRzdrt(3, 1) = ((-r * snt * rp.pR / drdx - cst * rp.pz / dvdz) / bot2
			- pz / bot2 * dbot2dtheta) * dvdz + snt / bot2 * rp.pz / pow_2(dvdz) * d2vdzn2;// dpz/dpsi
		dRzdrt(0, 2) = 0; dRzdrt(0, 3) = 0;//dR/dpr ppsi
		dRzdrt(1, 2) = 0; dRzdrt(1, 3) = 0;//dz/dpr ppsi
		dRzdrt(2, 2) = sq * snt / (bot1 * drdx); dRzdrt(2, 3) = r * cst / (bot1 * dvdz);//dpR/dpr ppsi
		dRzdrt(3, 2) = r * cst / (bot2 * drdx); dRzdrt(3, 3) = -snt / (bot2 * dvdz);//dpz/dpr ppsi
        if(dRzdP){
			dRzdP[0].R =dRzdP[0].phi=dRzdP[0].pphi=0; 
            dRzdP[0].z = cs.Delta / rt * cst;
			dRzdP[0].pR = 2 * cs.Delta * snt * pr/ bot1
				- pR / bot1 * 2 * cs.Delta * pow_2(snt);
			dRzdP[0].pz = -pz / bot2 * cs.Delta / rt * (-pow_2(r * cst) / sq
				+ pow_2(snt));
            dRzdP[1].R=snt*drndx0;
            dRzdP[1].z=r/rt*cst*drndx0;
            dRzdP[1].phi=dRzdP[1].pphi=0;
            dRzdP[1].pR=((2*r*snt*pr+cst * pv) / bot1-pR/bot1*dbot1dr)*drndx0-sq * snt * pr/(bot1*drdx)*ddrndrdx0;
            dRzdP[1].pz=(cst * pr / bot2-pz/bot2*dbot2dr)*drndx0-r * cst * pr/(bot2*drdx)*ddrndrdx0;
            dRzdP[2].R=r*cst*dvdz0;
            dRzdP[2].z=-rt*snt*dvdz0;
            dRzdP[2].phi= dRzdP[2].pphi=0;
            dRzdP[2].pR=((sq * cst * pr - r * snt * pv) / bot1-pR/bot1*dbot1dtheta)*dvdz0-r*cst*pv/(bot1*dvdz)*ddvdzdz0;
            dRzdP[2].pz=((-r*snt*pr-cst*pv) / bot2-pz/bot2*dbot2dtheta)*dvdz0+snt*pv/(bot2*dvdz)*ddvdzdz0;
            
        }
        if(dRzdFz)
            for(int i=0;i<N;i++){
                dRzdFz[i].R=r*cst*dvdP[i];
                dRzdFz[i].z=-rt*snt*dvdP[i];
                dRzdFz[i].phi=dRzdFz[i].pphi=0;
                dRzdFz[i].pR=((sq * cst * pr - r * snt * pv) / bot1-pR/bot1*dbot1dtheta)*dvdP[i]-r*cst*pv/(bot1*dvdz)*d2vdzdP[i];
                dRzdFz[i].pz=((-r*snt*pr-cst*pv) / bot2-pz/bot2*dbot2dtheta)*dvdP[i]+snt*pv/(bot2*dvdz)*d2vdzdP[i];
            }
        if(dRzdFr)
            for(int i=0;i<Nr;i++){
                dRzdFr[i].R=snt*drndP[i];
                dRzdFr[i].z=r/rt*cst*drndP[i];
                dRzdFr[i].phi=dRzdFr[i].pphi=0;
                dRzdFr[i].pR=((2*r*snt*pr+cst * pv) / bot1-pR/bot1*dbot1dr)*drndP[i]-sq * snt * pr/(bot1*drdx)*ddrndrdP[i];
                dRzdFr[i].pz=(cst * pr / bot2-pz/bot2*dbot2dr)*drndP[i]-r * cst * pr/(bot2*drdx)*ddrndrdP[i];
            }
		return coord::PosMomCyl(R, z, rp.phi, pR, pz, rp.pphi);
	}
	EXP PtrToyMap interpPtrToyMap(double x, const PtrToyMap& PtrTM0, const PtrToyMap& PtrTM1) {
		if (PtrTM0->getToyMapType() == PtrTM1->getToyMapType()) {
			double val1,val2;
			double valn1,valn2;
			PtrTM0->getParams(&val1,&val2);
			PtrTM1->getParams(&valn1,&valn2);
			if (PtrTM0->getToyMapType()==ToyPotType::Is){
				Isochrone Is0(val1,val2);
				Isochrone Is1(valn1,valn2);
                PtrPointTransform PtrPT0=PtrTM0->getPointTrans();
                PtrPointTransform PtrPT1=PtrTM1->getPointTrans();
                double par0[2];
                int Nr0=PtrPT0->FourierSizer(), Nz0=PtrPT0->FourierSizez();
                std::vector<double> p0(Nz0),pr0(Nr0);
                PtrPT0->getParams(par0,&pr0[0],&p0[0]);
                math::ScalingInfTh sc0(par0[1]);
                PTIso PTIs0(par0[0],sc0,p0,pr0);
                double par1[2];
                int Nr1=PtrPT1->FourierSizer(), Nz1=PtrPT1->FourierSizez();
                std::vector<double> p1(Nz1),pr1(Nr1);
                PtrPT1->getParams(par1,&pr1[0],&p1[0]);
                math::ScalingInfTh sc1(par1[1]);
                PTIso PTIs1(par1[0],sc1,p1,pr1);
				return PtrToyMap(new ToyMapIso(interpIsochrone(x, Is0, Is1),interpPTIso(x,PTIs0,PTIs1)));
			}
            else if(PtrTM0->getToyMapType()==ToyPotType::HO){
				HarmonicOscilattor HOs0(val1,val2);
			    HarmonicOscilattor HOs1(valn1,valn2);
                PtrPointTransform PtrPT0=interpPointTrans(x,PtrTM0->getPointTrans(),PtrTM1->getPointTrans());
                PtrPointTransform PtrPT1=PtrTM1->getPointTrans();
                double par0[3];
                int Nr0=PtrPT0->FourierSizer(), Nz0=PtrPT0->FourierSizez();
                std::vector<double> p0(Nz0),pr0(Nr0);
                PtrPT0->getParams(par0,&pr0[0],&p0[0]);
                math::ScalingInfTh sc0(par0[1]), scz0(par0[2]);
                PTHarm PTH0(par0[0],sc0,scz0,p0,pr0);
                double par1[3];
                int Nr1=PtrPT1->FourierSizer(), Nz1=PtrPT1->FourierSizez();
                std::vector<double> p1(Nz1),pr1(Nr1);
                PtrPT1->getParams(par1,&pr1[0],&p1[0]);
                math::ScalingInfTh sc1(par1[1]), scz1(par1[2]);
                PTHarm PTH1(par1[0],sc1,scz1,p1,pr1);
				return PtrToyMap(new ToyMapHarm(interpHarmonicOscillator(x, HOs0, HOs1),interpPTHarm(x,PTH0,PTH1)));
			}
		}
		printf("Error:cannot interpolate between a Harmonic oscillator and Isochrone Toy Map.\n");
		return PtrToyMap();
	}

    coord::PosMomCyl ToyMapIso::from_aaT(const ActionAngles& aaT, coord::PosMomCyl* dRzdPpot,
			coord::PosMomCyl* dRzdPPT, coord::PosMomCyl* dRzdFr, coord::PosMomCyl* dRzdFz) const {
		coord::PosMomSph drdJs, drdb;
		coord::PosMomSph rp(Is.aa2pq(aaT, drdJs, drdb));
		math::Matrix<double> dRzdrt(4, 4);
		coord::PosMomCyl Rz(PT.map(rp, dRzdrt, dRzdPPT, dRzdFr, dRzdFz));
        if(dRzdPpot){
            dRzdPpot[0].R = dRzdrt(0, 0) * drdJs.r + dRzdrt(0, 1) * drdJs.theta
                + dRzdrt(0, 2) * drdJs.pr + dRzdrt(0, 3) * drdJs.ptheta;
            dRzdPpot[0].z = dRzdrt(1, 0) * drdJs.r + dRzdrt(1, 1) * drdJs.theta
                + dRzdrt(1, 2) * drdJs.pr + dRzdrt(1, 3) * drdJs.ptheta;
            dRzdPpot[0].pR = dRzdrt(2, 0) * drdJs.r + dRzdrt(2, 1) * drdJs.theta
                + dRzdrt(2, 2) * drdJs.pr + dRzdrt(2, 3) * drdJs.ptheta;
            dRzdPpot[0].pz = dRzdrt(3, 0) * drdJs.r + dRzdrt(3, 1) * drdJs.theta
                + dRzdrt(3, 2) * drdJs.pr + dRzdrt(3, 3) * drdJs.ptheta;
            dRzdPpot[1].R = dRzdrt(0, 0) * drdb.r + dRzdrt(0, 1) * drdb.theta
                + dRzdrt(0, 2) * drdb.pr + dRzdrt(0, 3) * drdb.ptheta;
            dRzdPpot[1].z = dRzdrt(1, 0) * drdb.r + dRzdrt(1, 1) * drdb.theta
                + dRzdrt(1, 2) * drdb.pr + dRzdrt(1, 3) * drdb.ptheta;
            dRzdPpot[1].pR = dRzdrt(2, 0) * drdb.r + dRzdrt(2, 1) * drdb.theta
                + dRzdrt(2, 2) * drdb.pr + dRzdrt(2, 3) * drdb.ptheta;
            dRzdPpot[1].pz = dRzdrt(3, 0) * drdb.r + dRzdrt(3, 1) * drdb.theta
                + dRzdrt(3, 2) * drdb.pr + dRzdrt(3, 3) * drdb.ptheta;
        }
		return Rz;
	}
	coord::PosMomCyl ToyMapHarm::from_aaT(const ActionAngles& aaT, coord::PosMomCyl* dRzdPpot,
			coord::PosMomCyl* dRzdPPT, coord::PosMomCyl* dRzdFr, coord::PosMomCyl* dRzdFz) const {
		coord::PosMomCyl drdomegar, drdomegaz;
		coord::PosMomCyl rp(HOs.aa2pq(aaT, drdomegar, drdomegaz));
		math::Matrix<double> dRzdrt(4, 4);
		coord::PosMomCyl Rz(PT.map(rp, dRzdrt, dRzdPPT, dRzdFr, dRzdFz));
        if(dRzdPpot){
		dRzdPpot[0].R = dRzdrt(0, 0) * drdomegar.R + dRzdrt(0, 1) * drdomegar.z
			+ dRzdrt(0, 2) * drdomegar.pR + dRzdrt(0, 3) * drdomegar.pz;
		dRzdPpot[0].z = dRzdrt(1, 0) * drdomegar.R + dRzdrt(1, 1) * drdomegar.z
			+ dRzdrt(1, 2) * drdomegar.pR + dRzdrt(1, 3) * drdomegar.pz;
		dRzdPpot[0].pR = dRzdrt(2, 0) * drdomegar.R + dRzdrt(2, 1) * drdomegar.z
			+ dRzdrt(2, 2) * drdomegar.pR + dRzdrt(2, 3) * drdomegar.pz;
		dRzdPpot[0].pz = dRzdrt(3, 0) * drdomegar.R + dRzdrt(3, 1) * drdomegar.z
			+ dRzdrt(3, 2) * drdomegar.pR + dRzdrt(3, 3) * drdomegar.pz;
		dRzdPpot[1].R = dRzdrt(0, 0) * drdomegaz.R + dRzdrt(0, 1) * drdomegaz.z
			+ dRzdrt(0, 2) * drdomegaz.pR + dRzdrt(0, 3) * drdomegaz.pz;
		dRzdPpot[1].z = dRzdrt(1, 0) * drdomegaz.R + dRzdrt(1, 1) * drdomegaz.z
			+ dRzdrt(1, 2) * drdomegaz.pR + dRzdrt(1, 3) * drdomegaz.pz;
		dRzdPpot[1].pR = dRzdrt(2, 0) * drdomegaz.R + dRzdrt(2, 1) * drdomegaz.z
			+ dRzdrt(2, 2) * drdomegaz.pR + dRzdrt(2, 3) * drdomegaz.pz;
		dRzdPpot[1].pz = dRzdrt(3, 0) * drdomegaz.R + dRzdrt(3, 1) * drdomegaz.z
			+ dRzdrt(3, 2) * drdomegaz.pR + dRzdrt(3, 3) * drdomegaz.pz;
        }
		return Rz;
	}
	coord::PosMomCyl ToyMapIso::from_aaT(const ActionAngles& aaT, DerivAct<coord::Cyl>& dRzdJ) const {
		math::Matrix<double> dRzdrt(4, 4);
		DerivAct<coord::Sph> drdJ;
		coord::PosMomSph rp(Is.aa2pq(aaT, NULL, &drdJ));
		coord::PosMomCyl Rz = PT.map(rp, dRzdrt);
		dRzdJ.dbyJr.R = dRzdrt(0, 0) * drdJ.dbyJr.r + dRzdrt(0, 1) * drdJ.dbyJr.theta
			+ dRzdrt(0, 2) * drdJ.dbyJr.pr + dRzdrt(0, 3) * drdJ.dbyJr.ptheta;
		dRzdJ.dbyJr.z = dRzdrt(1, 0) * drdJ.dbyJr.r + dRzdrt(1, 1) * drdJ.dbyJr.theta
			+ dRzdrt(1, 2) * drdJ.dbyJr.pr + dRzdrt(1, 3) * drdJ.dbyJr.ptheta;
		dRzdJ.dbyJr.phi = drdJ.dbyJr.phi;
		dRzdJ.dbyJr.pR = dRzdrt(2, 0) * drdJ.dbyJr.r + dRzdrt(2, 1) * drdJ.dbyJr.theta
			+ dRzdrt(2, 2) * drdJ.dbyJr.pr + dRzdrt(2, 3) * drdJ.dbyJr.ptheta;
		dRzdJ.dbyJr.pz = dRzdrt(3, 0) * drdJ.dbyJr.r + dRzdrt(3, 1) * drdJ.dbyJr.theta
			+ dRzdrt(3, 2) * drdJ.dbyJr.pr + dRzdrt(3, 3) * drdJ.dbyJr.ptheta;
		dRzdJ.dbyJr.pphi = 0;
		dRzdJ.dbyJz.R = dRzdrt(0, 0) * drdJ.dbyJz.r + dRzdrt(0, 1) * drdJ.dbyJz.theta
			+ dRzdrt(0, 2) * drdJ.dbyJz.pr + dRzdrt(0, 3) * drdJ.dbyJz.ptheta;
		dRzdJ.dbyJz.z = dRzdrt(1, 0) * drdJ.dbyJz.r + dRzdrt(1, 1) * drdJ.dbyJz.theta
			+ dRzdrt(1, 2) * drdJ.dbyJz.pr + dRzdrt(1, 3) * drdJ.dbyJz.ptheta;
		dRzdJ.dbyJz.phi = drdJ.dbyJz.phi;
		dRzdJ.dbyJz.pR = dRzdrt(2, 0) * drdJ.dbyJz.r + dRzdrt(2, 1) * drdJ.dbyJz.theta
			+ dRzdrt(2, 2) * drdJ.dbyJz.pr + dRzdrt(2, 3) * drdJ.dbyJz.ptheta;
		dRzdJ.dbyJz.pz = dRzdrt(3, 0) * drdJ.dbyJz.r + dRzdrt(3, 1) * drdJ.dbyJz.theta
			+ dRzdrt(3, 2) * drdJ.dbyJz.pr + dRzdrt(3, 3) * drdJ.dbyJz.ptheta;
		dRzdJ.dbyJphi.R = dRzdrt(0, 0) * drdJ.dbyJphi.r + dRzdrt(0, 1) * drdJ.dbyJphi.theta
			+ dRzdrt(0, 2) * drdJ.dbyJphi.pr + dRzdrt(0, 3) * drdJ.dbyJphi.ptheta;
		dRzdJ.dbyJphi.z = dRzdrt(1, 0) * drdJ.dbyJphi.r + dRzdrt(1, 1) * drdJ.dbyJphi.theta
			+ dRzdrt(1, 2) * drdJ.dbyJphi.pr + dRzdrt(1, 3) * drdJ.dbyJphi.ptheta;
		dRzdJ.dbyJphi.phi = drdJ.dbyJphi.phi;
		dRzdJ.dbyJphi.pR = dRzdrt(2, 0) * drdJ.dbyJphi.r + dRzdrt(2, 1) * drdJ.dbyJphi.theta
			+ dRzdrt(2, 2) * drdJ.dbyJphi.pr + dRzdrt(2, 3) * drdJ.dbyJphi.ptheta;
		dRzdJ.dbyJphi.pz = dRzdrt(3, 0) * drdJ.dbyJphi.r + dRzdrt(3, 1) * drdJ.dbyJphi.theta
			+ dRzdrt(3, 2) * drdJ.dbyJphi.pr + dRzdrt(3, 3) * drdJ.dbyJphi.ptheta;
		dRzdJ.dbyJphi.pphi = 1;
		return Rz;
	}
	coord::PosMomCyl ToyMapHarm::from_aaT(const ActionAngles& aaT, DerivAct<coord::Cyl>& dRzdJ) const {
		math::Matrix<double> dRzdrt(4, 4);
		DerivAct<coord::Cyl> drdJ;
		coord::PosMomCyl rp(HOs.aa2pq(aaT, NULL, &drdJ));
		coord::PosMomCyl Rz = PT.map(rp, dRzdrt);
		dRzdJ.dbyJr.R = dRzdrt(0, 0) * drdJ.dbyJr.R + dRzdrt(0, 1) * drdJ.dbyJr.z
			+ dRzdrt(0, 2) * drdJ.dbyJr.pR + dRzdrt(0, 3) * drdJ.dbyJr.pz;
		dRzdJ.dbyJr.z = dRzdrt(1, 0) * drdJ.dbyJr.R + dRzdrt(1, 1) * drdJ.dbyJr.z
			+ dRzdrt(1, 2) * drdJ.dbyJr.pR + dRzdrt(1, 3) * drdJ.dbyJr.pz;
		dRzdJ.dbyJr.phi = drdJ.dbyJr.phi;
		dRzdJ.dbyJr.pR = dRzdrt(2, 0) * drdJ.dbyJr.R + dRzdrt(2, 1) * drdJ.dbyJr.z
			+ dRzdrt(2, 2) * drdJ.dbyJr.pR + dRzdrt(2, 3) * drdJ.dbyJr.pz;
		dRzdJ.dbyJr.pz = dRzdrt(3, 0) * drdJ.dbyJr.R + dRzdrt(3, 1) * drdJ.dbyJr.z
			+ dRzdrt(3, 2) * drdJ.dbyJr.pR + dRzdrt(3, 3) * drdJ.dbyJr.pz;
		dRzdJ.dbyJr.pphi = 0;
		dRzdJ.dbyJz.R = dRzdrt(0, 0) * drdJ.dbyJz.R + dRzdrt(0, 1) * drdJ.dbyJz.z
			+ dRzdrt(0, 2) * drdJ.dbyJz.pR + dRzdrt(0, 3) * drdJ.dbyJz.pz;
		dRzdJ.dbyJz.z = dRzdrt(1, 0) * drdJ.dbyJz.R + dRzdrt(1, 1) * drdJ.dbyJz.z
			+ dRzdrt(1, 2) * drdJ.dbyJz.pR + dRzdrt(1, 3) * drdJ.dbyJz.pz;
		dRzdJ.dbyJz.phi = drdJ.dbyJz.phi;
		dRzdJ.dbyJz.pR = dRzdrt(2, 0) * drdJ.dbyJz.R + dRzdrt(2, 1) * drdJ.dbyJz.z
			+ dRzdrt(2, 2) * drdJ.dbyJz.pR + dRzdrt(2, 3) * drdJ.dbyJz.pz;
		dRzdJ.dbyJz.pz = dRzdrt(3, 0) * drdJ.dbyJz.R + dRzdrt(3, 1) * drdJ.dbyJz.z
			+ dRzdrt(3, 2) * drdJ.dbyJz.pR + dRzdrt(3, 3) * drdJ.dbyJz.pz;
		dRzdJ.dbyJphi.R = dRzdrt(0, 0) * drdJ.dbyJphi.R + dRzdrt(0, 1) * drdJ.dbyJphi.z
			+ dRzdrt(0, 2) * drdJ.dbyJphi.pR + dRzdrt(0, 3) * drdJ.dbyJphi.pz;
		dRzdJ.dbyJphi.z = dRzdrt(1, 0) * drdJ.dbyJphi.R + dRzdrt(1, 1) * drdJ.dbyJphi.z
			+ dRzdrt(1, 2) * drdJ.dbyJphi.pR + dRzdrt(1, 3) * drdJ.dbyJphi.pz;
		dRzdJ.dbyJphi.phi = drdJ.dbyJphi.phi;
		dRzdJ.dbyJphi.pR = dRzdrt(2, 0) * drdJ.dbyJphi.R + dRzdrt(2, 1) * drdJ.dbyJphi.z
			+ dRzdrt(2, 2) * drdJ.dbyJphi.pR + dRzdrt(2, 3) * drdJ.dbyJphi.pz;
		dRzdJ.dbyJphi.pz = dRzdrt(3, 0) * drdJ.dbyJphi.R + dRzdrt(3, 1) * drdJ.dbyJphi.z
			+ dRzdrt(3, 2) * drdJ.dbyJphi.pR + dRzdrt(3, 3) * drdJ.dbyJphi.pz;
		dRzdJ.dbyJphi.pphi = 1;
		return Rz;
	}
	coord::PosMomCyl ToyMapIso::from_aaT(const ActionAngles& aaT, DerivAng<coord::Cyl>& dRzdT) const {
		DerivAng<coord::Sph> drdT;
		coord::PosMomSph rp(Is.aa2pq(aaT, NULL, NULL, &drdT));
		math::Matrix<double> dRzdrt(4, 4);
		coord::PosMomCyl Rz(PT.map(rp, dRzdrt));
		dRzdT.dbythetar.R = dRzdrt(0, 0) * drdT.dbythetar.r + dRzdrt(0, 1) * drdT.dbythetar.theta
			+ dRzdrt(0, 2) * drdT.dbythetar.pr + dRzdrt(0, 3) * drdT.dbythetar.ptheta;
		dRzdT.dbythetar.z = dRzdrt(1, 0) * drdT.dbythetar.r + dRzdrt(1, 1) * drdT.dbythetar.theta
			+ dRzdrt(1, 2) * drdT.dbythetar.pr + dRzdrt(1, 3) * drdT.dbythetar.ptheta;
		dRzdT.dbythetar.phi = drdT.dbythetar.phi;
		dRzdT.dbythetar.pR = dRzdrt(2, 0) * drdT.dbythetar.r + dRzdrt(2, 1) * drdT.dbythetar.theta
			+ dRzdrt(2, 2) * drdT.dbythetar.pr + dRzdrt(2, 3) * drdT.dbythetar.ptheta;
		dRzdT.dbythetar.pz = dRzdrt(3, 0) * drdT.dbythetar.r + dRzdrt(3, 1) * drdT.dbythetar.theta
			+ dRzdrt(3, 2) * drdT.dbythetar.pr + dRzdrt(3, 3) * drdT.dbythetar.ptheta;
		dRzdT.dbythetar.pphi = 0;
		dRzdT.dbythetaz.R = dRzdrt(0, 0) * drdT.dbythetaz.r + dRzdrt(0, 1) * drdT.dbythetaz.theta
			+ dRzdrt(0, 2) * drdT.dbythetaz.pr + dRzdrt(0, 3) * drdT.dbythetaz.ptheta;
		dRzdT.dbythetaz.z = dRzdrt(1, 0) * drdT.dbythetaz.r + dRzdrt(1, 1) * drdT.dbythetaz.theta
			+ dRzdrt(1, 2) * drdT.dbythetaz.pr + dRzdrt(1, 3) * drdT.dbythetaz.ptheta;
		dRzdT.dbythetaz.phi = drdT.dbythetaz.phi;
		dRzdT.dbythetaz.pR = dRzdrt(2, 0) * drdT.dbythetaz.r + dRzdrt(2, 1) * drdT.dbythetaz.theta
			+ dRzdrt(2, 2) * drdT.dbythetaz.pr + dRzdrt(2, 3) * drdT.dbythetaz.ptheta;
		dRzdT.dbythetaz.pz = dRzdrt(3, 0) * drdT.dbythetaz.r + dRzdrt(3, 1) * drdT.dbythetaz.theta
			+ dRzdrt(3, 2) * drdT.dbythetaz.pr + dRzdrt(3, 3) * drdT.dbythetaz.ptheta;
		dRzdT.dbythetaz.pphi = 0;
		dRzdT.dbythetaphi.R = dRzdrt(0, 0) * drdT.dbythetaphi.r + dRzdrt(0, 1) * drdT.dbythetaphi.theta
			+ dRzdrt(0, 2) * drdT.dbythetaphi.pr + dRzdrt(0, 3) * drdT.dbythetaphi.ptheta;
		dRzdT.dbythetaphi.z = dRzdrt(1, 0) * drdT.dbythetaphi.r + dRzdrt(1, 1) * drdT.dbythetaphi.theta
			+ dRzdrt(1, 2) * drdT.dbythetaphi.pr + dRzdrt(1, 3) * drdT.dbythetaphi.ptheta;
		dRzdT.dbythetaphi.phi = 1;
		dRzdT.dbythetaphi.pR = dRzdrt(2, 0) * drdT.dbythetaphi.r + dRzdrt(2, 1) * drdT.dbythetaphi.theta
			+ dRzdrt(2, 2) * drdT.dbythetaphi.pr + dRzdrt(2, 3) * drdT.dbythetaphi.ptheta;
		dRzdT.dbythetaphi.pz = dRzdrt(3, 0) * drdT.dbythetaphi.r + dRzdrt(3, 1) * drdT.dbythetaphi.theta
			+ dRzdrt(3, 2) * drdT.dbythetaphi.pr + dRzdrt(3, 3) * drdT.dbythetaphi.ptheta;
		dRzdT.dbythetaphi.pphi = 0;
		return Rz;
	}
	coord::PosMomCyl ToyMapHarm::from_aaT(const ActionAngles& aaT, DerivAng<coord::Cyl>& dRzdT) const {
		DerivAng<coord::Cyl> drdT;
		coord::PosMomCyl rp(HOs.aa2pq(aaT, NULL, NULL, &drdT));
		math::Matrix<double> dRzdrt(4, 4);
		coord::PosMomCyl Rz(PT.map(rp, dRzdrt));
		dRzdT.dbythetar.R = dRzdrt(0, 0) * drdT.dbythetar.R + dRzdrt(0, 1) * drdT.dbythetar.z
			+ dRzdrt(0, 2) * drdT.dbythetar.pR + dRzdrt(0, 3) * drdT.dbythetar.pz;
		dRzdT.dbythetar.z = dRzdrt(1, 0) * drdT.dbythetar.R + dRzdrt(1, 1) * drdT.dbythetar.z
			+ dRzdrt(1, 2) * drdT.dbythetar.pR + dRzdrt(1, 3) * drdT.dbythetar.pz;
		dRzdT.dbythetar.phi = drdT.dbythetar.phi;
		dRzdT.dbythetar.pR = dRzdrt(2, 0) * drdT.dbythetar.R + dRzdrt(2, 1) * drdT.dbythetar.z
			+ dRzdrt(2, 2) * drdT.dbythetar.pR + dRzdrt(2, 3) * drdT.dbythetar.pz;
		dRzdT.dbythetar.pz = dRzdrt(3, 0) * drdT.dbythetar.R + dRzdrt(3, 1) * drdT.dbythetar.z
			+ dRzdrt(3, 2) * drdT.dbythetar.pR + dRzdrt(3, 3) * drdT.dbythetar.pz;
		dRzdT.dbythetar.pphi = 0;
		dRzdT.dbythetaz.R = dRzdrt(0, 0) * drdT.dbythetaz.R + dRzdrt(0, 1) * drdT.dbythetaz.z
			+ dRzdrt(0, 2) * drdT.dbythetaz.pR + dRzdrt(0, 3) * drdT.dbythetaz.pz;
		dRzdT.dbythetaz.z = dRzdrt(1, 0) * drdT.dbythetaz.R + dRzdrt(1, 1) * drdT.dbythetaz.z
			+ dRzdrt(1, 2) * drdT.dbythetaz.pR + dRzdrt(1, 3) * drdT.dbythetaz.pz;
		dRzdT.dbythetaz.phi = drdT.dbythetaz.phi;
		dRzdT.dbythetaz.pR = dRzdrt(2, 0) * drdT.dbythetaz.R + dRzdrt(2, 1) * drdT.dbythetaz.z
			+ dRzdrt(2, 2) * drdT.dbythetaz.pR + dRzdrt(2, 3) * drdT.dbythetaz.pz;
		dRzdT.dbythetaz.pz = dRzdrt(3, 0) * drdT.dbythetaz.R + dRzdrt(3, 1) * drdT.dbythetaz.z
			+ dRzdrt(3, 2) * drdT.dbythetaz.pR + dRzdrt(3, 3) * drdT.dbythetaz.pz;
		dRzdT.dbythetaz.pphi = 0;
		dRzdT.dbythetaphi.R = dRzdrt(0, 0) * drdT.dbythetaphi.R + dRzdrt(0, 1) * drdT.dbythetaphi.z
			+ dRzdrt(0, 2) * drdT.dbythetaphi.pR + dRzdrt(0, 3) * drdT.dbythetaphi.pz;
		dRzdT.dbythetaphi.z = dRzdrt(1, 0) * drdT.dbythetaphi.R + dRzdrt(1, 1) * drdT.dbythetaphi.z
			+ dRzdrt(1, 2) * drdT.dbythetaphi.pR + dRzdrt(1, 3) * drdT.dbythetaphi.pz;
		dRzdT.dbythetaphi.phi = 1;
		dRzdT.dbythetaphi.pR = dRzdrt(2, 0) * drdT.dbythetaphi.R + dRzdrt(2, 1) * drdT.dbythetaphi.z
			+ dRzdrt(2, 2) * drdT.dbythetaphi.pR + dRzdrt(2, 3) * drdT.dbythetaphi.pR;
		dRzdT.dbythetaphi.pz = dRzdrt(3, 0) * drdT.dbythetaphi.R + dRzdrt(3, 1) * drdT.dbythetaphi.z
			+ dRzdrt(3, 2) * drdT.dbythetaphi.pR + dRzdrt(3, 3) * drdT.dbythetaphi.pz;
		dRzdT.dbythetaphi.pphi = 0;
		return Rz;
	}
}