/** \file    potential_composite.h
    \brief   Composite density and potential classes
    \author  Eugene Vasiliev
    \date    2014-2015
*/
#pragma once
#include "potential_base.h"
#include "potential_interpolators.h"
#include "smart.h"
#include <vector>

namespace potential{

/** A trivial collection of several density objects */
class EXP CompositeDensity: public BaseDensity{
public:
    /** construct from the provided array of components */
    CompositeDensity(const std::vector<PtrDensity>& _components);

    /** provides the 'least common denominator' for the symmetry degree */
    virtual coord::SymmetryType symmetry() const;

    virtual const char* name() const { return myName(); }
    static const char* myName() { static const char* text = "CompositeDensity"; return text; }

    unsigned int size() const { return components.size(); }
    PtrDensity component(unsigned int index) const { return components.at(index); }
private:
    std::vector<PtrDensity> components;
    virtual double densityCar(const coord::PosCar &pos) const;
    virtual double densityCyl(const coord::PosCyl &pos) const;
    virtual double densitySph(const coord::PosSph &pos) const;
};

/** A trivial collection of several potential objects, evaluated in cylindrical coordinates */
class EXP CompositeCyl: public BasePotentialCyl{
	private:
		std::vector<PtrPotential> components;
		PtrShellInterpolator PtrShellI;
		PtrPolarInterpolator PtrPolarI;
		virtual void evalCyl(const coord::PosCyl &pos,
				     double* potential, coord::GradCyl* deriv, coord::HessCyl* deriv2) const;
	public:
    /** construct from the provided array of components */
		CompositeCyl(const std::vector<PtrPotential>& _components,
			     PtrShellInterpolator _PtrShellI=NULL, PtrPolarInterpolator _PtrPolarI=NULL) :
		BasePotentialCyl(), components(_components), PtrShellI(_PtrShellI), PtrPolarI(_PtrPolarI)
		{
			if(_components.empty())
				throw std::invalid_argument("List of potential components cannot be empty");
		}

    /** provides the 'least common denominator' for the symmetry degree */
		virtual coord::SymmetryType symmetry() const;
		virtual const char* name() const { return myName(); };
		static const char* myName() { static const char* text = "CompositePotential"; return text; }

		unsigned int size() const { return components.size(); }
		PtrPotential component(unsigned int index) const { return components.at(index); }
		//Functions that extract data from ShellInterpolator
		virtual void getRshDelta(const double L, const double Xip,
					 double& Rsh, double& Delta) const{
			if(PtrShellI)
				PtrShellI->getRshDelta(L, Xip, Rsh, Delta);
			else
				printf("Error PtrShellI is NULL\n");
		}
		virtual void getRshDelta(const double E, const double Xip,
					 const double invPhi0, double& Rsh, double& Delta) const{
			if(PtrShellI)
				PtrShellI->getRshDelta(E, Xip, invPhi0, Rsh, Delta);
			else
				printf("Error PtrShellI is NULL\n");
		}
		virtual double getRsh(const double L, const double Xip) const{
			double Rsh, Delta;
			getRshDelta(L, Xip, Rsh, Delta);
			return Rsh;
		}
		virtual double getRsh(const double E, const double Xip, const double invPhi0) const{
			double Rsh, Delta;
			getRshDelta(E, Xip, invPhi0, Rsh, Delta);
			return Rsh;
		}
		virtual double getDelta(const double L, const double Xip) const{
			double Rsh, Delta;
			getRshDelta(L, Xip, Rsh, Delta);
			return Delta;
		}
		virtual double getDelta(const double E, const double Xip, const double invPhi0) const{
			double Rsh, Delta;
			getRshDelta(E, Xip, invPhi0, Rsh, Delta);
			return Delta;
		}
		//Now extraction from PolarInterpolator
		virtual double getJzcrit(const double Jf) const{
			if(PtrPolarI)
				return PtrPolarI->getJz(Jf);
			else{
				printf("PtrPolar is NULL\n");
				return 0;
			}
		}
		virtual void getBestDeltaI3(const double E, const double invPhi0,
					    double& Delta, double& I3) const{
			if(PtrPolarI)
				PtrPolarI->getDeltaI3(E, invPhi0, Delta, I3);
			else
				printf("Error PtrPolarI is NULL\n");
		}
		virtual double getBestDelta(const double E, const double invPhi0) const{
			if(PtrPolarI)
				return PtrPolarI->getDelta(E, invPhi0);
			else{
				printf("Error PtrPolarI is NULL\n");
				return 0;
			}
		}
		virtual double getI3(const double E, const double invPhi0) const{
			if(PtrPolarI)
				return PtrPolarI->getI3(E, invPhi0);
			else{
				printf("Error PtrPolarI is NULL\n");
				return 0;
			}
		}
		virtual std::pair<double,double> getXiChi(const double E,const double Xip, const double invPhi0) const{
			if(PtrShellI)
				return PtrShellI->getXiChi(E, Xip, invPhi0);
			else{
				printf("Error PtrPolarI is NULL\n");
				return std::make_pair(0,0);
			}
		}


};

}  // namespace potential
