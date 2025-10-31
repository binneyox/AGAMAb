/** \file    df_halo.h
    \brief   Distribution functions for the spheroidal component (halo)
    \author  Eugene Vasiliev, James Binney
    \date    2015-2019
*/
#pragma once
#include "df_base.h"
#include "units.h"
#include "math_core.h"
#include "math_ode.h"
#include "potential_utils.h"
#include "potential_factory.h"
#include "actions_newtorus.h"
#include <fstream>
#include <cmath>

namespace df{

/// \name   Classes for action-based double power-law distribution function (DF)
///@{

/// Parameters that describe a double power law distribution function.
struct DoublePowerLawParam{
	double
			norm,      ///< normalization factor with the dimension of mass
			mass,	   ///< desired total mass
			J0,        ///< break action (defines the transition between inner and outer regions)
			Jcutoff,   ///< cutoff action (sets exponential suppression at J>Jcutoff, 0 to disable)
			Jphi0,     ///< controls the steepness of rotation and the size of non-rotating core
			Jcore,     ///< central core size for a Cole&Binney-type modified double-power-law halo
			epsilonJ,    ///< parameter setting size of region where df/dJz=df/dJphi et 
			slopeIn,   ///< power-law index for actions below the break action (Gamma)
			slopeOut,  ///< power-law index for actions above the break action (Beta)
			steepness, ///< steepness of the transition between two asymptotic regimes (eta)
			cutoffStrength, ///< steepness of exponential suppression at J>Jcutoff (zeta)
			coefJrIn,  ///< contribution of radial   action to h(J), controlling anisotropy below J_0 (h_r)
			coefJzIn,  ///< contribution of vertical action to h(J), controlling anisotropy below J_0 (h_z)
			coefJrOut, ///< contribution of radial   action to g(J), controlling anisotropy above J_0 (g_r)
			coefJzOut, ///< contribution of vertical action to g(J), controlling anisotropy above J_0 (g_z)
			rotFrac;   ///< relative amplitude of the odd-Jphi component (-1 to 1, 0 means no rotation)
	std::string Fname;
	DoublePowerLawParam() :  ///< set default values for all fields (NAN means that it must be set manually)
	    norm(NAN), J0(NAN), Jcutoff(1e6), Jphi0(NAN), Jcore(0), epsilonJ(NAN), 
	    slopeIn(NAN), slopeOut(NAN), steepness(1), cutoffStrength(2),
	    coefJrIn(1.5), coefJzIn(1), coefJrOut(1.3), coefJzOut(1), rotFrac(0) {}
};

// Parameters of DSph DF introduced by Pascale+
struct DwarfSpheroidParam{
	double
			norm,
			mass,
			J0,
			Jphi0,
			epsilonJ,
			alpha,
			coefJr,
			coefJz,
			rotFrac;
	std::string Fname;
	DwarfSpheroidParam() :
	    norm(1), J0(NAN), Jphi0(NAN), epsilonJ(NAN), alpha(.25),
	    coefJr(1.5), coefJz(1), rotFrac(0) {} 
};

/** General double power-law model.
    The distribution function is given by
    \f$  f(J) = norm / (2\pi J_0)^3
         (1 + (J_0 /h(J))^\eta )^{\Gamma / \eta}
         (1 + (g(J)/ J_0)^\eta )^{-B / \eta }
         \exp[ - (g(J) / J_{cutoff})^\zeta ]
         ( [J_{core} / h(J)]^2 - \beta J_{core} / h(J) + 1)^{\Gamma/2}   \f$,  where
    \f$  g(J) = g_r J_r + g_z J_z + g_\phi |J_\phi|  \f$,
    \f$  h(J) = h_r J_r + h_z J_z + h_\phi |J_\phi|  \f$.
    Gamma is the power-law slope of DF at small J (slopeIn), and Beta -- at large J (slopeOut),
    the transition occurs around J=J0, and its steepness is adjusted by the parameter eta.
    h_r, h_z and h_phi control the anisotropy of the DF at small J (their sum is always taken
    to be equal to 3, so that there are two free parameters -- coefJrIn = h_r, coefJzIn = h_z),
    and g_r, g_z, g_phi do the same for large J (coefJrOut = g_r, coefJzOut = g_z).
    Jcutoff is the threshold of an optional exponential suppression, and zeta measures its strength.
    Jcore is the size of the central core (if nonzero, f(J) tends to a constant limit as J --> 0
    even when the power-law slope Gamma is positive), and the auxiliary coefficient beta
    is assigned automatically from the requirement that the introduction of the core (almost)
    doesn't change the overall normalization (eq.5 in Cole&Binney 2017).
*/
class EXP DoublePowerLaw: public BaseDistributionFunction{
	const DoublePowerLawParam par;  ///< parameters of DF
	const double beta;              ///< auxiliary coefficient for the case of a central core
	public:
    /** Create an instance of double-power-law distribution function with given parameters
        \param[in] params  are the parameters of DF
        \throws std::invalid_argument exception if parameters are nonsense
    */
		DoublePowerLaw(const DoublePowerLawParam &params);

    /** return value of DF for the given set of actions.
        \param[in] J are the actions  */
		virtual double value(const actions::Actions &J, const double Jzcrit) const;
};

/// Parameters of f(J_r,L) that delivers anisotropic double-power-law
/// spheres by solving an ODE for f
struct NewOxfordParam{
	double
			norm,
			mass,
			J0,
			Jcutoff,
			Jphi0,
			Jcore,
			slopeIn,
			slopeOut,
			steepness,
			beta,
			rotFrac,
			cutoffStrength;
	std::string Fname;
	NewOxfordParam() :  ///< set default values for all fields (NAN means that it must be set manually)
	    norm(1), mass(NAN), J0(NAN), Jcutoff(1e6), Jphi0(0), Jcore(0),
	    slopeIn(NAN), slopeOut(NAN), steepness(1), beta(0),
	    rotFrac(0), cutoffStrength(2) {}
};


/// Creates anisotropic flattened double-powe-law systems
class EXP NewDoublePowerLaw: public BaseDistributionFunction{
	const double epsilonJ;
	const df::DoublePowerLaw DF0;
	double wt(const actions::Actions& J) const;
	public:
		NewDoublePowerLaw(const DoublePowerLawParam& params,
				  const potential::BasePotential& pot);
		virtual double value(const actions::Actions& J, const double Jzcrit) const;
};

// Creates anisotropic spheres
class EXP NewOxford : public BaseDistributionFunction{
	const NewOxfordParam par;  ///< parameters of DF
	const potential::Interpolator freq;  ///< interface providing the epicyclic frequencies and Rcirc
	public:
		NewOxford(const NewOxfordParam& _params,
			  const potential::Interpolator& _freq) :
		    par(_params), freq(_freq) {}
		virtual double value(const actions::Actions &, const double Jzcrit) const;
		double h(const actions::Actions&) const;
};

//Creates Pascale-type dwarf spheroids: use by calling NewDwarfSpheroid
class EXP DwarfSpheroid: public BaseDistributionFunction{
	const DwarfSpheroidParam par;
	public:
		DwarfSpheroid(const DwarfSpheroidParam& _par);
		virtual double value(const actions::Actions& J, const double Jzcrit) const;
};

//Creates anistropic, flattened dSph systems
class EXP NewDwarfSpheroid: public BaseDistributionFunction{
	const DwarfSpheroidParam par;
	const df::DwarfSpheroid DF0;
	const double epsilonJ;
	double wt(const actions::Actions& J) const;
	public:
		NewDwarfSpheroid(const DwarfSpheroidParam& _par,
				 const potential::BasePotential& pot);
		virtual double value(const actions::Actions& J, const double Jzcrit) const;
};

struct PlummerParam{
	double mass;
	double scaleRadius;
	double scaleAction;
	double mu, nu;
	PlummerParam(): mass(1), scaleRadius(1), scaleAction(1), mu(0), nu(0){}
};
class EXP PlummerDF : public BaseDistributionFunction{
	const PlummerParam par;
	double norm,Etop,Ebot,Jrtop,cLtop,Jrbot,cLbot;
	math::LinearInterpolator jrmax;
	math::LinearInterpolator cls;
	public:
		PlummerDF(const PlummerParam&);
		virtual double value(const actions::Actions& J, const double Jzcrit) const;
		virtual void write_params(std::ofstream&, const units::InternalUnits&) const;
};

struct IsochroneParam{
	double mass;
	double scaleRadius;
	double mu, nu;
	IsochroneParam() : mass(1), scaleRadius(1), mu(0), nu(0) {}
};
class EXP IsochroneDF : public BaseDistributionFunction{
	const IsochroneParam par;
	double norm;
	public:
		IsochroneDF(const IsochroneParam& params);
		virtual double value(const actions::Actions& J, const double Jzcrit) const;
		virtual void write_params(std::ofstream& stream,const units::InternalUnits& intUnits) const;
};

///@}
}  // namespace df
