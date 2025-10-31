#include "df_factory.h"
#include "df_disk.h"
#include "df_halo.h"
#include "utils.h"
#include <cassert>
#include <stdexcept>

namespace df {

EXP DoublePowerLawParam parseDoublePowerLawParam(
	const utils::KeyValueMap& kvmap,
	const units::ExternalUnits& conv)
{
	DoublePowerLawParam par;
	par.norm      = kvmap.getDouble("norm",      par.norm)    * conv.massUnit;
	par.mass      = kvmap.getDouble("mass",      par.mass)    * conv.massUnit;
	par.J0        = kvmap.getDouble("J0",        par.J0)      * conv.lengthUnit * conv.velocityUnit;
	par.Jcutoff   = kvmap.getDouble("Jcutoff",   par.Jcutoff) * conv.lengthUnit * conv.velocityUnit;
	par.Jphi0     = kvmap.getDouble("Jphi0",     par.Jphi0)   * conv.lengthUnit * conv.velocityUnit;
	par.Jcore     = kvmap.getDouble("Jcore",     par.Jcore)   * conv.lengthUnit * conv.velocityUnit;
	par.epsilonJ  = kvmap.getDouble("epsilonJ",  par.epsilonJ)  * conv.lengthUnit * conv.velocityUnit;
	par.slopeIn   = kvmap.getDouble("slopeIn",   par.slopeIn);
	par.slopeOut  = kvmap.getDouble("slopeOut",  par.slopeOut);
	par.steepness = kvmap.getDouble("steepness", par.steepness);
	par.coefJrIn  = kvmap.getDouble("coefJrIn",  par.coefJrIn);
	par.coefJzIn  = kvmap.getDouble("coefJzIn",  par.coefJzIn);
	par.coefJrOut = kvmap.getDouble("coefJrOut", par.coefJrOut);
	par.coefJzOut = kvmap.getDouble("coefJzOut", par.coefJzOut);
	par.rotFrac   = kvmap.getDouble("rotFrac",   par.rotFrac);
	par.cutoffStrength = kvmap.getDouble("cutoffStrength", par.cutoffStrength);
	par.Fname     = kvmap.getString("PopFile", par.Fname);
	if(std::isnan(par.J0+par.epsilonJ+par.slopeIn+par.slopeOut))
		printf("parseDoublePowerLawParam: vital parametr not set\n");
	return par;
}

EXP DwarfSpheroidParam parseDwarfSpheroidParam(
	const utils::KeyValueMap& kvmap,
	const units::ExternalUnits& conv)
{
	DwarfSpheroidParam par;
	par.norm	= kvmap.getDouble("norm", par.norm)	* conv.massUnit;
	par.mass        = kvmap.getDouble("mass",      par.mass)    * conv.massUnit;
	par.J0		= kvmap.getDouble("J0", par.J0)		* conv.lengthUnit * conv.velocityUnit;
	par.Jphi0	= kvmap.getDouble("Jphi0", par.Jphi0)	* conv.lengthUnit * conv.velocityUnit;
	par.epsilonJ	= kvmap.getDouble("epsilonJ", par.epsilonJ) * conv.lengthUnit * conv.velocityUnit;
	par.alpha	= kvmap.getDouble("alpha", par.alpha);
	par.coefJr      = kvmap.getDouble("coefJr",  par.coefJr);
	par.coefJz      = kvmap.getDouble("coefJz",  par.coefJz);
	par.rotFrac	= kvmap.getDouble("rotFrac", par.rotFrac);
	par.Fname       = kvmap.getString("PopFile", par.Fname);
	if(std::isnan(par.J0+par.epsilonJ))
		printf("parseDwarfSpheroidParam: vital parametr not set\n");
	return par;
}

EXP QuasiIsothermalParam parseQuasiIsothermalParam(
    const utils::KeyValueMap& kvmap,
    const units::ExternalUnits& conv)
{
    QuasiIsothermalParam par;
    par.Sigma0  = kvmap.getDouble("Sigma0",  par.Sigma0)  * conv.massUnit / pow_2(conv.lengthUnit);
    par.Rdisk   = kvmap.getDouble("Rdisk",   par.Rdisk)   * conv.lengthUnit;
    par.Hdisk   = kvmap.getDouble("Hdisk",   par.Hdisk)   * conv.lengthUnit;
    par.sigmar0 = kvmap.getDouble("sigmar0", par.sigmar0) * conv.velocityUnit;
    par.sigmaz0 = kvmap.getDouble("sigmaz0", par.sigmaz0) * conv.velocityUnit;
    par.sigmamin= kvmap.getDouble("sigmamin",par.sigmamin)* conv.velocityUnit;
    par.Rsigmar = kvmap.getDouble("Rsigmar", par.Rsigmar) * conv.lengthUnit;
    par.Rsigmaz = kvmap.getDouble("Rsigmaz", par.Rsigmaz) * conv.lengthUnit;
    par.coefJr  = kvmap.getDouble("coefJr",  par.coefJr);
    par.coefJz  = kvmap.getDouble("coefJz",  par.coefJz);
    par.Jmin    = kvmap.getDouble("Jmin",    par.Jmin)    * conv.lengthUnit * conv.velocityUnit;
    par.beta    = kvmap.getDouble("beta",    par.beta);
    par.Tsfr    = kvmap.getDouble("Tsfr",    par.Tsfr);  // dimensionless! in units of galaxy age
    par.sigmabirth = kvmap.getDouble("sigmabirth", par.sigmabirth);  // dimensionless ratio
    return par;
}
/*
ExponentialParam parseExponentialParam(
				       const utils::KeyValueMap& kvmap,
				       const units::ExternalUnits& conv)
{
	ExponentialParam par;
	par.norm   = kvmap.getDouble("norm",   par.norm)   * conv.massUnit;
	par.Jr0    = kvmap.getDouble("Jr0",    par.Jr0)    * conv.lengthUnit * conv.velocityUnit;
	par.Jz0    = kvmap.getDouble("Jz0",    par.Jz0)    * conv.lengthUnit * conv.velocityUnit;
	par.Jphi0  = kvmap.getDouble("Jphi0",  par.Jphi0)  * conv.lengthUnit * conv.velocityUnit;
	par.addJden= kvmap.getDouble("addJden")* conv.lengthUnit * conv.velocityUnit;
	par.addJvel= kvmap.getDouble("addJvel")* conv.lengthUnit * conv.velocityUnit;
	par.coefJr = kvmap.getDouble("coefJr", par.coefJr);
	par.coefJz = kvmap.getDouble("coefJz", par.coefJz);
	par.beta   = kvmap.getDouble("beta",   par.beta);
	par.Tsfr   = kvmap.getDouble("Tsfr",   par.Tsfr);  // dimensionless! in units of Hubble time
	par.sigmabirth = kvmap.getDouble("sigmabirth", par.sigmabirth);  // dimensionless ratio
	return par;
}*/

EXP ExponentialParam parseExponentialParam(
				       const utils::KeyValueMap& kvmap,
				       const units::ExternalUnits& conv)
{
	ExponentialParam par;
	par.norm   = kvmap.getDouble("norm",	par.norm)   * conv.massUnit;
	par.mass   = kvmap.getDouble("mass",	par.mass)   * conv.massUnit;
	par.Jr0    = kvmap.getDouble("Jr0",	par.Jr0)    * conv.lengthUnit * conv.velocityUnit;
	par.Jz0    = kvmap.getDouble("Jz0",	par.Jz0)    * conv.lengthUnit * conv.velocityUnit;
	par.Jphi0  = kvmap.getDouble("Jphi0",	par.Jphi0)  * conv.lengthUnit * conv.velocityUnit;
	par.pr     = kvmap.getDouble("pr",	par.pr);
	par.pz     = kvmap.getDouble("pz",	par.pz);
	par.addJden= kvmap.getDouble("addJden",	par.addJden)* conv.lengthUnit * conv.velocityUnit;
	par.addJvel= kvmap.getDouble("addJvel",	par.addJvel)* conv.lengthUnit * conv.velocityUnit;
	par.Fname  = kvmap.getString("PopFile",	par.Fname);
	return par;
}

EXP taperExpParam parsetaperExpParam(
				  const utils::KeyValueMap& kvmap,
				  const units::ExternalUnits& conv)
{
	taperExpParam par;
	par.norm   = kvmap.getDouble("norm",	par.norm)   * conv.massUnit;
	par.mass   = kvmap.getDouble("mass",	par.mass)   * conv.massUnit;
	par.Jr0    = kvmap.getDouble("Jr0",	par.Jr0)    * conv.lengthUnit * conv.velocityUnit;
	par.Jz0    = kvmap.getDouble("Jz0",	par.Jz0)    * conv.lengthUnit * conv.velocityUnit;
	par.Jtaper = kvmap.getDouble("Jtaper",	par.Jtaper) * conv.lengthUnit * conv.velocityUnit;
	par.Jtrans = kvmap.getDouble("Jtrans",	par.Jtrans) * conv.lengthUnit * conv.velocityUnit;
	par.Jcut   = kvmap.getDouble("Jcut",	par.Jcut) * conv.lengthUnit * conv.velocityUnit;
	par.Delta  = kvmap.getDouble("Delta",	par.Delta) * conv.lengthUnit * conv.velocityUnit;
	par.Jphi0  = kvmap.getDouble("Jphi0",	par.Jphi0)  * conv.lengthUnit * conv.velocityUnit;
	par.pr	   = kvmap.getDouble("pr",	par.pr);
	par.pz	   = kvmap.getDouble("pz",	par.pz);
	par.addJden= kvmap.getDouble("addJden",	par.addJden)* conv.lengthUnit * conv.velocityUnit;
	par.addJvel= kvmap.getDouble("addJvel", par.addJvel)* conv.lengthUnit * conv.velocityUnit;
	par.Fname  = kvmap.getString("PopFile", par.Fname);
	return par;
}
EXP PlummerParam parsePlummerParam(
				   const utils::KeyValueMap& kvmap,
				   const units::ExternalUnits& conv)
{
	PlummerParam par;
	par.mass        = kvmap.getDouble("mass"       ,par.mass)       * conv.massUnit;
	par.scaleRadius = kvmap.getDouble("scaleRadius",par.scaleRadius)* conv.lengthUnit;
	par.scaleAction = kvmap.getDouble("scaleAction",par.scaleAction)* conv.lengthUnit*conv.velocityUnit;
	par.mu          = kvmap.getDouble("mu"         ,par.mu);
	par.nu          = kvmap.getDouble("nu"         ,par.nu);
	printf("loading Plummer M %f b %f L0 %f mu %f nu %f\n",
	       par.mass,par.scaleRadius,par.scaleAction,par.mu,par.nu);
	return par;
}
EXP IsochroneParam parseIsochroneParam(
				       const utils::KeyValueMap& kvmap,
				       const units::ExternalUnits& conv)
{
	IsochroneParam par;
	par.mass        = kvmap.getDouble("mass"       ,par.mass)       * conv.massUnit;
	par.scaleRadius = kvmap.getDouble("scaleRadius",par.scaleRadius)* conv.lengthUnit;
	par.mu          = kvmap.getDouble("mu"         ,par.mu);
	par.nu          = kvmap.getDouble("nu"         ,par.nu);
	printf("loading Isochrone M %f b %f mu %f nu %f\n",par.mass,par.scaleRadius,par.mu,par.nu);
	return par;
}

EXP NewOxfordParam parseNewOxfordParam(
				  const utils::KeyValueMap& kvmap,
				  const units::ExternalUnits& conv)
{
	NewOxfordParam par;
	par.norm      = kvmap.getDouble("norm", par.norm)    * conv.massUnit;
	par.mass      = kvmap.getDouble("mass", par.mass)    * conv.massUnit;
	par.J0        = kvmap.getDouble("J0", par.J0)        * conv.lengthUnit * conv.velocityUnit;
	par.Jcutoff   = kvmap.getDouble("Jcutoff", par.Jcutoff) * conv.lengthUnit * conv.velocityUnit;
	par.Jphi0     = kvmap.getDouble("Jphi0", par.Jcutoff)* conv.lengthUnit * conv.velocityUnit;
	par.Jcore     = kvmap.getDouble("Jcore", par.Jcore)  * conv.lengthUnit * conv.velocityUnit;
	par.slopeIn   = kvmap.getDouble("slopeIn",   par.slopeIn);
	par.slopeOut  = kvmap.getDouble("slopeOut",  par.slopeOut);
	par.steepness = kvmap.getDouble("steepness", par.steepness);
	par.beta     = kvmap.getDouble("beta",	     par.beta);
	par.rotFrac   = kvmap.getDouble("rotFrac",   par.rotFrac);
	par.cutoffStrength = kvmap.getDouble("cutoffStrength", par.cutoffStrength);
	par.Fname  = kvmap.getString("PopFile", par.Fname);
	return par;
}


inline void checkNonzero(const potential::BasePotential* potential, const std::string& type)
{
    if(potential == NULL)
        throw std::invalid_argument("Need an instance of potential to initialize "+type+" DF");
}

EXP PtrDistributionFunction createDistributionFunction(
    const utils::KeyValueMap& kvmap,
    const potential::BasePotential* potential,
    const potential::BaseDensity* density,
    const units::ExternalUnits& converter)
{
	std::string type = kvmap.getString("type");
    // for some DF types, there are two alternative ways of specifying the normalization:
    // either directly as norm, Sigma0, etc., or as the total mass, from which the norm is computed
    // by creating a temporary instance of a corresponding DF class, and computing its mass
    double mass = kvmap.getDouble("mass", NAN) * converter.massUnit;
    if(utils::stringsEqual(type, "DoublePowerLaw")) {
	    DoublePowerLawParam par = parseDoublePowerLawParam(kvmap, converter);
	    if(mass>0) {
		    par.norm = 1.0;
		    par.norm = mass / DoublePowerLaw(par).totalMass();
	    }
	    return PtrDistributionFunction(new DoublePowerLaw(par));
    }
    else if(utils::stringsEqual(type, "NewDoublePowerLaw")) {
	    DoublePowerLawParam par = parseDoublePowerLawParam(kvmap, converter);
	    if(mass>0) {
		    par.norm = 1.0;
		    par.norm = mass / DoublePowerLaw(par).totalMass();
	    }
	    return PtrDistributionFunction(new NewDoublePowerLaw(par, *potential));
    }
    else if(utils::stringsEqual(type, "DwarfSpheroid")) {
	    DwarfSpheroidParam par = parseDwarfSpheroidParam(kvmap, converter);
	    if(mass>0) {
		    par.norm = 1;
		    par.norm = mass/DwarfSpheroid(par).totalMass();
	    }
	    return PtrDistributionFunction(new DwarfSpheroid(par));
    }
    else if(utils::stringsEqual(type, "Exponential")) {
	    ExponentialParam par = parseExponentialParam(kvmap, converter);
	    if(mass>0) {
		    par.norm = 1.0;
		    par.norm = mass / Exponential(par).totalMass();
	    }
	    return PtrDistributionFunction(new Exponential(par));
    }
    else if(utils::stringsEqual(type, "taperExp")) {
	    taperExpParam par = parsetaperExpParam(kvmap, converter);
	    if(mass>0) {
		    par.norm = 1.0;
		    par.norm = mass / taperExp(par).totalMass();
	    }
	    return PtrDistributionFunction(new taperExp(par));
    }
    else if(utils::stringsEqual(type, "QuasiIsothermal")) {
	    checkNonzero(potential, type);
	    potential::Interpolator pot_interp(*potential);
        QuasiIsothermalParam par = parseQuasiIsothermalParam(kvmap, converter);
        if(mass>0) {
            par.Sigma0 = 1.0;
            par.Sigma0 = mass / QuasiIsothermal(par, pot_interp).totalMass();
        }
        return PtrDistributionFunction(new QuasiIsothermal(par, pot_interp));
    }
    else if(utils::stringsEqual(type, "QuasiSpherical")) {
        checkNonzero(potential, type);
        if(density == NULL)
            density = potential;
        double beta0 = kvmap.getDoubleAlt("beta", "beta0", 0);
        double r_a   = kvmap.getDoubleAlt("anisotropyRadius", "r_a", INFINITY) * converter.lengthUnit;
        return PtrDistributionFunction(new QuasiSphericalCOM(
            potential::DensityWrapper(*density), potential::PotentialWrapper(*potential), beta0, r_a));
    }
    else  if(utils::stringsEqual(type, "Plummer")) {
	    checkNonzero(potential, type);
	    potential::Interpolator pot_interp(*potential);
	    return PtrDistributionFunction(new PlummerDF(parsePlummerParam(kvmap, converter)));
    }
    else  if(utils::stringsEqual(type, "Isochrone")) {
	    checkNonzero(potential, type);
	    potential::Interpolator pot_interp(*potential);
	    return PtrDistributionFunction(new IsochroneDF(parseIsochroneParam(kvmap, converter)));
    }
    else if(utils::stringsEqual(type, "NewOxford")) {
	    checkNonzero(potential, type);
	    NewOxfordParam par = parseNewOxfordParam(kvmap, converter);
	    potential::Interpolator pot_interp(*potential);
	    if(mass>0) {
		    par.norm = 1.0;
		    par.norm = mass / NewOxford(par,pot_interp).totalMass();
	    }
	    return PtrDistributionFunction(new NewOxford(par, pot_interp));
    }
    else{
	    printf("Unknown type of distribution function: %s\n",type.c_str()); exit(1);
    }
}

}  // namespace df
