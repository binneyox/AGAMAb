/** \file    potential_factory.h
    \brief   Creation and input/output of Potential instances
    \author  Eugene Vasiliev
    \date    2010-2015

    This file provides several utility function to manage instances of BaseDensity and BasePotential: 
    creating a density or potential model from parameters provided in ConfigPotential, 
    creating a potential from a set of point masses or from an N-body snapshot file,
    loading potential coefficients from a text file, 
    writing expansion coefficients to a text file,
    converting between potential parameters in `potential::ConfigPotential` structure and 
    a text array of key=value pairs (`utils::KeyValueMap`).
*/

#pragma once
#include "particles_base.h"
#include "potential_base.h"
#include "math_spline.h"
#include "smart.h"
#include "units.h"
#include "utils_config.h"
#include <string>

namespace potential {

/** Create an instance of density according to the parameters contained in the key-value map.
    \param[in] params is the list of parameters ("density=..." is a required one);
    \param[in] converter is the unit converter for transforming the dimensional quantities 
    in parameters (such as mass and radii) into internal units; can be a trivial converter.
    \return    a new instance of PtrDensity on success.
    \throw     std::invalid_argument or std::runtime_error or other density-specific exception on failure
*/
EXP PtrDensity createDensity(
    const utils::KeyValueMap& params,
    const units::ExternalUnits& converter = units::ExternalUnits());

/** Create an instance of potential according to the parameters contained in the key-value map.
    \param[in] params is the list of parameters ("type=..." is a required one);
    \param[in] converter is the unit converter for transforming the dimensional quantities 
    in parameters (such as mass and radii) into internal units; can be a trivial converter;
    \return    a new instance of PtrPotential on success.
    \throw     std::invalid_argument or std::runtime_error or other potential-specific exception
    on failure (e.g., if some of the parameters are invalid or missing, or refer to a non-existent file).
*/
EXP PtrPotential createPotential(
    const utils::KeyValueMap& params,
    const units::ExternalUnits& converter = units::ExternalUnits());

/** Create an instance of potential expansion for the user-provided density model,
    with parameters contained in the key-value map.
    \param[in] params is the list of parameters
    ("type=..." should specify one of the potential expansions);
    \param[in] dens   is the density model which will serve as the source to the potential;
    \param[in] converter (optional) is the unit converter for transforming dimensional quantities
    in parameters (essentially the grid sizes) into internal units;
    \return    a new instance of PtrPotential on success.
    \throw     std::invalid_argument if the requested potential is not of an expansion type,
    or any potential-specific exception on failure (if some parameters are missing or invalid).
*/
EXP PtrPotential createPotential(
    const utils::KeyValueMap& params,
    const BaseDensity& dens,
    const units::ExternalUnits& converter = units::ExternalUnits());

/** Create an instance of potential expansion approximating the user-provided potential model,
    with parameters contained in the key-value map; this is useful if the original potential
    is expensive to compute.
    \param[in] params is the list of parameters ("type=..." specifies the type of expansion);
    \param[in] pot    is the potential model which will be approximated with the potential expansion;
    \param[in] converter (optional) is the unit converter for transforming dimensional quantities
    in parameters (essentially the grid sizes) into internal units;
    \return    a new instance of PtrPotential on success.
    \throw     std::invalid_argument if the requested potential is not of an expansion type,
    or any potential-specific exception on failure (if some parameters are missing or invalid).
*/
EXP PtrPotential createPotential(
    const utils::KeyValueMap& params,
    const BasePotential& pot,
    const units::ExternalUnits& converter = units::ExternalUnits());

/** Create an instance of composite potential according to the parameters contained in the 
    array of key-value maps for each component. Note that there is no 1:1 correspondence 
    between the parameters and the components of the resulting potential, since the parameters 
    may contain several density models that would be used for initializing a single potential 
    expansion object (like in the case of GalPot).
    \param[in] params is the array of parameter lists, one per component (since vectors
    of references cannot exist, these are actual KeyValueMap objects rather than references).
    \param[in] converter is the unit converter for transforming the dimensional quantities 
    in parameters (such as mass and radii) into internal units; can be a trivial converter.
    \return    a new instance of PtrPotential on success.
    \throw     std::invalid_argument or std::runtime_error or other potential-specific exception
    on failure (e.g., if some of the parameters are invalid or missing, or refer to a non-existent file).
*/
EXP PtrPotential createPotential(
    const std::vector<utils::KeyValueMap>& params,
    const units::ExternalUnits& converter = units::ExternalUnits());

/** Create an instance of potential according to the parameters contained in an INI file.
    \param[in] iniFileName is the name of an INI file that contains one or more sections 
    with potential parameters, named as [Potential], [Potential1], ...
    \param[in] converter is the unit converter for transforming the dimensional quantities 
    in parameters (such as mass and radii) into internal units; can be a trivial converter.
    \return    a new instance of PtrPotential on success (if there are several components
    in the INI file, the returned potential is composite).
    \throw     std::invalid_argument or std::runtime_error or other potential-specific exception
    on failure (e.g., if some of the parameters are invalid or missing, or refer to a non-existent file).
*/
EXP PtrPotential createPotential(
    const std::string& iniFileName,
    const units::ExternalUnits& converter = units::ExternalUnits());

/** Create an instance of potential expansion from the provided array of particles.
    \param[in] params  is the list of required parameters (e.g., the type of potential expansion,
    number of terms, prescribed symmetry, etc.).
    \param[in] particles  is the array of particle positions and masses.
    \param[in] converter  is the unit converter for transforming the dimensional parameters 
    (min/max radii of grid) into internal units; can be a trivial converter. 
    Coordinates and masses of particles are _not_ transformed: if they are loaded from an external 
    N-body snapshot file, the conversion is applied at that stage, and if they come from 
    other routines in the library, they are already in internal units.
    \return    a new instance of PtrPotential on success.
    \throw     std::invalid_argument or std::runtime_error or other potential-specific exception
    on failure (e.g., if some of the parameters are invalid or missing).
*/
EXP PtrPotential createPotential(
    const utils::KeyValueMap& params, 
    const particles::ParticleArray<coord::PosCyl>& particles,
    const units::ExternalUnits& converter = units::ExternalUnits());


/** Construct an interpolated spherical density profile from two arrays -- radii and
    enclosed mass M(<r).
    First a suitably scaled interpolator is constructed for M(r);
    if it is found to have a finite limiting value at r --> infinity, the asymptotic power-law
    behaviour of density at large radii will be correctly represented.
    Then the density at each point of the radial grid is computed from the derivative of
    this interpolator. The returned array may be used to construct a LogLogSpline interpolator
    or a DensitySphericalHarmonic object (obviously, with only one harmonic).
    \param[in]  gridr  is the grid in radius (must have positive values sorted in order of increase);
    typically the radial grid should be exponentially spaced with r[i+1]/r[i] ~ 1.2 - 2.
    \param[in]  gridm  is the array of enclosed mass at each radius (must be positive and monotonic);
    \return  an array of density values at the given radii.
    \throw   std::invalid_argument if the input arrays were incorrect
    (incompatible sizes, non-monotinic or negative values), or
    std::runtime_error if the interpolator failed to produce a positive-definite density.
*/
EXP std::vector<double> densityFromCumulativeMass(
    const std::vector<double>& gridr,
    const std::vector<double>& gridm);


/** Read a file with the cumulative mass profile and construct a density model from it.
    The text file should be a whitespace- or comma-separated table with at least two columns
    (the rest is ignored) -- radius and the enclosed mass within this radius,
    both must be in increasing order. Lines not starting with a number are ignored.
    The enclosed mass profile should not include the central black hole (if present),
    because it could not be represented in terms of a density profile anyway.
    \param[in]  fileName  is the input file name.
    \return  an interpolated density profile, represented by a LogLogSpline class.
    \throw  std::runtime_error if the file does not exist, or the mass profile is not monotonic.
*/
EXP math::LogLogSpline readMassProfile(const std::string& fileName);


/** Utility function providing a legacy interface compatible with the original GalPot (deprecated).
    It reads the parameters from a text file and converts them into the internal unit system,
    using the conversion factors in the provided `units::ExternalUnits` object,
    then constructs the potential using `createGalaxyPotential()` routine.
    Standard GalPot units are Kpc and Msun, so to simplify matters, one may instead use 
    the overloaded function `readGalaxyPotential()` that takes the instance of internal units
    as the second argument, and creates a converter from standard GalPot to these internal units.
    \param[in]  filename  is the name of parameter file;
    \param[in]  converter provides the conversion from GalPot to internal units;
    \return     a new instance of PtrPotential;
    \throw      a std::runtime_error exception if file is not readable or does not contain valid parameters.
*/
EXP PtrPotential readGalaxyPotential(
    const std::string& filename,
    const units::ExternalUnits& converter);

/** Utility function providing a legacy interface compatible with the original GalPot (deprecated).
    It reads the parameters from a text file and converts them into the internal unit system, 
    then constructs the potential using `createGalaxyPotential()` routine.
    This function creates the instance of unit converter from standard GalPot units (Kpc and Msun)
    into the provided internal units, and calls the overloaded function `readGalaxyPotential()` 
    with this converter object as the second argument. 
    \param[in]  filename is the name of parameter file;
    \param[in]  unit     is the specification of internal unit system;
    \return     a new instance of PtrPotential;
    \throw      a std::runtime_error exception if file is not readable or does not contain valid parameters.
*/
inline PtrPotential readGalaxyPotential(
    const std::string& filename,
    const units::InternalUnits& unit) 
{   // create a temporary converter; velocity unit is not used
    return readGalaxyPotential(filename, units::ExternalUnits(unit, units::Kpc, units::kms, units::Msun));
}


/** Create a density expansion from coefficients stored in a text file.
    The file must contain coefficients for DensitySphericalHarmonic or DensityAzimuthalHarmonic;
    the density type is determined automatically from the first line of the file.
    \param[in] coefFileName specifies the file to read;
    \param[in] converter is the unit converter for transforming the density coefficients;
    from dimensional into internal units; can be a trivial converter;
    \return    a new instance of PtrDensity on success;
    \throw     std::invalid_argument or std::runtime_error or other density-specific exception
    on failure (e.g., if the file does not exist, or does not contain valid coefficients).
*/
EXP PtrDensity readDensity(
    const std::string& coefFileName,
    const units::ExternalUnits& converter = units::ExternalUnits());

/** Create a potential expansion from coefficients stored in a text file.
    The file must contain coefficients for BasisSetExp, SplineExp, CylSpline, or Multipole;
    the potential type is determined automatically from the first line of the file.
    \param[in] coefFileName specifies the file to read;
    \param[in] converter is the unit converter for transforming the potential coefficients;
    from dimensional into internal units; can be a trivial converter;
    \return    a new instance of PtrPotential on success;
    \throw     std::invalid_argument or std::runtime_error or other potential-specific exception
    on failure (e.g., if the file does not exist, or does not contain valid coefficients).
*/
EXP PtrPotential readPotential(
    const std::string& coefFileName,
    const units::ExternalUnits& converter = units::ExternalUnits());


/** Write density or potential expansion coefficients to a text file.
    The potential must be one of the following expansion classes: 
    `BasisSetExp`, `SplineExp`, `CylSpline`, `Multipole`,
    or the density may be `DensitySphericalHarmonic` or `DensityAzimuthalHarmonic`.
    The coefficients stored in a file may be later loaded by `readPotential()` or
    `readDensity()` routines.
    If the potential or density is composite, each component is saved into a separate file
    with suffixes "_0", "_1", etc. attached to the name, and the list of these files is
    stored in the main file.
    \param[in] fileName is the output file;
    \param[in] density is the reference to density or potential object;
    \param[in] converter is the unit converter for transforming the density or potential
    coefficients from internal into dimensional units; can be a trivial converter;
    \return    success or failure (the latter may also mean that export is 
    not available for this type of potential/density).
*/
EXP bool writeDensity(
    const std::string& fileName,
    const BaseDensity& density,
    const units::ExternalUnits& converter = units::ExternalUnits());

/// alias to writeDensity
inline bool writePotential(
    const std::string& fileName,
    const BasePotential& potential,
    const units::ExternalUnits& converter = units::ExternalUnits()) {
    return writeDensity(fileName, potential, converter); }


/// return file extension for writing the coefficients of a given potential type,
/// or empty string if it is neither one of the expansion types nor a composite potential
EXP const char* getCoefFileExtension(const std::string& potName);

/// return file extension for writing the coefficients of a given potential object,
/// or empty string if the potential type is not one of the expansion types
inline const char* getCoefFileExtension(const BasePotential& p) {
    return getCoefFileExtension(p.name()); }

/** return the symmetry type encoded in the string.
    Spherical, Axisymmetric, Triaxial and None are recognized by the first letter,
    whereas other types must be given by their numerical code.
    If the string is empty, the default value ST_TRIAXIAL is returned.
*/
EXP coord::SymmetryType getSymmetryTypeByName(const std::string& name);

/** return the name of symmetry encoded in SymmetryType.
    Spherical, Axisymmetric, Triaxial and None are returned as symbolic names,
    other types as their numerical code.
*/
EXP std::string getSymmetryNameByType(coord::SymmetryType type);

}  // namespace potential
