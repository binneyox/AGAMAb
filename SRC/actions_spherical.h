/** \file    actions_spherical.h
    \brief   Action-angle finders for a generic spherical potential
    \author  Eugene Vasiliev
    \date    2015-2016
*/
#pragma once
#include "actions_base.h"
#include "potential_utils.h"
#include "math_spline.h"
#include "smart.h"

namespace actions {

/** Compute actions in a given spherical potential.
    \param[in]  potential is the arbitrary spherical potential;
    \param[in]  point     is the position/velocity point;
    \return     actions for the given point, or Jr=NAN if the energy is positive;
    \throw      std::invalid_argument exception if the potential is not spherical
    or some other error occurs.
*/
EXP Actions actionsSpherical(
			     const potential::BasePotential& potential,
			     const coord::PosVelCyl& point);
EXP Actions actionsSpherical(
			     const potential::BasePotential& potential,
			     const coord::PosMomCyl& point);

/** Compute actions, angles and frequencies in a spherical potential.
    \param[in]  potential is the arbitrary spherical potential;
    \param[in]  point     is the position/velocity point;
    \param[out] freq      if not NULL, output the frequencies in this variable;
    \return     actions and angles for the given point, or Jr=NAN if the energy is positive;
    \throw      std::invalid_argument exception if the potential is not spherical
    or some other error occurs.
*/ 
EXP ActionAngles actionAnglesSpherical(
    const potential::BasePotential& pot, const coord::PosVelCyl& point, Frequencies* freqout);
    
EXP ActionAngles actionAnglesSpherical(
    const potential::BasePotential& potential,
    const coord::PosMomCyl& point,
    Frequencies* freq=NULL);


/** Compute the total energy for an orbit in a spherical potential from the given values of actions.
    \param[in]  potential  is the arbitrary spherical potential;
    \param[in]  acts       are the actions;
    \return     the value of Hamiltonian (total energy) corresponding to the given actions;
    \throw      std::invalid_argument exception if the potential is not spherical
    or Jr/Jz actions are negative.
*/
EXP double computeHamiltonianSpherical(const potential::BasePotential& potential, const Actions& acts);


/** Compute position/velocity from actions/angles in an arbitrary spherical potential.
    \param[in]  potential   is the instance of a spherical potential;
    \param[in]  actAng  is the action/angle point
    \param[out] freq    if not NULL, store the frequencies for these actions.
    \return     position and velocity point
*/
EXP coord::PosMomCyl mapSpherical(
    const potential::BasePotential &potential,
    const ActionAngles &actAng, Frequencies* freq=NULL);


/** Class for performing transformations between action/angle and coordinate/momentum for
    an arbitrary spherical potential, using 2d interpolation tables */
class EXP ActionFinderSpherical: public BaseActionFinder {
public:
    /// Initialize the internal interpolation tables; the potential itself is not used later on
    explicit ActionFinderSpherical(const potential::BasePotential& potential);

    virtual Actions actions(const coord::PosVelCyl& point) const;
    virtual ActionAngles actionAngles(const coord::PosVelCyl& point, Frequencies* freq=NULL) const;
    virtual coord::PosMomSphMod map(
        const ActionAngles& actAng,
        Frequencies* freq=NULL,
        DerivAct<coord::SphMod>* derivAct=NULL,
        DerivAng<coord::SphMod>* derivAng=NULL,
        coord::PosMomSphMod* derivParam=NULL) const;

    /** return the interpolated value of radial action as a function of energy and angular momentum;
        also return the frequencies in Omegar and Omegaz if these arguments are not NULL */
    double Jr(double E, double L, double *Omegar=NULL, double *Omegaz=NULL) const;

    /** return the energy corresponding to the given actions */
    double E(const Actions& act) const;
private:
    const double invPhi0;                 ///< 1/(value of potential at r=0)
    const potential::Interpolator2d pot;  ///< interpolator for potential and peri/apocenter radii
    const math::QuinticSpline2d intJr;    ///< interpolator for Jr(E,L)
    const math::QuinticSpline2d intE;     ///< interpolator for E(Jr,L)
};

typedef ActionFinderSpherical ToyMapSpherical;
EXP void mapHJr(const potential::BasePotential &pot,math::QuinticSpline2d& intJr,math::QuinticSpline2d& intE);
}  // namespace actions
