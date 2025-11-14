#!/usr/bin/python
"""
This example demonstrates the machinery for constructing multicomponent self-consistent models
specified by distribution functions (DFs) in terms of actions.
We create a four-component galaxy with disk, bulge and halo components defined by their DFs,
and a static density profile of gas disk.
Then we perform several iterations of recomputing the density profiles of components from their DFs
and recomputing the total potential.
Finally, we create N-body representations of all mass components: dark matter halo,
stars (bulge, thin and thick disks and stellar halo combined), and gas disk.
A modification of this script that creates a self-consistent three-component model
(disk, bulge and halo) is given in example_self_consistent_model3.py
This example is the Python counterpart of tests/example_self_consistent_model.cpp
"""
import agama as ag
import numpy, sys, os
intUnits=ag.IntUnits(ag.Kpc, ag.Myr)
extUnits=ag.ExtUnits(intUnits, 1.*ag.Kpc, 1.*ag.kms, 1.*ag.Msun)
solarRadius=0
# write out the rotation curve (separately for each component, and the total one)
def writeRotationCurve(filename, potentials):
    radii = numpy.logspace(-2., 2., 81)*intUnits.from_Kpc
    vcomp=numpy.column_stack([ag.Vcirc(potential,radii) for potential in potentials])
    vtot=numpy.sum(vcomp**2,axis=1)**0.5
    numpy.savetxt(filename, numpy.column_stack((radii, vtot, vcomp)), fmt="%.6g", delimiter="\t", \
        header="radius[Kpc]\tv_circ,total[km/s]\tdisk\tbulge\thalo")

# print surface density profiles to a file
def writeSurfaceDensityProfile(filename, model):
    print("Writing surface density profile")
    radii = numpy.hstack(([1./8, 1./4], numpy.linspace(0.5, 16, 32), numpy.linspace(18, 30, 7)))*intUnits.from_Kpc
    
    nc=model.distrFunc().numValues()
    surfDen=numpy.zeros((len(radii),nc))
    for ir in range(len(radii)):
        (surfDen1,rmsHeight,rmsVel)=ag.computeProjectedMoments(model, radii[ir],seperate=True)
        if(nc>1):
            for i in range(nc):
                surfDen[ir,i]=surfDen1[i]
        else:
            surfDen[ir,0]=surfDen1

    numpy.savetxt(filename, numpy.column_stack((radii*intUnits.to_Kpc, surfDen*intUnits.to_Msun_per_pc2)), fmt="%.6g", delimiter="\t", \
        header="Radius[Kpc]\tThinDisk\tThickDisk\tStellarHalo:SurfaceDensity[Msun/pc^2]")

# print vertical density profile for several sub-components of the stellar DF
def writeVerticalDensityProfile(filename,pot,af, DFcomponents):
    print("Writing vertical density profile")
    heights = numpy.hstack((numpy.linspace(0, 1.5, 13), numpy.linspace(2, 8, 13)))*intUnits.from_Kpc
    print(heights)
    R = solarRadius * intUnits.from_Kpc
    nh = len(heights)
    nc=DFcomponents.numValues()
    M=numpy.zeros((nh,nc+1))
    for ih in range(nh):
        M[ih,0]=heights[ih]*intUnits.to_Kpc
    for i in range(nh*nc):
        ih = int(numpy.floor(i/nc))
        ic = i%nc
        M[ih,ic]=ag.computeMoments(ag.GalaxyModel(pot, af, DFcomponents.component(ic)),\
            ag.PosCyl(R,heights[ih],0),Dens=True)*intUnits.to_Msun_per_pc3
    numpy.savetxt(filename, M, fmt="%.6g", delimiter="\t", \
        header="z[Kpc]\tThinDisk\tThickDisk\tStellarHalo:Density[Msun/pc^3]")

# print velocity dispersion profiles in the equatorial plane as functions of radius to a file
def writeVelocityDispersionProfile(filename, model):
    print("Writing velocity dispersion profile")
    radii = numpy.hstack(([1./8, 1./4], numpy.linspace(0.5, 16, 32), numpy.linspace(18, 30, 7)))*intUnits.from_Kpc
    nr=len(radii)
    nc=model.distrFunc().numValues()
    vel=numpy.zeros((nr,nc,3))
    vel2=numpy.zeros((nr,nc,3))
    for ir in range(nr):
        (dens,vel[ir,:,1],veloc2)=ag.computeMoments(model,ag.PosCyl(radii[ir],0,0),True,False,True,True,True)
        for i in range(nc):
            vel2[ir,i,0]=veloc2[i].vR2
            vel2[ir,i,1]=veloc2[i].vphi2
            vel2[ir,i,2]=veloc2[i].vz2
    vel2[:,:,1] -= vel[:,:,1]**2
    A=(vel2[:,:,0:3]**0.5, vel[:,:,1:2])*intUnits.to_kms
    numpy.savetxt(filename, numpy.column_stack((radii, numpy.dstack(A).reshape(len(radii),-1))),
        fmt="%.6g", delimiter="\t", header="Radius[Kpc]\t"
        "ThinDisk:sigma_r\tsigma_phi\tsigma_z\tv_phi\t"
        "ThickDisk:sigma_r\tsigma_phi\tsigma_z\tv_phi\t"
        "StellarHalo:sigma_r\tsigma_phi\tsigma_z\tv_phi[km/s]")

# print velocity distributions at the given point to a file
def writeVelocityDistributions(filename, model):
    point = ag.PosCyl(solarRadius, 0.1* intUnits.from_Kpc,0)
    print("Writing velocity distributions at (x=%g, z=%g)" % (point.R, point.z))
    # create grids in velocity space for computing the spline representation of VDF
    v_max = 360.0*intUnits.from_kms    # km/s
    gridv = numpy.linspace(-v_max, v_max, 75) # use the same grid for all dimensions
    interp=ag.BsplineInterpolator1d3(gridv)
    (dens,amplvR,amplvphi,amplvz)=ag.computeVelocityDistributionO3(model, point,False,\
        gridv, gridv, gridv)
    # compute the distributions (represented as cubic splines)
    # output f(v) at a different grid of velocity values
    gridv = numpy.linspace(-v_max, v_max, 201)
    A=numpy.zeros((len(gridv),4))
    for i in range(len(gridv)):
        A[i,0]=gridv[i]*intUnits.to_kms
        A[i,1]=interp.interpolate(gridv[i],amplvR,0)/ intUnits.to_kms
        A[i,2]=interp.interpolate(gridv[i],amplvz,0)/ intUnits.to_kms
        A[i,3]=interp.interpolate(gridv[i],amplvphi,0)/ intUnits.to_kms
    numpy.savetxt(filename, A,
        fmt="%.6g", delimiter="\t", header="V\tf(V_R)\tf(V_z)\tf(V_phi) [1/(km/s)]")
# display some information after each iteration
def printoutInfo(model, iteration):
    compDisk = model.components[0].getDensity()
    compBulge= model.components[1].getDensity()
    compHalo = model.components[2].getDensity()
    pt0 = ag.PosCar(solarRadius, 0, 0)
    pt1 = ag.PosCar(solarRadius, 0, 1)
    print("Disk total mass=%g Msun, rho(Rsolar,z=0)=%g, rho(Rsolar,z=1kpc)=%g Msun/pc^3" % \
        (compDisk.totalMass()* intUnits.to_Msun, compDisk.density(pt0)* intUnits.to_Msun_per_pc3, compDisk.density(pt1)* intUnits.to_Msun_per_pc3))  # per pc^3, not kpc^3
    print("Halo total mass=%g Msun, rho(Rsolar,z=0)=%g, rho(Rsolar,z=1kpc)=%g Msun/pc^3" % \
        (compHalo.totalMass()* intUnits.to_Msun, compHalo.density(pt0)* intUnits.to_Msun_per_pc3, compHalo.density(pt1)* intUnits.to_Msun_per_pc3))
    print("Potential at origin=-(%g km/s)^2, total mass=%g Msun" % \
        ((-model.totalPotential.value(ag.PosCar(0,0,0)))**0.5* intUnits.to_kms, model.totalPotential.totalMass()* intUnits.to_Msun))
    potentials=[]
    potentials.append(ag.createCylSpline(compDisk,model.mmaxAngularCyl,
            model.sizeRadialCyl,   model.RminCyl, model.RmaxCyl,
            model.sizeVerticalCyl, model.zminCyl, model.zmaxCyl, True))
    potentials.append(ag.createMultipole(compBulge,6,0,100))
    potentials.append(ag.createMultipole(compHalo,6,0,100))
    writeRotationCurve("rotcurve_"+iteration, potentials)


if __name__ == "__main__":
    
    ini=ag.ConfigFile("SCM.ini")
    iniPotenThinDisk = ini.findSection("Potential thin disk")
    iniPotenThickDisk= ini.findSection("Potential thick disk")
    iniPotenGasDisk  = ini.findSection("Potential gas disk")
    iniPotenBulge    = ini.findSection("Potential bulge")
    iniPotenDarkHalo = ini.findSection("Potential dark halo")
    iniDFThinDisk    = ini.findSection("DF thin disk")
    iniDFThickDisk   = ini.findSection("DF thick disk")
    iniDFStellarHalo = ini.findSection("DF stellar halo")
    iniDFBulge       = ini.findSection("DF bulge")
    iniDFDarkHalo    = ini.findSection("DF dark halo")
    iniSCMDisk       = ini.findSection("SelfConsistentModel disk")
    iniSCMBulge      = ini.findSection("SelfConsistentModel bulge")
    iniSCMHalo       = ini.findSection("SelfConsistentModel halo")
    iniSCM           = ini.findSection("SelfConsistentModel")
    solarRadius = ini.findSection("Data").getDouble("SolarRadius", solarRadius)
    model=ag.SelfConsistentModel()
    model.rminSph         = iniSCM.getDouble("rminSph") * extUnits.lengthUnit
    model.rmaxSph         = iniSCM.getDouble("rmaxSph") * extUnits.lengthUnit
    model.sizeRadialSph   = iniSCM.getInt("sizeRadialSph")
    model.lmaxAngularSph  = iniSCM.getInt("lmaxAngularSph")
    model.RminCyl         = iniSCM.getDouble("RminCyl") * extUnits.lengthUnit
    model.RmaxCyl         = iniSCM.getDouble("RmaxCyl") * extUnits.lengthUnit
    model.zminCyl         = iniSCM.getDouble("zminCyl") * extUnits.lengthUnit
    model.zmaxCyl         = iniSCM.getDouble("zmaxCyl") * extUnits.lengthUnit
    model.sizeRadialCyl   = iniSCM.getInt("sizeRadialCyl")
    model.sizeVerticalCyl = iniSCM.getInt("sizeVerticalCyl")
    model.useActionInterpolation = iniSCM.getBool("useActionInterpolation")

#initialize density profiles of various components
    densityStellarDisk=[]
    densityBulge    = ag.createDensity(iniPotenBulge,    extUnits)
    densityDarkHalo =ag.createDensity(iniPotenDarkHalo, extUnits)
    densityStellarDisk.append(ag.createDensity(iniPotenThinDisk, extUnits))
    densityStellarDisk.append(ag.createDensity(iniPotenThickDisk, extUnits))
    densityGasDisk  =ag.createDensity(iniPotenGasDisk,  extUnits)

#add components to SCM - at first, all of them are static density profiles
    A= model.components
    A.append(ag.ComponentStatic(ag.CompositeDensity(densityStellarDisk),True))
    A.append(ag.ComponentStatic(densityBulge,False))
    A.append(ag.ComponentStatic(densityDarkHalo,False))
    A.append(ag.ComponentStatic(densityGasDisk,True))
    model.components=A
#initialize total potential of the model (first guess)
    ag.updateTotalPotential(model)
    printoutInfo(model, "init")
    Text="\033[1;33m**** STARTING MODELLING ****\033[0m\nInitial masses of density components: "
    Text+= "Mdisk="+str(model.components[0].getDensity().totalMass() * intUnits.to_Msun)+" Msun, "
    Text+= "Mbulge="+str(densityBulge.totalMass() * intUnits.to_Msun)+" Msun, "
    Text+= "Mhalo="+str(densityDarkHalo.totalMass() * intUnits.to_Msun)+" Msun, "
    Text+= "Mgas="+str(densityGasDisk.totalMass() * intUnits.to_Msun)+" Msun, "
    print(Text)
#create the dark halo DF
    dfHalo = ag.createDistributionFunction(iniDFDarkHalo, model.totalPotential,extUnits)
#same for the bulge
    dfBulge = ag.createDistributionFunction(iniDFBulge, model.totalPotential,extUnits)
#same for the stellar components (thin/thick disks and stellar halo)
    dfStellarArray=[]
    dfStellarArray.append(ag.createDistributionFunction(iniDFThinDisk,model.totalPotential,extUnits))
    dfStellarArray.append(ag.createDistributionFunction(iniDFThickDisk, model.totalPotential,extUnits))
    dfStellarArray.append(ag.createDistributionFunction(iniDFStellarHalo,model.totalPotential,extUnits))
#composite DF of all stellar components except the bulge
    dfStellar=ag.CompositeDF(dfStellarArray)
#replace the static disk density component of SCM with a DF-based disk component
    initDens=ag.createDensity(ag.KeyValueMap("type=Isochrone"),extUnits)
    model.components[0] = ag.ComponentWithDisklikeDF(dfStellar, initDens,iniSCMDisk.getInt("mmaxAngularCyl"), 
                                                     iniSCMDisk.getInt("sizeRadialCyl"), iniSCMDisk.getDouble("RminCyl") * extUnits.lengthUnit,iniSCMDisk.getDouble("RmaxCyl") * extUnits.lengthUnit,iniSCMDisk.getInt("sizeVerticalCyl"),iniSCMDisk.getDouble("zminCyl") * extUnits.lengthUnit,iniSCMDisk.getDouble("zmaxCyl") * extUnits.lengthUnit)
    #model.components[0] = ag.ComponentWithDisklikeDF(dfBulge, initDens,iniSCMDisk.getInt("mmaxAngularCyl"), iniSCMDisk.getInt("sizeRadialCyl"), iniSCMDisk.getDouble("RminCyl") * extUnits.lengthUnit,iniSCMDisk.getDouble("RmaxCyl") * extUnits.lengthUnit,iniSCMDisk.getInt("sizeVerticalCyl"),iniSCMDisk.getDouble("zminCyl") * extUnits.lengthUnit,iniSCMDisk.getDouble("zmaxCyl") * extUnits.lengthUnit)
#same for the bulge
    model.components[1] = ag.ComponentWithSpheroidalDF(dfBulge,initDens, iniSCMBulge.getInt("lmaxAngularSph"),iniSCMBulge.getInt("mmaxAngularSph"),
                                                       iniSCMBulge.getInt("sizeRadialSph"),iniSCMBulge.getDouble("rminSph") * extUnits.lengthUnit,iniSCMBulge.getDouble("rmaxSph") * extUnits.lengthUnit)
#same for the halo
    model.components[2] = ag.ComponentWithSpheroidalDF(dfHalo,initDens,iniSCMHalo.getInt("lmaxAngularSph"),iniSCMHalo.getInt("mmaxAngularSph"),iniSCMHalo.getInt("sizeRadialSph"),iniSCMHalo.getDouble("rminSph") * extUnits.lengthUnit,iniSCMHalo.getDouble("rmaxSph") * extUnits.lengthUnit)

    Text="Masses of DF components:"
    Text=" Mdisk=" +str((dfStellar.totalMass() * intUnits.to_Msun))
    Text+= " Msun (Mthin="+str(dfStellarArray[0].totalMass() * intUnits.to_Msun)
    Text +=", Mthick="+str(dfStellarArray[1].totalMass() * intUnits.to_Msun)
    Text+=", Mstel.halo=" +str(dfStellarArray[2].totalMass() * intUnits.to_Msun)
    Text+="); Mbulge="+str(dfBulge.totalMass() * intUnits.to_Msun)+" Msun"
    Text+="; Mdarkhalo="+str(dfHalo.totalMass() * intUnits.to_Msun)+" Msun\n"
#do a few more iterations to obtain the self-consistent density profile for both disks
    for iteration in range(5):
        ag.doIteration(model)

#output various profiles (only for stellar components)
    print("\033[1;33mComputing density profiles and velocity distribution\033[0m")
    modelStars=ag.GalaxyModel(model.totalPotential, model.actionFinder, dfStellar)
    writeSurfaceDensityProfile("model_stars_final.surfdens", modelStars)
    writeVerticalDensityProfile("model_stars_final.vertical",model.totalPotential, model.actionFinder,dfStellar)
    #writeVelocityDispersionProfile("model_stars_final.veldisp",  modelStars)
    writeVelocityDistributions("model_stars_final.veldist",  modelStars)

#export model to an N-body snapshot
    print("\033[1;33mCreating an N-body representation of the model\033[0m")
    format = "text"

#first create a representation of density profiles without velocities
#(just for demonstration), by drawing samples from the density distribution
    print("Writing N-body sampled density profile for the dark matter halo")
    ag.writeSnapshot("dens_dm_final", ag.sampleDensity(
        model.components[2].getDensity(), 800000), format, extUnits)
    print("Writing N-body sampled density profile for the stellar bulge, disk and halo\n")
    densityStars=[]
    densityStars.append(model.components[0].getDensity())
    densityStars.append(model.components[1].getDensity())
    ag.writeSnapshot("dens_stars_final", ag.sampleDensity(
        ag.CompositeDensity(densityStars), 200000), format, extUnits)

#now create genuinely self-consistent models of both components,
#by drawing positions and velocities from the DF in the given (self-consistent) potential
    print("Writing a complete DF-based N-body model for the dark matter halo\n")
    ag.writeSnapshot("model_dm_final", ag.samplePosVel(
        ag.GalaxyModel(model.totalPotential, model.actionFinder, dfHalo), 800000),
        format, extUnits)
    print("Writing a complete DF-based N-body model for the stellar bulge, disk and halo\n")
    dfStellarArray.append(dfBulge)
    dfStellar=ag.CompositeDF(dfStellarArray)
    ag.writeSnapshot("model_stars_final", ag.samplePosVel(
        ag.GalaxyModel(model.totalPotential, model.actionFinder, dfStellar), 200000),
        format, extUnits)
#we didn't use an action-based DF for the gas disk, leaving it as a static component;
# to create an N-body representation, we sample the density profile and assign velocities
#from the axisymmetric Jeans equation with equal velocity dispersions in R,z,phi
    print("Writing an N-body model for the gas disk\n")
    ag.writeSnapshot("model_gas_final",ag.assignVelocity(
        ag.sampleDensity(model.components[3].getDensity(), 24000),
        model.components[3].getDensity(),model.totalPotential, 0., 1.),
        format, extUnits)