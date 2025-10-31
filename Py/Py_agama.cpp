#include <variant>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "actions_base.h"
#include "actions_newtorus.h"
#include "actions_staeckel.h"
#include "potential_factory.h"
#include "potential_composite.h"
#include "potential_cylspline.h"
#include "potential_multipole.h"
#include "potential_analytic.h"
#include "particles_io.h"
#include "math_spline.h"
#include "orbit.h"
#include "coord.h"
#include "obs_base.h"
#include "units.h"
#include "utils.h"
#include "utils_config.h"
#include "df_base.h"
#include "df_factory.h"
#include "df_spherical.h"
#include "df_halo.h"
#include "galaxymodel_base.h"
#include "galaxymodel_selfconsistent.h"
#include "galaxymodel_velocitysampler.h"
#include <stdio.h>
namespace py = pybind11;
using namespace pybind11::literals;
namespace actions{
EXP class Jzfind{
public:
    std::vector<double> params;
    Jzfind(const potential::BasePotential &_pot):params(mapJcrit(_pot)){};
    double value(const double Jr){return math::evalPoly(params,Jr);}
};
}
PYBIND11_MODULE(Py_agama, m) {
    using Bspl13=math::BsplineInterpolator1d<3>;
    py::class_<Bspl13>(m,"BsplineInterpolator1d3")
        .def(py::init<const std::vector<double>&>())
        .def("integrate",&Bspl13::integrate)
        .def("interpolate",&Bspl13::interpolate)
        .def("antideriv",&Bspl13::antideriv)
        .def("xmax",&Bspl13::xmax)
        .def("xmin",&Bspl13::xmin)
        .def("xvalues",&Bspl13::xvalues)
        .def("deriv",&Bspl13::deriv)
        .def("numValues",&Bspl13::numValues)
        .def("numVars",&Bspl13::numVars);
    py::class_<utils::KeyValueMap>(m,"KeyValueMap")
        .def(py::init<>())
        .def(py::init<const std::string &,const std::string &>(),"params"_a,"whitespace"_a=", ")
        .def("add",&utils::KeyValueMap::add)
        .def("contains",&utils::KeyValueMap::contains)
        .def("dump",&utils::KeyValueMap::dump)
        .def("dumpSingleLine",&utils::KeyValueMap::dumpSingleLine)
        .def("getBool",&utils::KeyValueMap::getBool,"key"_a,"defaultValue"_a=false)
        .def("getDouble",&utils::KeyValueMap::getDouble,"key"_a,"defaultValue"_a=0.0)
        .def("getDoubleAlt",&utils::KeyValueMap::getDoubleAlt,"key1"_a,"key2"_a,"defaultValue"_a=0.0)
        .def("getDoubleVector",&utils::KeyValueMap::getDoubleVector)
        .def("getInt",&utils::KeyValueMap::getInt,"key"_a,"defaultValue"_a=0)
        .def("getIntAlt",&utils::KeyValueMap::getIntAlt,"key1"_a,"key2"_a,"defaultValue"_a=0)
        .def("getString",&utils::KeyValueMap::getString,"key"_a,"defaultValue"_a="")
        .def("getStringAlt",&utils::KeyValueMap::getStringAlt,"key1"_a,"key2"_a,"defaultValye"_a="")
        .def("isModified",&utils::KeyValueMap::isModified)
        .def("keys",&utils::KeyValueMap::keys)
        .def("set",[](utils::KeyValueMap &self,const std::string &key,const bool value){return self.set(key,value);})
        .def("set",[](utils::KeyValueMap &self,const std::string &key,const unsigned int value){return self.set(key,value);})
        .def("set",[](utils::KeyValueMap &self,const std::string &key,const int value){return self.set(key,value);})
        .def("set",[](utils::KeyValueMap &self,const std::string &key,const char* value){return self.set(key,value);})
        .def("set",[](utils::KeyValueMap &self,const std::string &key,const std::string &value){return self.set(key,value);})
        .def("set",[](utils::KeyValueMap &self,const std::string &key,const double value){return self.set(key,value);})
        .def("unset",&utils::KeyValueMap::unset);
    py::class_<utils::ConfigFile>(m,"ConfigFile")
        .def(py::init<const std::string &,bool>(),"fileName"_a,"mustExist"_a=true)
        .def("findSection",&utils::ConfigFile::findSection)
        .def("listSections",[](utils::ConfigFile &self){return self.listSections();});
    py::class_<units::InternalUnits>(m,"IntUnits")
        .def(py::init<double, double>())
        .def_readonly("from_Gev_per_cm3", &units::InternalUnits::from_Gev_per_cm3)
        .def_readonly("from_Gyr", &units::InternalUnits::from_Gyr)
        .def_readonly("from_kms", &units::InternalUnits::from_kms)
        .def_readonly("from_Kpc_kms", &units::InternalUnits::from_Kpc_kms)
        .def_readonly("from_ly", &units::InternalUnits::from_ly)
        .def_readonly("from_mas_per_yr", &units::InternalUnits::from_mas_per_yr)
        .def_readonly("from_Msun", &units::InternalUnits::from_Msun)
        .def_readonly("from_Mpc", &units::InternalUnits::from_Mpc)
        .def_readonly("from_Msun_per_Kpc3", &units::InternalUnits::from_Msun_per_Kpc3)
        .def_readonly("from_Msun_per_Kpc2", &units::InternalUnits::from_Msun_per_Kpc2)
        .def_readonly("from_Msun_per_pc3", &units::InternalUnits::from_Msun_per_pc3)
        .def_readonly("from_Msun_per_pc2", &units::InternalUnits::from_Msun_per_pc2)
	.def_readonly("from_pc", &units::InternalUnits::from_pc)
	.def_readonly("from_Kpc", &units::InternalUnits::from_Kpc)
        .def_readonly("from_yr", &units::InternalUnits::from_yr)
        .def_readonly("from_Myr", &units::InternalUnits::from_Myr)
        .def_readonly("to_Gev_per_cm3", &units::InternalUnits::to_Gev_per_cm3)
        .def_readonly("to_Gyr", &units::InternalUnits::to_Gyr)
        .def_readonly("to_kms", &units::InternalUnits::to_kms)
        .def_readonly("to_Kpc_kms", &units::InternalUnits::to_Kpc_kms)
        .def_readonly("to_ly", &units::InternalUnits::to_ly)
        .def_readonly("to_mas_per_yr", &units::InternalUnits::to_mas_per_yr)
        .def_readonly("to_Msun", &units::InternalUnits::to_Msun)
        .def_readonly("to_Mpc", &units::InternalUnits::to_Mpc)
        .def_readonly("to_Msun_per_Kpc3", &units::InternalUnits::to_Msun_per_Kpc3)
        .def_readonly("to_Msun_per_Kpc2", &units::InternalUnits::to_Msun_per_Kpc2)
        .def_readonly("to_Msun_per_pc3", &units::InternalUnits::to_Msun_per_pc3)
        .def_readonly("to_Msun_per_pc2", &units::InternalUnits::to_Msun_per_pc2)
        .def_readonly("to_pc", &units::InternalUnits::to_pc)
        .def_readonly("to_Kpc", &units::InternalUnits::to_Kpc)
        .def_readonly("to_yr", &units::InternalUnits::to_yr)
        .def_readonly("to_Myr", &units::InternalUnits::to_Myr);
    m.attr("Kpc")=units::Kpc;
    m.attr("Myr")=units::Myr;
    m.attr("Msun")=units::Msun;
    m.attr("kms")=units::kms;
    m.attr("galactic_kms") = units::galactic_kms;
    m.attr("galactic_Myr") = units::galactic_Myr;
    py::class_<units::ExternalUnits>(m,"ExtUnits")
        .def(py::init<const units::InternalUnits &,double, double,double>())
        .def_readonly("lengthUnit", &units::ExternalUnits::lengthUnit)
        .def_readonly("massUnit", &units::ExternalUnits::massUnit)
        .def_readonly("velocityUnit", &units::ExternalUnits::velocityUnit);
    py::class_<coord::PosCar>(m,"PosCar")
        .def(py::init<double, double, double >(),"x"_a,"y"_a,"z"_a)
        .def_readwrite("x", &coord::PosCar::x)
        .def_readwrite("y", &coord::PosCar::y)
        .def_readwrite("z", &coord::PosCar::z)
        .def("__repr__",[](const coord::PosCar &x) {
            char buffer[100];
            std::snprintf(buffer,100,"PosCar: [%.8g,%.8g,%.8g]",x.x,x.y,x.z);
            return std::string(buffer);
        })
        .doc()="Class PosCar is used to store position in cartesian coordinates.\nTakes in 3 aguments: x,y,z.";
    py::class_<coord::VelCar>(m,"VelCar")
        .def(py::init<double, double, double >(),"vx"_a,"vy"_a,"vz"_a)
        .def_readwrite("vx", &coord::VelCar::vx)
        .def_readwrite("vy", &coord::VelCar::vy)
        .def_readwrite("vz", &coord::VelCar::vz)
        .def("__repr__",[](const coord::VelCar &v) {
            char buffer[100];
            std::snprintf(buffer,100,"VelCar: (%.8g,%.8g,%.8g)",v.vx,v.vy,v.vz);
            return std::string(buffer);
        })
        .doc()="Class VelCar is used to store velocity in cartesian coordinates.\nTakes in 3 aguments: vx,vy,vz as floats.";
    py::class_<coord::PosVelCar>(m,"PosVelCar")
        .def(py::init<double, double, double, double, double, double>(),"x"_a,"y"_a,"z"_a,"vx"_a,"vy"_a,"vz"_a)
        .def_readwrite("x", &coord::PosVelCar::x)
        .def_readwrite("y", &coord::PosVelCar::y)
        .def_readwrite("z", &coord::PosVelCar::z)
        .def_readwrite("vx", &coord::PosVelCar::vx)
        .def_readwrite("vy", &coord::PosVelCar::vy)
        .def_readwrite("vz", &coord::PosVelCar::vz)
        .def("__repr__",[](const coord::PosVelCar &xv) {
            char buffer[100];
            std::snprintf(buffer,100,"PosVelCar: ( %.8g, %.8g, %.8g, %.8g, %.8g, %.8g)",
                xv.x,xv.y,xv.z,xv.vx,xv.vy,xv.vz);
            return std::string(buffer);
        })
        .doc()="Class PosVelCar is used to store position and velocity in cartesian coordinates.\n"
        "Takes in 6 aguments: x, y, z, vx, vy, vz as floats.";
    py::class_<coord::PosMomCar>(m,"PosMomCar")
        .def(py::init<double, double, double, double, double, double>(),"x"_a,"y"_a,"z"_a,"px"_a,"py"_a,"pz"_a)
        .def_readwrite("x", &coord::PosMomCar::x)
        .def_readwrite("y", &coord::PosMomCar::y)
        .def_readwrite("z", &coord::PosMomCar::z)
        .def_readwrite("px", &coord::PosMomCar::px)
        .def_readwrite("py", &coord::PosMomCar::py)
        .def_readwrite("pz", &coord::PosMomCar::pz)
        .def("__repr__",[](const coord::PosMomCar &xp) {
            char buffer[100];
            std::snprintf(buffer,100,"PosMomCar: ( %.8g, %.8g, %.8g, %.8g, %.8g, %.8g)",
                xp.x,xp.y,xp.z,xp.px,xp.py,xp.pz);
            return std::string(buffer);
        })
        .doc()="Class PosMomCar is used to store position and momentum in cartesian coordinates.\n"
        "Takes in 6 aguments: x, y, z, px, py, pz as floats.";
    py::class_<coord::PosCyl>(m,"PosCyl")
        .def(py::init<double, double, double >(),"R"_a,"z"_a,"phi"_a)
        .def_readwrite("R", &coord::PosCyl::R)
        .def_readwrite("z", &coord::PosCyl::z)
        .def_readwrite("phi", &coord::PosCyl::phi)
        .def("__repr__",[](const coord::PosCyl &R) {
            char buffer[100];
            std::snprintf(buffer,100,"PosCyl: ( %.8g, %.8g, %.8g)",R.R,R.z,R.phi);
            return std::string(buffer);
        })
        .doc()="Class PosCyl is used to store position in cylindrical coordinates.\nTakes in 3 aguments: R,z,phi.";
    py::class_<coord::VelCyl>(m,"VelCyl")
        .def(py::init<double, double, double >(),"vR"_a,"vz"_a,"vphi"_a)
        .def_readwrite("vR", &coord::VelCyl::vR)
        .def_readwrite("vz", &coord::VelCyl::vz)
        .def_readwrite("vphi", &coord::VelCyl::vphi)
        .def("__repr__",[](const coord::VelCyl &v) {
            char buffer[100];
            std::snprintf(buffer,100,"VelCyl: ( %.8g, %.8g, %.8g)",v.vR,v.vz,v.vphi);
            return std::string(buffer);
        })
        .doc()="Class VelCyl is used to store velocity in cylindrical coordinates.\nTakes in 3 aguments: vR,vz,vphi.";
    py::class_<coord::PosVelCyl>(m,"PosVelCyl")
        .def(py::init<double, double, double, double, double, double>(),"R"_a,"z"_a,"phi"_a,"vR"_a,"vz"_a,"vphi"_a)
        .def_readwrite("R", &coord::PosVelCyl::R)
        .def_readwrite("z", &coord::PosVelCyl::z)
        .def_readwrite("phi", &coord::PosVelCyl::phi)
        .def_readwrite("vR", &coord::PosVelCyl::vR)
        .def_readwrite("vz", &coord::PosVelCyl::vz)
        .def_readwrite("vphi", &coord::PosVelCyl::vphi)
        .def("__repr__",[](const coord::PosVelCyl &Rv) {
            char buffer[100];
            std::snprintf(buffer,100,"PosVelCyl: ( %.8g, %.8g, %.8g, %.8g, %.8g, %.8g)",
                Rv.R,Rv.z,Rv.phi,Rv.vR,Rv.vz,Rv.vphi);
            return std::string(buffer);
        })
        .doc()="Class PosVelCyl is used to store position and momentum in cylindrical coordinates.\n"
        "Takes in 6 aguments: R, z, phi, vR, vz, vphi as floats.";
    py::class_<coord::PosMomCyl>(m,"PosMomCyl")
        .def(py::init<double, double, double, double, double, double>(),"R"_a,"z"_a,"phi"_a,"pR"_a,"pz"_a,"pphi"_a)
        .def_readwrite("R", &coord::PosMomCyl::R)
        .def_readwrite("z", &coord::PosMomCyl::z)
        .def_readwrite("phi", &coord::PosMomCyl::phi)
        .def_readwrite("pR", &coord::PosMomCyl::pR)
        .def_readwrite("pz", &coord::PosMomCyl::pz)
        .def_readwrite("pphi", &coord::PosMomCyl::pphi)
        .def("__repr__",[](const coord::PosMomCyl &Rp) {
            char buffer[100];
            std::snprintf(buffer,100,"PosMomCyl: ( %.8g, %.8g, %.8g, %.8g, %.8g, %.8g)",
                Rp.R,Rp.z,Rp.phi,Rp.pR,Rp.pz,Rp.pphi);
            return std::string(buffer);
        })
        .doc()="Class PosMomCyl is used to store position and momentum in cylindrical coordinates.\n"
        "Takes in 6 aguments: R, z, phi, pR, pz, pphi as floats.";
    py::class_<coord::PosSph>(m,"PosSph")
        .def(py::init<double, double, double >(),"r"_a,"theta"_a,"phi"_a)
        .def_readwrite("r", &coord::PosSph::r)
        .def_readwrite("theta", &coord::PosSph::theta)
        .def_readwrite("phi", &coord::PosSph::phi)
        .def("__repr__",[](const coord::PosSph &r) {
            char buffer[100];
            std::snprintf(buffer,100,"PosSph: ( %.8g, %.8g, %.8g)",r.r,r.theta,r.phi);
            return std::string(buffer);
        })
        .doc()="Class PosSph is used to store position in spherical coordinates.\nTakes in 3 aguments: r, theta, phi.";
    py::class_<coord::VelSph>(m,"VelSph")
        .def(py::init<double, double, double >(),"vr"_a,"vtheta"_a,"vphi"_a)
        .def_readwrite("vr", &coord::VelSph::vr)
        .def_readwrite("vtheta", &coord::VelSph::vtheta)
        .def_readwrite("vphi", &coord::VelSph::vphi)
        .def("__repr__",[](const coord::VelSph &vr) {
            char buffer[100];
            std::snprintf(buffer,100,"VelSph: ( %.8g, %.8g, %.8g)", vr.vr,vr.vtheta,vr.vphi);
            return std::string(buffer);
        })
        .doc()="Class VelSph is used to store velocity in spherical coordinates.\nTakes in 3 aguments: vr, vtheta, vphi.";
    py::class_<coord::PosVelSph>(m,"PosVelSph")
        .def(py::init<double, double, double, double, double, double>(),"r"_a,"theta"_a,"phi"_a,"vr"_a,"vtheta"_a,"vphi"_a)
        .def_readwrite("r", &coord::PosVelSph::r)
        .def_readwrite("theta", &coord::PosVelSph::theta)
        .def_readwrite("phi", &coord::PosVelSph::phi)
        .def_readwrite("vr", &coord::PosVelSph::vr)
        .def_readwrite("vtheta", &coord::PosVelSph::vtheta)
        .def_readwrite("vphi", &coord::PosVelSph::vphi)
        .def("__repr__",[](const coord::PosVelSph &rv) {
            char buffer[100];
            std::snprintf(buffer,100,"PosVelSph: ( %.8g, %.8g, %.8g, %.8g, %.8g, %.8g)",
                rv.r,rv.theta,rv.phi,rv.vr,rv.vtheta,rv.vphi);
            return std::string(buffer);
        })
        .doc()="Class PosVelSph is used to store position and momentum in spherical coordinates.\n" 
        "Takes in 6 aguments: r, theta, phi, pr, ptheta, pphi as floats.";
    py::class_<coord::PosMomSph>(m,"PosMomSph")
        .def(py::init<double, double, double, double, double, double>(),"r"_a,"theta"_a,"phi"_a,"pr"_a,"ptheta"_a,"pphi"_a)
        .def_readwrite("r", &coord::PosMomSph::r)
        .def_readwrite("theta", &coord::PosMomSph::theta)
        .def_readwrite("phi", &coord::PosMomSph::phi)
        .def_readwrite("pr", &coord::PosMomSph::pr)
        .def_readwrite("ptheta", &coord::PosMomSph::ptheta)
        .def_readwrite("pphi", &coord::PosMomSph::pphi)
        .def("__repr__",[](const coord::PosMomSph &rp) {
            char buffer[100];
            std::snprintf(buffer,100,"PosMomSph: ( %.8g, %.8g, %.8g, %.8g, %.8g, %.8g)",
                rp.r,rp.theta,rp.phi,rp.pr,rp.ptheta,rp.pphi);
            return std::string(buffer);
        })
        .doc()="Class PosMomSph is used to store position and momentum in spherical coordinates.\n"
        "Takes in 6 aguments: R, z, phi, pR, pz, pphi as floats.";
    py::class_<coord::GradCar>(m,"GradCar")
        .def(py::init<double, double, double>(),"dx"_a,"dy"_a,"dz"_a)
        .def_readwrite("dx", &coord::GradCar::dx)
        .def_readwrite("dy", &coord::GradCar::dy)
        .def_readwrite("dz", &coord::GradCar::dz)
        .def("__repr__",[](const coord::GradCar &dx) {
            char buffer[100];
            std::snprintf(buffer,100,"GradCar: ( %.8g, %.8g, %.8g)",
                dx.dx,dx.dy,dx.dz);
            return std::string(buffer);
        })
        .doc()="Class GradCar stores gradient of a scalar function in cartesian coordinates\n"
        "Takes in 3 arguments: dx, dy, dz as floats.";
    py::class_<coord::GradCyl>(m,"GradCyl")
        .def(py::init<double, double, double>(),"dR"_a,"dz"_a,"dphi"_a)
        .def_readwrite("dR", &coord::GradCyl::dR)
        .def_readwrite("dphi", &coord::GradCyl::dphi)
        .def_readwrite("dz", &coord::GradCyl::dz)
        .def("__repr__",[](const coord::GradCyl &dR) {
            char buffer[100];
            std::snprintf(buffer,100,"GradCyl: ( %.8g, %.8g, %.8g)",
                dR.dR,dR.dz,dR.dphi);
            return std::string(buffer);
        })
        .doc()="Class GradCyl stores gradient of a scalar function in cylindrical coordinates\n"
        "Takes in 3 arguments: dR, dz, dphi as floats.";
    py::class_<coord::GradSph>(m,"GradSph")
        .def(py::init<double, double, double>(),"dr"_a,"dtheta"_a,"dphi"_a)
        .def_readwrite("dr", &coord::GradSph::dr)
        .def_readwrite("dtheta", &coord::GradSph::dtheta)
        .def_readwrite("dphi", &coord::GradSph::dphi)
        .def("__repr__",[](const coord::GradSph &dr) {
            char buffer[100];
            std::snprintf(buffer,100,"GradSph: ( %.8g, %.8g, %.8g)",
                dr.dr,dr.dtheta,dr.dphi);
            return std::string(buffer);
        })
        .doc()="Class GradSph stores gradient of a scalar function in spherical coordinates\n"
        "Takes in 3 arguments: dr, dtheta, dphi as floats.";
    py::class_<coord::HessCar>(m,"HessCar")
        .def(py::init<>())
        .def_readwrite("dx2", &coord::HessCar::dx2)
        .def_readwrite("dxdy", &coord::HessCar::dxdy)
        .def_readwrite("dxdz", &coord::HessCar::dxdz)
        .def_readwrite("dy2", &coord::HessCar::dy2)
        .def_readwrite("dydz", &coord::HessCar::dydz)
        .def_readwrite("dz2", &coord::HessCar::dz2)
        .doc()="Class HessCar stores Hessian of a scalar function in cartesian coordinates - d2Fdx^2, d2Fdydx etc.";
    py::class_<coord::HessCyl>(m,"HessCyl")
        .def(py::init<>())
        .def_readwrite("dR2", &coord::HessCyl::dR2)
        .def_readwrite("dRdphi", &coord::HessCyl::dRdphi)
        .def_readwrite("dRdz", &coord::HessCyl::dRdz)
        .def_readwrite("dphi2", &coord::HessCyl::dphi2)
        .def_readwrite("dzdphi", &coord::HessCyl::dzdphi)
        .def_readwrite("dz2", &coord::HessCyl::dz2)
        .doc()="Class HessCyl stores Hessian of a scalar function in cylindrical coordinates - d2FdR^2, d2FdRdz etc.";
    py::class_<coord::HessSph>(m,"Hesssph")
        .def(py::init<>())
        .def_readwrite("dr2", &coord::HessSph::dr2)
        .def_readwrite("drdtheta", &coord::HessSph::drdtheta)
        .def_readwrite("drdphi", &coord::HessSph::drdphi)
        .def_readwrite("dtheta2", &coord::HessSph::dtheta2)
        .def_readwrite("dthetadphi", &coord::HessSph::dthetadphi)
        .def_readwrite("dphi2", &coord::HessSph::dphi2)
        .doc()="Class HessCyl stores Hessian of a scalar function in spherical coordinates - d2Fdr^2, d2Fdrdtheta etc.";
    py::class_<obs::PosSky>(m,"PosSky")
        .def(py::init<double, double, bool>(),"_l"_a,"_b"_a,"_is_ra"_a=false)
        .def_readwrite("b", &obs::PosSky::b)
        .def_readwrite("l", &obs::PosSky::l)
        .def_readwrite("is_ra", &obs::PosSky::is_ra);
    py::class_<obs::VelSky>(m,"VelSky")
        .def(py::init<double, double, bool>(),"_mul"_a,"_mub"_a,"_is_ra"_a=false)
        .def_readwrite("mul", &obs::VelSky::mul)
        .def_readwrite("mub", &obs::VelSky::mub)
        .def_readwrite("is_ra", &obs::VelSky::is_ra);
    py::class_<obs::PosVelSky>(m,"PosVelSky")
        .def(py::init<obs::PosSky, obs::VelSky>())
        .def_readwrite("pm", &obs::PosVelSky::pm)
        .def_readwrite("pos", &obs::PosVelSky::pos)
        .def_readwrite("is_ra", &obs::PosVelSky::is_ra);
    py::class_<coord::Vel2Cyl>(m,"Vel2Cyl")
        .def(py::init<>())
        .def_readwrite("vphi2",&coord::Vel2Cyl::vphi2)
        .def_readwrite("vR2",&coord::Vel2Cyl::vR2)
        .def_readwrite("vz2",&coord::Vel2Cyl::vz2)
        .def_readwrite("vRvphi",&coord::Vel2Cyl::vRvphi)
        .def_readwrite("vRvz",&coord::Vel2Cyl::vRvz)
        .def_readwrite("vzvphi",&coord::Vel2Cyl::vzvphi);
    py::class_<coord::Vel2Car>(m,"Vel2Car")
        .def(py::init<>())
        .def_readwrite("vx2",&coord::Vel2Car::vx2)
        .def_readwrite("vy2",&coord::Vel2Car::vy2)
        .def_readwrite("vz2",&coord::Vel2Car::vz2)
        .def_readwrite("vxvy",&coord::Vel2Car::vxvy)
        .def_readwrite("vxvz",&coord::Vel2Car::vxvz)
        .def_readwrite("vyvz",&coord::Vel2Car::vyvz);
    py::class_<coord::Vel2Sph>(m,"Vel2Sph")
        .def(py::init<>())
        .def_readwrite("vphi2",&coord::Vel2Sph::vphi2)
        .def_readwrite("vr2",&coord::Vel2Sph::vr2)
        .def_readwrite("vtheta2",&coord::Vel2Sph::vtheta2)
        .def_readwrite("vrvtheta",&coord::Vel2Sph::vrvtheta)
        .def_readwrite("vrvphi",&coord::Vel2Sph::vrvphi)
        .def_readwrite("vthetavphi",&coord::Vel2Sph::vthetavphi);
    py::class_<obs::solarShifter>ss (m,"solarShifter");
    ss.def(py::init([](const units::InternalUnits &intUnits, coord::PosVelCar Vsun=coord::PosVelCar(NAN,NAN,NAN,NAN,NAN,NAN))
        { 
            if(Vsun.vx==Vsun.vx&&Vsun.vy==Vsun.vy&&Vsun.vz==Vsun.vz) return obs::solarShifter(intUnits,&Vsun);
            return obs::solarShifter(intUnits);
        }),"intUnits"_a,"Vsun"_a=coord::PosVelCar(NAN,NAN,NAN,NAN,NAN,NAN));
    ss.def("sKpc",[](obs::solarShifter &self,coord::PosCar p){return self.sKpc(p);});
    ss.def("sKpc",[](obs::solarShifter &self,coord::PosCyl p){return self.sKpc(p);});
    ss.def("toCar",[](obs::solarShifter &self,obs::PosSky pos,double sKpc){return self.toCar(pos,sKpc);});
    ss.def("toCar",[](obs::solarShifter &self,const obs::PosSky pos,double sKpc,
            const obs::VelSky pm,double Vlos_kms){return self.toCar(pos,sKpc,pm,Vlos_kms);});
    ss.def("toCyl",[](obs::solarShifter &self,obs::PosSky pos,double sKpc,obs::VelSky pm,double Vlos_kms){
            return self.toCyl(pos,sKpc,pm,Vlos_kms);});
    ss.def("toCyl",[](obs::solarShifter &self,const obs::PosSky pos,double sKpc,
            const obs::VelSky pm,double Vlos_kms){return self.toCyl(pos,sKpc,pm,Vlos_kms);});
    ss.def("toPM",[](obs::solarShifter &self,const coord::PosVelCyl pv,double &Vlos_kms)
        {return self.toPM(pv,Vlos_kms);});
    ss.def("toPM",[](obs::solarShifter &self,const coord::PosVelCar pv,double &Vlos_kms)
        {return self.toPM(pv,Vlos_kms);});
    ss.def("toSky",[](obs::solarShifter &self,const coord::PosVelCar pv,double &sKpc,double &Vlos_kms)
        {return self.toSky(pv,sKpc,Vlos_kms);});
    ss.def("toSky",[](obs::solarShifter &self,const coord::PosVelCyl pv,double &sKpc,double &Vlos_kms)
        {return self.toSky(pv,sKpc,Vlos_kms);});
    ss.def("toSky",[](obs::solarShifter &self,const coord::PosCar p,double &sKpc)
        {return self.toSky(p,sKpc);});
    ss.def("toSky",[](obs::solarShifter &self,const coord::PosCyl p,double &sKpc)
        {return self.toSky(p,sKpc);});
    ss.def("Vxyz",&obs::solarShifter::Vxyz);
    ss.def("xyz",&obs::solarShifter::xyz);
    ss.def_readonly("from_Kpc", &obs::solarShifter::from_Kpc);
    ss.def_readonly("from_kms", &obs::solarShifter::from_kms);
    ss.def_readonly("from_mas_per_yr", &obs::solarShifter::from_mas_per_yr);
    ss.def_readonly("torad", &obs::solarShifter::torad);
    py::class_<potential::BaseDensity,std::shared_ptr<potential::BaseDensity>>(m,"BaseDensity")
        .def("enclosedMass",&potential::BaseDensity::enclosedMass)
        .def("name",&potential::BaseDensity::name)
        .def("totalMass",&potential::BaseDensity::totalMass)
        .def("density",[](potential::BaseDensity &self,coord::PosCar pos){return self.density(pos);})
        .def("density",[](potential::BaseDensity &self,coord::PosCyl pos){return self.density(pos);})
        .def("density",[](potential::BaseDensity &self,coord::PosSph pos){return self.density(pos);})
        .doc()="Class BaseDensity defines a density profile without a corresponding potential.\n"
        "It enables computation of the density of a point in cartesian, cylindrical and spherical coordinates.";
    py::class_<potential::BasePotential,std::shared_ptr<potential::BasePotential>>(m, "BasePotential")
        .def("totalMass",&potential::BasePotential::totalMass)
        .def("enclosedMass",&potential::BasePotential::enclosedMass)
        .def("name",&potential::BasePotential::name)
        .def("value",&potential::BasePotential::value<coord::Car>)
        .def("value",&potential::BasePotential::value<coord::Cyl>)
        .def("eval",[](potential::BasePotential &self, coord::PosCar x,bool pot=false,bool der=false,bool hess=false)
        ->py::object{
            if(!pot&!der&&!hess)pot=true;
            double pot0;
            coord::GradCar ders;
            coord::HessCar hess1;
            self.eval(x,pot?&pot0:NULL,der?&ders:NULL,hess?&hess1:NULL);
            if((pot&&der)||(pot&&hess)||(der&&hess)){
                py::list ls2;
                if(pot)ls2.append(pot0);
                if(der)ls2.append(ders);
                if(hess)ls2.append(hess1);
                return ls2;
            }
            if(der) return py::cast(ders);
            if(hess) return py::cast(hess1); 
            return py::cast(pot0);
        },"x"_a,"pot"_a=false,"der"_a=false,"hess"_a=false)
        .def("eval",[](potential::BasePotential &self, coord::PosCyl x,bool pot=false,bool der=false,bool hess=false)
        -> py::object{
            if(!pot&!der&&!hess)pot=true;
            double pot0;
            coord::GradCyl ders;
            coord::HessCyl hess1;
            self.eval(x,pot?&pot0:NULL,der?&ders:NULL,hess?&hess1:NULL);
            if((pot&&der)||(pot&&hess)||(der&&hess)){
                py::list ls2;
                if(pot)ls2.append(pot0);
                if(der)ls2.append(ders);
                if(hess)ls2.append(hess1);
                return ls2;
            }
            if(der) return py::cast(ders);
            if(hess) return py::cast(hess1); 
            return py::cast(pot0);
        },"x"_a,"pot"_a=false,"der"_a=false,"hess"_a=false)
        .def("eval",[](potential::BasePotential &self, coord::PosSph x,bool pot=false,bool der=false,bool hess=false)
        -> py::object{
            if(!pot&!der&&!hess)pot=true;
            double pot0;
            coord::GradSph ders;
            coord::HessSph hess1;
            self.eval(x,pot?&pot0:NULL,der?&ders:NULL,hess?&hess1:NULL);
            if((pot&&der)||(pot&&hess)||(der&&hess)){
                py::list ls2;
                if(pot)ls2.append(pot0);
                if(der)ls2.append(ders);
                if(hess)ls2.append(hess1);
                return ls2;
            }
            if(der) return py::cast(ders);
            if(hess) return py::cast(hess1); 
            return py::cast(pot0);
        },"x"_a,"pot"_a=false,"der"_a=false,"hess"_a=false)
        .doc()="Class BasePotential defines a gravitational potential.\n"
        "It enables computation of the potential and derivative of a point in cartesian, cylindrical and spherical coordinates.";
    py::class_<actions::Actions>(m,"Actions")
        .def(py::init<double,double,double>(),"Jr"_a,"Jz"_a,"Jphi"_a)
        .def_readwrite("Jr", &actions::Actions::Jr)
        .def_readwrite("Jz", &actions::Actions::Jz)
        .def_readwrite("Jphi", &actions::Actions::Jphi)
        .def("__repr__",[](const actions::Actions &J) {
            char buffer[100];
            std::snprintf(buffer,100,"Actions: ( %.8g, %.8g, %.8g)",
                J.Jr,J.Jz,J.Jphi);
            return std::string(buffer);
        })
        .doc()="Class Actions stores value of the actions.\n"
        "Takes in 3 arguments: Jr, Jz, Jphi as floats.";
    py::class_<actions::Angles>(m,"Angles")
        .def(py::init<double,double,double>(),"thetar"_a,"thetaz"_a,"thetaphi"_a)
        .def_readwrite("thetar", &actions::ActionAngles::thetar)
        .def_readwrite("thetaz", &actions::ActionAngles::thetaz)
        .def_readwrite("thetaphi", &actions::ActionAngles::thetaphi)
        .def("__repr__",[](const actions::Angles &theta) {
            char buffer[100];
            std::snprintf(buffer,100,"Angles: ( %.8g, %.8g, %.8g)",
                theta.thetar,theta.thetaz,theta.thetaphi);
            return std::string(buffer);
        })
        .doc()="Class Angles stores value of the angle coordinates.\n"
        "Takes in 3 arguments: thetar, thetaz, thetaphi as floats.";
    py::class_<actions::ActionAngles>(m,"ActionAngles")
        .def(py::init<actions::Actions,actions::Angles>(),"Actions"_a,"Angles"_a)
        .def(py::init<double,double,double,double,double,double>(),"Jr"_a,"Jz"_a,"Jphi"_a,"thetar"_a,"thetaz"_a,"thetaphi"_a)
        .def_readwrite("Jr", &actions::ActionAngles::Jr)
        .def_readwrite("Jz", &actions::ActionAngles::Jz)
        .def_readwrite("Jphi", &actions::ActionAngles::Jphi)
        .def_readwrite("thetar", &actions::ActionAngles::thetar)
        .def_readwrite("thetaz", &actions::ActionAngles::thetaz)
        .def_readwrite("thetaphi", &actions::ActionAngles::thetaphi)
        .def("__repr__",[](const actions::ActionAngles &aa) {
            char buffer[200];
            std::snprintf(buffer,200,"ActionAngles: ( %.8g, %.8g, %.8g, %.8g, %.8g, %.8g)",
                aa.Jr,aa.Jz,aa.Jphi,aa.thetar,aa.thetaz,aa.thetaphi);
            return std::string(buffer);
        })
        .doc()="Class ActionAngles stores value of the actions and angle coordinates.\n"
        "Takes in 6 arguments: Jr, Jz, Jphi, thetar, thetaz, thetaphi as floats.";
    py::class_<actions::Frequencies>(m,"Frequencies")
        .def(py::init<double,double,double>(), "Omegar"_a, "Omegaz"_a,"Omegaphi"_a)
        .def_readwrite("Omegar", &actions::Frequencies::Omegar)
        .def_readwrite("Omegaz", &actions::Frequencies::Omegaz)
        .def_readwrite("Omegaphi", &actions::Frequencies::Omegaphi)
        .def("__repr__",[](const actions::Frequencies &fr) {
            char buffer[100];
            std::snprintf(buffer,100,"Frequencies: ( %.8g, %.8g, %.8g)",
                fr.Omegar,fr.Omegaz,fr.Omegaphi);
            return std::string(buffer);
        })
        .doc()="Class Frequencies stores value of the frequencies.\n"
        "Takes in 3 arguments: Omegar, Omegaz, Omegaphi as floats.";
    py::class_<actions::TorusGenerator>(m,"TorusGenerator")
        .def(py::init([](potential::PtrPotential pot, const double  tol=1e-9) 
        { return actions::TorusGenerator(*pot,tol);}),"Potential"_a,"tolerance"_a=1e-9)
        .def("fitTorus",[] (actions::TorusGenerator& self,actions::Actions J, double tighten = 1, int type = 0){
            int typer=type;
            if(type>2)typer=2;
            if(type<0)typer=0;
            return self.fitTorus(J,tighten,actions::ToyPotType(type));
        },"J"_a,"tighten"_a=1,"type"_a=0)
        .doc()="Class Torus Generator is used to make Tori with fit torus function.\n"
        "Initialised with a potential as well as optionally the tolerance which sets how much hamiltonian can vary along torus.";
    py::class_<actions::Jzfind>(m,"Jzfind")
        .def(py::init<const potential::BasePotential&>(),"Potential"_a)
        .def_readwrite("parameters",&actions::Jzfind::params)
        .def("value",&actions::Jzfind::value)
        .doc()="Class Jzfind is used to find the value of Jz at the box loop orbit transition. This is done by a polynomial fit.\n"
        "Initialised with a potential. Function value return critical Jz given Jr and parameters is parameters of the polynomial.";
    py::class_<actions::Torus>(m,"Torus")
		.def("from_true",&actions::Torus::from_true,"Given true angle coordinates, gives position and momentum.")
	    .def("from_toy",&actions::Torus::from_toy,"Given toy angle coordinates, gives position and momentum.")
	    .def("Omega",&actions::Torus::Omega, "Gives frequency.")
	    .def("density",&actions::Torus::density, "Gives Density")
	    .def("orbit",&actions::Torus::orbit,"Given true angle coordinates, returns an array containing the position and momentum at each time from the constructed torus.")
        .doc()="Class Torus is used to represent a constructed torus. Initialised with TorusGenerator's fit torus function\n" ;
    py::class_<actions::BaseActionFinder,std::shared_ptr<actions::BaseActionFinder>>(m,"BaseActionFinder")
        .def("actionAngles",[] (actions::BaseActionFinder &self, coord::PosVelCyl xv,bool freq=false)
        -> py::object 
        { 
            if(!freq)return py::cast(self.actionAngles(xv));
            actions::Frequencies freqs;
            py::list ls;
            ls.append(self.actionAngles(xv,&freqs));
            ls.append(freqs);
            return ls;
        },"xv"_a,"freq"_a=false )
        .def("actions",&actions::BaseActionFinder::actions);
     py::class_<actions::ActionFinderSpherical,std::shared_ptr<actions::ActionFinderSpherical>,actions::BaseActionFinder>(m,"ActionFinderSpherical")
        .def("actionAngles",[] (actions::ActionFinderSpherical &self, coord::PosVelCyl xv,bool freq=false)
        -> py::object
        { 
            if(!freq)return py::cast(self.actionAngles(xv));
            actions::Frequencies freqs;
            py::list ls;
            ls.append(self.actionAngles(xv,&freqs));
            ls.append(freqs);
            return ls;
        },"xv"_a,"freq"_a=false )
        .def("actions",&actions::ActionFinderSpherical::actions);
    py::class_<actions::ActionFinderAxisymFudge,std::shared_ptr<actions::ActionFinderAxisymFudge>,actions::BaseActionFinder>(m,"ActionFinderAxisymFudge")
		    .def(py::init<const potential::PtrPotential&,bool>(),"potential"_a,"interpolate"_a=false)
        .def("actionAngles",[] (actions::ActionFinderAxisymFudge &self, coord::PosVelCyl xv,bool freq=false)->py::object
        { 
            if(!freq)return py::cast(self.actionAngles(xv));
            actions::Frequencies freqs;
            py::list ls;
            ls.append(self.actionAngles(xv,&freqs));
            ls.append(freqs);
            return ls;
        },"xv"_a,"freq"_a=false)
        .def("actions",&actions::ActionFinderAxisymFudge::actions);
    py::class_<actions::ActionFinderTG,std::shared_ptr<actions::ActionFinderTG>,actions::BaseActionFinder>(m,"ActionFinderTG")
	    .def(py::init<const potential::PtrPotential&,
			 const actions::TorusGenerator&>())
        .def("actionAngles",[] (actions::ActionFinderTG &self, coord::PosVelCyl xv,bool freq=false)->py::object
        { 
            if(!freq)return py::cast(self.actionAngles(xv));
            actions::Frequencies freqs;
            py::list ls;
            ls.append(self.actionAngles(xv,&freqs));
            ls.append(freqs);
            return ls;
        },"xv"_a,"freq"_a=false )
        .def("actions",&actions::ActionFinderTG::actions);
    py::class_<potential::CompositeDensity,std::shared_ptr<potential::CompositeDensity>,potential::BaseDensity>(m,"CompositeDensity")
        .def(py::init<const std::vector<potential::PtrDensity>&>())
        .def("enclosedMass",&potential::BaseDensity::enclosedMass)
        .def("name",&potential::BaseDensity::name)
        .def("totalMass",&potential::BaseDensity::totalMass)
        .def("density",[](potential::BaseDensity &self,coord::PosCar pos){return self.density(pos);})
        .def("density",[](potential::BaseDensity &self,coord::PosCyl pos){return self.density(pos);})
        .def("density",[](potential::BaseDensity &self,coord::PosSph pos){return self.density(pos);})
        .def("component",&potential::CompositeDensity::component);
    py::class_<potential::Multipole,std::shared_ptr<potential::Multipole>,potential::BasePotential>(m,"Multipole")
        .def(py::init<std::vector<double>&,std::vector<std::vector<double>>&,std::vector<std::vector<double>>&>())
        .def("density",[](potential::Multipole &self,coord::PosCar pos){return self.density(pos);})
        .def("density",[](potential::Multipole &self,coord::PosCyl pos){return self.density(pos);})
        .def("density",[](potential::Multipole &self,coord::PosSph pos){return self.density(pos);})
        .def("enclosedMass",&potential::Multipole::enclosedMass)
        .def("name",&potential::Multipole::name)
        .def("totalMass",&potential::Multipole::totalMass);
    py::class_<potential::CylSpline,std::shared_ptr<potential::CylSpline>,potential::BasePotential>(m,"CylSpline");
    py::class_<df::BaseDistributionFunction,std::shared_ptr<df::BaseDistributionFunction>>(m,"BaseDistributionFunction")
        .def("value",&df::BaseDistributionFunction::value)
        .def("numValues",&df::BaseDistributionFunction::numValues)
        .def("LF",&df::BaseDistributionFunction::LF)
        .def("readBrighterThan",&df::BaseDistributionFunction::readBrighterThan)
        .def("setNorm",&df::BaseDistributionFunction::set_norm)
        .def("selectMag",&df::BaseDistributionFunction::selectMag)
        .def("totalMass",[](df::BaseDistributionFunction &self,double reqRelError=1e-6,int MaxNumEval=1000000)
        {return self.totalMass(reqRelError,MaxNumEval);},"reqRelError"_a=1e-6,"MaxNumEval"_a=1000000)
        .def("epicycle_ratios",[](df::BaseDistributionFunction &self,const actions::Actions J){
            double r;
            self.epicycle_ratios(J,&r);
            return r;
        })
        .def("eval",[](df::BaseDistributionFunction &self,const double Jzcr, const actions::Actions J){
            int n=self.numValues();
            std::vector<double> vals(n);
            self.eval(J,Jzcr,&vals[0]);
            return vals;
        })
        .def("tab_params",[](df::BaseDistributionFunction &self,std::string filename,const units::InternalUnits &units){
            std::ofstream os(filename);
            self.tab_params(os,units);
        })
        .def("write_params",[](df::BaseDistributionFunction &self,std::string filename,const units::InternalUnits &units){
            std::ofstream os(filename);
            self.write_params(os,units);
        });
    py::class_<df::CompositeDF,std::shared_ptr<df::CompositeDF>,df::BaseDistributionFunction>(m,"CompositeDF")
        .def(py::init<const std::vector<df::PtrDistributionFunction>>())
        .def("component",&df::CompositeDF::component);
    py::class_<df::DoublePowerLawParam>(m,"DoublePowerLawParam")
        .def(py::init<>())
        .def_readwrite("coefJrIn",&df::DoublePowerLawParam::coefJrIn)
        .def_readwrite("coefJrOut",&df::DoublePowerLawParam::coefJrOut)
        .def_readwrite("coefJzIn",&df::DoublePowerLawParam::coefJzIn)
        .def_readwrite("coefJzOut",&df::DoublePowerLawParam::coefJzOut)
        .def_readwrite("cutoffStrength",&df::DoublePowerLawParam::cutoffStrength)
        .def_readwrite("Fname",&df::DoublePowerLawParam::Fname)
        .def_readwrite("J0",&df::DoublePowerLawParam::J0)
        .def_readwrite("Jcore",&df::DoublePowerLawParam::Jcore)
        .def_readwrite("Jcutoff",&df::DoublePowerLawParam::Jcutoff)
        .def_readwrite("Jphi0",&df::DoublePowerLawParam::Jphi0)
        .def_readwrite("norm",&df::DoublePowerLawParam::norm)
        .def_readwrite("rotFrac",&df::DoublePowerLawParam::rotFrac)
        .def_readwrite("slopeIn",&df::DoublePowerLawParam::slopeIn)
        .def_readwrite("slopeOut",&df::DoublePowerLawParam::slopeOut)
        .def_readwrite("steepness",&df::DoublePowerLawParam::steepness);
    py::class_<df::DoublePowerLaw,std::shared_ptr<df::DoublePowerLaw>,df::BaseDistributionFunction>(m,"DoublePowerLawDF")
        .def(py::init<df::DoublePowerLawParam>());
    py::class_<df::IsochroneParam>(m,"IsochroneParam")
        .def(py::init<>())
        .def_readwrite("scaleRadius",&df::IsochroneParam::scaleRadius)
        .def_readwrite("mass",&df::IsochroneParam::mass)
        .def_readwrite("nu",&df::IsochroneParam::nu)
        .def_readwrite("mu",&df::IsochroneParam::mu);
    py::class_<df::IsochroneDF,std::shared_ptr<df::IsochroneDF>,df::BaseDistributionFunction>(m,"IsochroneDF")
        .def(py::init<df::IsochroneParam>());
    py::class_<particles::ParticleArrayCar>(m,"ParticleArrayCar")
        .def(py::init<>())
        .def("add",&particles::ParticleArrayCar::add)
        .def("mass",&particles::ParticleArrayCar::mass)
        .def("point",&particles::ParticleArrayCar::point)
        .def("totalMass",&particles::ParticleArrayCar::totalMass)
        .def("size",&particles::ParticleArrayCar::size);
    py::class_<particles::ParticleArrayCyl>(m,"ParticleArrayCyl")
        .def(py::init<>())
        .def("add",&particles::ParticleArrayCyl::add)
        .def("mass",&particles::ParticleArrayCyl::mass)
        .def("point",&particles::ParticleArrayCyl::point)
        .def("totalMass",&particles::ParticleArrayCyl::totalMass)
        .def("size",&particles::ParticleArrayCyl::size);
    py::class_<particles::ParticleArraySph>(m,"ParticleArraySph")
        .def(py::init<>())
        .def("add",&particles::ParticleArraySph::add)
        .def("mass",&particles::ParticleArraySph::mass)
        .def("point",&particles::ParticleArraySph::point)
        .def("totalMass",&particles::ParticleArraySph::totalMass)
        .def("size",&particles::ParticleArraySph::size);
    py::class_<particles::ParticleArray<coord::PosCyl>>(m,"ParticleArrayPosCyl")
        .def(py::init<>())
        .def("add",&particles::ParticleArray<coord::PosCyl>::add)
        .def("mass",&particles::ParticleArray<coord::PosCyl>::mass)
        .def("point",&particles::ParticleArray<coord::PosCyl>::point)
        .def("totalMass",&particles::ParticleArray<coord::PosCyl>::totalMass)
        .def("size",&particles::ParticleArray<coord::PosCyl>::size);
    py::class_<particles::ParticleArray<coord::PosCar>>(m,"ParticleArrayPosCar")
        .def(py::init<>())
        .def("add",&particles::ParticleArray<coord::PosCar>::add)
        .def("mass",&particles::ParticleArray<coord::PosCar>::mass)
        .def("point",&particles::ParticleArray<coord::PosCar>::point)
        .def("totalMass",&particles::ParticleArray<coord::PosCar>::totalMass)
        .def("size",&particles::ParticleArray<coord::PosCar>::size);
    py::class_<particles::ParticleArray<coord::PosSph>>(m,"ParticleArrayPosSph")
        .def(py::init<>())
        .def("add",&particles::ParticleArray<coord::PosSph>::add)
        .def("mass",&particles::ParticleArray<coord::PosSph>::mass)
        .def("point",&particles::ParticleArray<coord::PosSph>::point)
        .def("totalMass",&particles::ParticleArray<coord::PosSph>::totalMass)
        .def("size",&particles::ParticleArray<coord::PosSph>::size);
    py::class_<galaxymodel::BaseComponent,std::shared_ptr<galaxymodel::BaseComponent>>(m,"BaseComponent")
        .def("getPotential",&galaxymodel::BaseComponent::getPotential)
        .def("getDensity",&galaxymodel::BaseComponent::getDensity)
        .def("update",&galaxymodel::BaseComponent::update)
        .def_readonly("isDensityDisklike",&galaxymodel::BaseComponent::isDensityDisklike);
    py::class_<galaxymodel::ComponentStatic,std::shared_ptr<galaxymodel::ComponentStatic>,galaxymodel::BaseComponent>(m,"ComponentStatic")
        .def(py::init<const potential::PtrPotential&>())
        .def(py::init<const potential::PtrDensity&,bool>());
    py::class_<galaxymodel::ComponentWithSpheroidalDF,std::shared_ptr<galaxymodel::ComponentWithSpheroidalDF>,galaxymodel::BaseComponent>(m,"ComponentWithSpheroidalDF")
        .def(py::init<const df::PtrDistributionFunction &,const potential::PtrDensity &,unsigned int,unsigned int,
            unsigned int,double,double,double,unsigned int>(),"df"_a,"initDensity"_a,"lmax"_a,"mmax"_a,"gridSizeR"_a,
            "rmin"_a,"rmax"_a,"relError"_a=(0.001),"maxNumEval"_a=100000U);
     py::class_<galaxymodel::ComponentWithDisklikeDF,std::shared_ptr<galaxymodel::ComponentWithDisklikeDF>,galaxymodel::BaseComponent>(m,"ComponentWithDisklikeDF")
        .def(py::init<const df::PtrDistributionFunction &,const potential::PtrDensity &,unsigned int,unsigned int,double,double,unsigned int,
           double,double,double,unsigned int>(),"df"_a,"initDensity"_a,"mmax"_a,"gridSizeR"_a,
            "rmin"_a,"rmax"_a,"gridSizez"_a,"zmin"_a,"zmax"_a,"relError"_a=(0.001),"maxNumEval"_a=100000U);
    py::class_<galaxymodel::GalaxyModel>(m,"GalaxyModel")
        .def(py::init<potential::BasePotential &, actions::BaseActionFinder &,df::BaseDistributionFunction &>())
        .def("potential",[](galaxymodel::GalaxyModel &self){return &self.potential;})
        .def("actFinder",[](galaxymodel::GalaxyModel &self){return &self.actFinder;})
        .def("distrFunc",[](galaxymodel::GalaxyModel &self){return &self.distrFunc;});
    py::class_<galaxymodel::SelfConsistentModel>(m,"SelfConsistentModel")
        .def(py::init<>())
        .def_readwrite("actionFinder",&galaxymodel::SelfConsistentModel::actionFinder)
        .def_readwrite("components",&galaxymodel::SelfConsistentModel::components)
        .def_readwrite("rminSph",&galaxymodel::SelfConsistentModel::rminSph)
        .def_readwrite("rmaxSph",&galaxymodel::SelfConsistentModel::rmaxSph)
        .def_readwrite("RminCyl",&galaxymodel::SelfConsistentModel::RminCyl)
        .def_readwrite("RmaxCyl",&galaxymodel::SelfConsistentModel::RmaxCyl)
        .def_readwrite("lmaxAngularSph",&galaxymodel::SelfConsistentModel::lmaxAngularSph)
        .def_readwrite("mmaxAngularSph",&galaxymodel::SelfConsistentModel::mmaxAngularSph)
        .def_readwrite("mmaxAngularCyl",&galaxymodel::SelfConsistentModel::mmaxAngularCyl)
        .def_readwrite("zmaxCyl",&galaxymodel::SelfConsistentModel::zmaxCyl)
        .def_readwrite("zminCyl",&galaxymodel::SelfConsistentModel::zminCyl)
        .def_readwrite("sizeRadialCyl",&galaxymodel::SelfConsistentModel::sizeRadialCyl)
        .def_readwrite("sizeRadialSph",&galaxymodel::SelfConsistentModel::sizeRadialSph)
        .def_readwrite("sizeVerticalCyl",&galaxymodel::SelfConsistentModel::sizeVerticalCyl)
        .def_readwrite("totalPotential",&galaxymodel::SelfConsistentModel::totalPotential)
        .def_readwrite("useActionInterpolation",&galaxymodel::SelfConsistentModel::useActionInterpolation);

    m.def("Vcirc",[](potential::PtrPotential pot, const double R)
    {return potential::v_circ(*pot, R);});
    m.def("Vcirc",[](potential::PtrPotential pot, std::vector<double> R)
    {
        py::list ls;
        for(int i=0;i<R.size();i++){
            ls.append(potential::v_circ(*pot, R[i]));
        }
        return ls;
    });

    m.def("toPosCar",&coord::toPosCyl<coord::Cyl>);
    m.def("toPosCar",&coord::toPosCyl<coord::Sph>);

    m.def("toPosCyl",&coord::toPosCyl<coord::Car>);
    m.def("toPosCyl",&coord::toPosCyl<coord::Sph>);

    m.def("toPosSph",&coord::toPosSph<coord::Car>);
    m.def("toPosSph",&coord::toPosSph<coord::Cyl>);

    m.def("toPosVelCar",[](coord::PosVelCyl xv){return coord::toPosVelCar(xv);});
    m.def("toPosVelCar",[](coord::PosVelSph xv){return coord::toPosVelCar(xv);});
    m.def("toPosVelCar",[](coord::PosMomCar xp){return coord::PosVelCar(xp.x,xp.y,xp.z,xp.px,xp.py,xp.pz);});

    m.def("toPosVelCyl",[](coord::PosVelCar xv){return coord::toPosVelCyl(xv);});
    m.def("toPosVelCyl",[](coord::PosVelSph xv){return coord::toPosVelCyl(xv);});
    m.def("toPosVelCyl",[](coord::PosMomCyl xp){return coord::toPosVelCyl(xp);});

    m.def("toPosVelSph",[](coord::PosVelCar xv){return coord::toPosVelSph(xv);});
    m.def("toPosVelSph",[](coord::PosVelCyl xv){return coord::toPosVelSph(xv);});
    //m.def("toPosVelSph",[](coord::PosMomSph xp){return coord::toPosVelSph(xp);});

    //m.def("toPosMomCar",[](coord::PosMomCyl xp){return coord::toPosMomCar(xp);});
    //m.def("toPosMomCar",[](coord::PosMomSph xp){return coord::toPosMomCar(xp);});
    m.def("toPosMomCar",[](coord::PosVelCar xv){return coord::toPosMomCar(xv);});

    //m.def("toPosMomCyl",[](coord::PosMomCar xp){return coord::toPosMomCyl(xp);});
    //m.def("toPosMomCyl",[](coord::PosMomSph xp){return coord::toPosMomCyl(xp);});
    m.def("toPosMomCyl",[](coord::PosVelCyl xv){return coord::toPosMomCyl(xv);});

    //m.def("toPosMomSph",[](coord::PosVelCar xp){return coord::toPosMomSph(xp);});
    //m.def("toPosMomSph",[](coord::PosVelCyl xp){return coord::toPosMomSph(xp);});
    //m.def("toPosMomSph",[](coord::PosMomSph xv){return coord::toPosMomSph(xv);});
    m.def("from_RAdec",[](obs::PosVelSky p){return obs::from_RAdec(p);});
    m.def("from_RAdec",[](obs::PosSky p){return obs::from_RAdec(p);});

//    m.def("createDensity",[](const std::string vals){return potential::PtrDensity(potential::createDensity(utils::KeyValueMap(vals)));});
    m.def("createPotential",[](const std::string vals){return potential::PtrPotential(potential::createPotential(utils::KeyValueMap(vals)));});
    m.def("createPotential",[](const utils::KeyValueMap &params,const particles::ParticleArray<coord::PosCyl> &Particles,const units::ExternalUnits &converter=units::ExternalUnits())
    {return potential::createPotential(params,Particles,converter);},"params"_a,"particles"_a,"converter"_a=units::ExternalUnits());
    m.def("createPotential",[](const utils::KeyValueMap &params,const units::ExternalUnits &converter=units::ExternalUnits())
    {return potential::createPotential(params,converter);},"params"_a,"converter"_a=units::ExternalUnits());
    m.def("createPotential",[](const utils::KeyValueMap &params,const potential::BasePotential &pot, const units::ExternalUnits &converter=units::ExternalUnits())
    {return potential::createPotential(params,pot,converter);},"params"_a,"pot"_a,"converter"_a=units::ExternalUnits());
    m.def("createPotential",[](const utils::KeyValueMap &params,const potential::BaseDensity &dens, const units::ExternalUnits &converter=units::ExternalUnits())
    {return potential::createPotential(params,dens,converter);},"params"_a,"dens"_a,"converter"_a=units::ExternalUnits());
    m.def("createPotential",[](const std::vector<utils::KeyValueMap> &params,const units::ExternalUnits &converter=units::ExternalUnits())
    {return potential::createPotential(params,converter);},"params"_a,"converter"_a=units::ExternalUnits());
    m.def("createMultipole",[](const potential::BasePotential &src,int lmax,int mmax,int gridSizeR,double rmin=(0.0),double rmax=(0.0))
    {return potential::Multipole::create(src,lmax,mmax,gridSizeR,rmin,rmax);}
    ,"src"_a,"lmax"_a,"mmax"_a,"gridSizeR"_a,"rmin"_a=(0.0),"rmax"_a=(0.0));
    m.def("createMultipole",[](const potential::BaseDensity &src,int lmax,int mmax,int gridSizeR,double rmin=(0.0),double rmax=(0.0))
    {return potential::Multipole::create(src,lmax,mmax,gridSizeR,rmin,rmax);}
    ,"src"_a,"lmax"_a,"mmax"_a,"gridSizeR"_a,"rmin"_a=(0.0),"rmax"_a=(0.0));
    m.def("createDensity",&potential::createDensity,"params"_a,"converter"_a=units::ExternalUnits());
    m.def("integrateTraj",[] (coord::PosVelCyl xv,double T,double dt,potential::PtrPotential pot){
        return orbit::integrateTraj(xv,T,dt,*pot);
    } );
    m.def("interpTorus",[](double x,actions::Torus T0,actions::Torus T1, actions::TorusGenerator TG){
        return actions::interpTorus(x,T0,T1,&TG);
    },"x"_a,"T0"_a,"T1"_a,"TorusGenerator"_a,"function to interpolate between two tori, T=x*T0+(1-x)*T1\n");  
    m.def("createDistributionFunction",&df::createDistributionFunction,"params"_a,"potential"_a,"density"_a=NULL,"converter"_a=units::ExternalUnits());
    m.def("createDistributionFunction",[](const utils::KeyValueMap &params,potential::BasePotential *potential,units::ExternalUnits converter=units::ExternalUnits())
      {return df::createDistributionFunction(params,potential,NULL,converter);}
      ,"params"_a,"potential"_a,"converter"_a=units::ExternalUnits());
    m.def("sampleDensity",[](const potential::BaseDensity &dens,int numPoints)
    {return galaxymodel::sampleDensity(dens,numPoints);});
    m.def("samplePosVel",[](const galaxymodel::GalaxyModel &model,int numPoints)
    {return galaxymodel::samplePosVel(model,numPoints);});
    m.def("sampleVelocity",[](const galaxymodel::GalaxyModel &model,int numPoints,coord::PosCyl x)
    {return galaxymodel::sampleVelocity(model,x,numPoints);});
    m.def("writeSnapshot",[](const std::string &fileName, const particles::ParticleArray<coord::PosCar> &particles,
        const std::string &fileFormat="Text",
        const units::ExternalUnits &unitConverter=units::ExternalUnits(),const std::string &header="", 
        const double time=-((float)((float)((1E300)*(1E300)*(0.0F)))), const bool append=false)
        {particles::writeSnapshot(fileName,particles,fileFormat,unitConverter,header,time,append);},"fileName"_a,"particles"_a,"fileFormat"_a="Text",
        "unitConverter"_a=units::ExternalUnits(),"header"_a="","time"_a=-((float)((float)((1E300)*(1E300)*(0.0F)))),"append"_a=false);
    m.def("writeSnapshot",[](const std::string &fileName, const particles::ParticleArray<coord::PosCyl> &particles,
        const std::string &fileFormat="Text",
        const units::ExternalUnits &unitConverter=units::ExternalUnits(),const std::string &header="", 
        const double time=-((float)((float)((1E300)*(1E300)*(0.0F)))), const bool append=false)
        {particles::writeSnapshot(fileName,particles,fileFormat,unitConverter,header,time,append);},"fileName"_a,"particles"_a,"fileFormat"_a="Text",
        "unitConverter"_a=units::ExternalUnits(),"header"_a="","time"_a=-((float)((float)((1E300)*(1E300)*(0.0F)))),"append"_a=false);
    m.def("writeSnapshot",[](const std::string &fileName, const particles::ParticleArray<coord::PosSph> &particles,
        const std::string &fileFormat="Text",
        const units::ExternalUnits &unitConverter=units::ExternalUnits(),const std::string &header="", 
        const double time=-((float)((float)((1E300)*(1E300)*(0.0F)))), const bool append=false)
        {particles::writeSnapshot(fileName,particles,fileFormat,unitConverter,header,time,append);},"fileName"_a,"particles"_a,"fileFormat"_a="Text",
        "unitConverter"_a=units::ExternalUnits(),"header"_a="","time"_a=-((float)((float)((1E300)*(1E300)*(0.0F)))),"append"_a=false);
     m.def("writeSnapshot",[](const std::string &fileName, const particles::ParticleArray<coord::PosVelCar> &particles,
        const std::string &fileFormat="Text",
        const units::ExternalUnits &unitConverter=units::ExternalUnits(),const std::string &header="", 
        const double time=-((float)((float)((1E300)*(1E300)*(0.0F)))), const bool append=false)
        {particles::writeSnapshot(fileName,particles,fileFormat,unitConverter,header,time,append);},"fileName"_a,"particles"_a,"fileFormat"_a="Text",
        "unitConverter"_a=units::ExternalUnits(),"header"_a="","time"_a=-((float)((float)((1E300)*(1E300)*(0.0F)))),"append"_a=false);
     m.def("writeSnapshot",[](const std::string &fileName, const particles::ParticleArray<coord::PosVelCyl> &particles,
        const std::string &fileFormat="Text",
        const units::ExternalUnits &unitConverter=units::ExternalUnits(),const std::string &header="", 
        const double time=-((float)((float)((1E300)*(1E300)*(0.0F)))), const bool append=false)
        {particles::writeSnapshot(fileName,particles,fileFormat,unitConverter,header,time,append);},"fileName"_a,"particles"_a,"fileFormat"_a="Text",
        "unitConverter"_a=units::ExternalUnits(),"header"_a="","time"_a=-((float)((float)((1E300)*(1E300)*(0.0F)))),"append"_a=false);
    m.def("writeSnapshot",[](const std::string &fileName, const particles::ParticleArray<coord::PosVelSph> &particles,
        const std::string &fileFormat="Text",
        const units::ExternalUnits &unitConverter=units::ExternalUnits(),const std::string &header="", 
        const double time=-((float)((float)((1E300)*(1E300)*(0.0F)))), const bool append=false)
        {particles::writeSnapshot(fileName,particles,fileFormat,unitConverter,header,time,append);},"fileName"_a,"particles"_a,"fileFormat"_a="Text",
        "unitConverter"_a=units::ExternalUnits(),"header"_a="","time"_a=-((float)((float)((1E300)*(1E300)*(0.0F)))),"append"_a=false);
    m.def("updateTotalPotential",&galaxymodel::updateTotalPotential);
    m.def("doIteration",&galaxymodel::doIteration);
    m.def("createCylSpline",[](const potential::BaseDensity &src, int mmax,unsigned int gridSizeR,double Rmin,double Rmax,unsigned int gridSizez,
    double zmin,double zmax,bool useDerivs=true){return potential::CylSpline::create(src,mmax,gridSizeR,Rmin,Rmax,gridSizez,zmin,zmax,useDerivs);},
    "src"_a,"mmax"_a,"gridSizeR"_a,"Rmin"_a,"Rmax"_a,"gridSizez"_a,"zmin"_a,"zmax"_a,"useDerivs"_a=true);
    m.def("assignVelocity",&galaxymodel::assignVelocity);
    m.def("computeMoments",[](const galaxymodel::GalaxyModel &model,const coord::PosCyl &pos,bool Dens=true, 
        bool freqs=false, bool vel=false,bool vel2=false,const bool seperate=false, const double reqRelerror=0.001,
        const int maxNumEval=100000)->py::object{
        if(!Dens&&!freqs&&!vel&&!vel2)Dens=true;
        int dflen=seperate?model.distrFunc.numValues():1;
        if(dflen>1){
            std::vector<double> density(dflen),vels(dflen);
            std::vector<coord::Vel2Cyl> vel2s(dflen);
            std::vector<actions::Frequencies> freq(dflen);
            galaxymodel::computeMoments(model,pos,Dens?&density[0]:NULL,vel?&vels[0]:NULL,
                vel2?&vel2s[0]:NULL,freqs?&freq[0]:NULL,NULL,NULL,NULL,seperate,reqRelerror,maxNumEval);
            if((Dens&&vel)||(Dens&&vel2)||(Dens&&freqs)||(freqs&&vel)||(freqs&&vel2)||(vel&&vel2)){
                py::list ls;
                if(Dens)ls.append(density);
                if(vel)ls.append(vels);
                if(vel2)ls.append(vel2s);
                if(freqs)ls.append(freq);
                return ls;
            }
            if(Dens) return py::cast(density);
            if(vel) return py::cast(vels);
            if(vel2) return py::cast(vel2s);
            if(freqs) return py::cast(freq);
        }
        double density,vels;
        coord::Vel2Cyl vel2s;
        actions::Frequencies freq;
        galaxymodel::computeMoments(model,pos,Dens?&density:NULL,vel?&vels:NULL,
                vel2?&vel2s:NULL,freqs?&freq:NULL,NULL,NULL,NULL,seperate,reqRelerror,maxNumEval);
        if((Dens&&vel)||(Dens&&vel2)||(Dens&&freqs)||(freqs&&vel)||(freqs&&vel2)||(vel&&vel2)){
            py::list ls;
            if(Dens)ls.append(density);
            if(vel)ls.append(vels);
            if(vel2)ls.append(vel2s);
            if(freqs)ls.append(freq);
             return ls;
        }
        if(Dens) return py::cast(density);
        if(vel) return py::cast(vels);
        if(vel2) return py::cast(vel2s);
        return py::cast(freq);

    },"model"_a,"pos"_a,"Dens"_a=true,"freqs"_a=false,"vel"_a=false,"vel2"_a=false,"seperate"_a=false,
    "reqRelError"_a=0.001,"maxNumEval"_a=100000);

    m.def("computeProjectedMoments",[](const galaxymodel::GalaxyModel &model,const double R,
        const bool seperate=false, const double reqRelerror=0.001, const int maxNumEval=100000){
        py::list ls;
        int dflen=seperate?model.distrFunc.numValues():1;
        if(dflen>1){
            std::vector<double> density(dflen),rmsH(dflen),rmsV(dflen);
            galaxymodel::computeProjectedMoments(model,R,&density[0],&rmsH[0],&rmsV[0]
                ,NULL, NULL,NULL,seperate,reqRelerror,maxNumEval);
            ls.append(density);
            ls.append(rmsH);
            ls.append(rmsV);
        }else{
            double density,rmsH,rmsV;
            galaxymodel::computeProjectedMoments(model,R,&density,&rmsH,&rmsV
                    ,NULL, NULL,NULL,seperate,reqRelerror,maxNumEval);
            ls.append(density);
            ls.append(rmsH);
            ls.append(rmsV);
        }
        return ls;

    },"model"_a,"R"_a,"seperate"_a=false,"reqRelError"_a=0.001,"maxNumEval"_a=100000);
    m.def("computeVelocityDistributionO3",[](const galaxymodel::GalaxyModel &model,const coord::PosCyl &point,
        const bool projected, const std::vector<double>& gridVR,const std::vector<double>& gridVphi,
        const std::vector<double>& gridVz,const bool seperate=false,const double reqRelError=(0.01),
    const int maxNumEval=1000000){
        py::list ls;
        int dflen=seperate?model.distrFunc.numValues():1;
        std::vector<double> density,amplVR,amplVz,amplVphi;
        galaxymodel::computeVelocityDistributionO3(model,point,projected,gridVR,gridVz,gridVphi,&density[0],
            &amplVR,&amplVz,&amplVphi,seperate,reqRelError,maxNumEval);
        ls.append(density);
        ls.append(amplVR);
        ls.append(amplVphi); 
        ls.append(amplVz);
        return ls;

    },"model"_a,"point"_a,"projected"_a,"gridVR"_a,"gridVphi"_a,"gridVz"_a,"seperate"_a=false,"reqRelError"_a=(0.01),"maxNumEval"_a=1000000);
}