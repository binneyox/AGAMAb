import numpy as np
import agama as gam
import matplotlib.pyplot as plt
IU=gam.galactic_kms # Conversions to InternalUnits
sun=gam.solarShifter(IU)
h=IU.from_Kpc_kms
p=gam.createPotential("type=spheroid, gamma=1, beta=3, alpha=1, scaleradius=18, densityNorm=170, q=0.5")
print("Vcirc(8kpc) ",IU.to_kms*gam.Vcirc(p,8*IU.from_Kpc))
TG=gam.TorusGenerator(p,5e-5)
J=gam.Actions(50*h,800*h,2800*h)
theta0=gam.Angles(0,0,4)
T=TG.fitTorus(J)
s=0
vr=0
l=[]
b=[]
traj=T.orbit(theta0,1*IU.from_Myr,.2*IU.from_Gyr)
for i in range(len(traj)):
    astrom=sun.toSky(traj[i][0],s,vr)
    l.append(astrom.pos.l)
    b.append(astrom.pos.b)
plt.xlabel("longitude")
plt.ylabel("b")
plt.plot(l,b,color="b")
plt.show()