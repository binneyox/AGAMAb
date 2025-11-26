import numpy as np
import agama as AG
import matplotlib.pyplot as plt
IU=AG.galactic_kms # Conversions to InternalUnits
sun=AG.solarShifter(IU)
h=IU.from_Kpc_kms
p=AG.createPotential("type=spheroid, gamma=1, beta=3, alpha=1, scaleradius=18, densityNorm=170, q=0.5")
print("Vcirc(8kpc) ",IU.to_kms*AG.Vcirc(p,8*IU.from_Kpc))
TG=AG.TorusGenerator(p,5e-5)
J=AG.Actions(50*h,800*h,2800*h)
theta0=AG.Angles(0,0,4)
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
