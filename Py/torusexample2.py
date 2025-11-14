import Py_agama as AG
import matplotlib.pyplot as plt


p=AG.createPotential("type=spheroid, gamma=1, beta=3, alpha=1, scaleradius=1.0, densityNorm=0.1, q=0.5")
#p.eval([x,y,z], bool pot=False, bool acc=False, bool der=False). Where if acc is true returns acceleration 
# (minus gradient). If deriv true returns hessian of potential and if pot true returns potential.
#e.g. acc,deriv=p.eval([1.1,2.3,0.3],acc=True,deriv=True) will get acc as 3d array containing acceleration, 
# der 6d array containing hessian
#Torus and TorusGenerator class. Torus(Potential p, tol=1e-4). tolerance set to 1e-4 but can change
TG=AG.TorusGenerator(p,1e-4)
J=AG.Actions(0.1,0.3,3.0)
theta0=AG.Angles(1.2,0.3,0.4)
T=TG.fitTorus(J)
AF=AG.ActionFinderAxisymFudge(p)
AFTG=AG.ActionFinderTG(p,AF,TG)
xp0=T.from_true(theta0)
xv0=AG.PosVelCyl(xp0.R,xp0.z,xp0.phi,xp0.pR,xp0.pz,xp0.pphi/xp0.R)
print(f"R:{xp0.R} z:{xp0.z} phi:{xp0.phi} pR:{xp0.pR} pz:{xp0.pz} pphi:{xp0.pphi}")
f=AG.Omegas(1,2,3)
aa1=AF.actionAngles(xv0)
print(f"AF Jr:{aa1.Jr} Jz:{aa1.Jz} thetar:{aa1.thetar} thetaz:{aa1.thetaz} thetaphi:{aa1.thetaphi}")
aa2=AFTG.actionAngles(xv0)
print(f"AFTG Jr:{aa2.Jr} Jz:{aa2.Jz} thetar:{aa2.thetar} thetaz:{aa2.thetaz} thetaphi:{aa2.thetaphi}")
R=[]
z=[]
time=[]
R2=[]
z2=[]
A=T.orbit(theta0,1,400)
print(len(A))
for i in range(len(A)):
    R.append(A[i][0].R)
    z.append(A[i][0].z)
    time.append(A[i][1])
A2=AG.integrateTraj(xv0,400,1,p)
for i in range(len(A2)):
    R2.append(A2[i][0].R)
    z2.append(A2[i][0].z)
#numerical integration from intiial position and velocity to get phase space point along trajectory
plt.plot(R,z)
plt.plot(R2,z2)
plt.show()