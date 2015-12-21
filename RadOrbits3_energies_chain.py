import numpy as np, pylab as pl, pyfits as py
from scipy.special import gamma as Gamma
from scipy.interpolate import splrep, splev, splint
from scipy.optimize import brentq
from scipy import integrate
import GetSigmaML
import cPickle

### USING RADIAL ORBIT MODEL TO CALCULATE ENERGIES OF THE SATELLITES AND LOOK FOR ANY CLUSTERING


G = 4.3e-6

######
### POWER LAW MODEL
#psi0 = 10**9.24
#alpha = 1.05
#psi0 = 0.25*10**9.49
psi0 = 10**8.89
alpha = 1.13
r0 = 1.
#F0 = (10**-2.80)/0.6 # correcting for the factor of 4 in psi0
F0 = 10**-2.59
gamma = 2.41

# load up satellite data
ra,dec,cz,dist2 = np.loadtxt('/data/ljo31/CFHT/Catalogues/EDD_EVCC_clipped_20.cat',unpack=True,usecols=(0,1,2,3))
vlos = cz - 1284. 
# calculate 3D radial separations from M87
ra_m87, dec_m87 =  187.7059304, 12.3911231
ra,dec = (ra-ra_m87)*np.cos(0.5*(dec+dec_m87)), dec-dec_m87
R = np.sqrt(ra**2 + dec**2)*0.083*3600 # in kpc
# R<800
ii = np.where((R<800) & (R>1))
R,vlos = R[ii],vlos[ii]
ra,dec=ra[ii],dec[ii]
ii = np.where(R>200)
R,vlos=R[ii],vlos[ii]
ra,dec=ra[ii],dec[ii]

# construct density profile as a function of r
r = np.linspace(100,2800,5000)
N0 = np.pi**1.5 * Gamma((2.*gamma + alpha - 4.)/(2.*alpha)) / Gamma(gamma/2.) / Gamma((gamma+alpha-2.)/alpha) / (np.sqrt(2)*r0**2)
N0 *= F0 * psi0**((gamma-2.)/alpha)
nu = N0 * (r/r0)**-gamma

# and energy
psi = lambda r: psi0*(r/r0)**-alpha
energies = []
things = []
pl.figure()
Emaxes = []
rmaxes = []
Ehalfs = []
rhalfs = []
for i in range(R.size):
    #print R[i],vlos[i]
    E = lambda r: 0.5*r**2. * vlos[i]**2./(r**2 - R[i]**2.) - psi0*(r/r0)**-alpha
    rr = r[r>R[i]]
    energies.append(np.column_stack((rr,E(rr))))
    model = splrep(rr,E(rr))
    dEdr = splev(rr,model,der=1)
    prob = nu[r>R[i]]/dEdr
    # return E,r,prob
    things.append(np.column_stack((E(rr),prob)))
    #pl.plot(E(r),abs(prob)/np.sum(abs(prob)))
    #ii = np.where(prob>0)
    #pl.plot(E(rr),abs(prob)/np.sum(abs(prob)))
    #pl.plot(rr,abs(prob)/np.sum(abs(prob)))
    prob = abs(prob)/np.sum(abs(prob))
    Emaxind = np.argmax(prob)
    Emax = E(rr)[Emaxind]
    Emaxes.append(Emax)
    rmax = rr[Emaxind]
    rmaxes.append(rmax)
    #pl.figure()
    pl.plot(rr,prob)
    probmax=prob[Emaxind]
    pl.plot(rmax,probmax,'*')
    probhalf = 0.5*probmax
    rc = rr-rmax
    ii=np.where(rc>=0)
    mod = splrep(prob[ii][0:-1][::-1],rc[ii][0:-1][::-1],k=1)
    rhalf = splev(probhalf,mod)
    #p#l.plot(splev(prob[ii][1:-1],mod),prob[ii][1:-1])
    print i,rmax,rhalf,probmax,probhalf
    pl.plot(rmax+rhalf,probhalf,'*')
    pl.plot(rmax-rhalf,probhalf,'*') # generous limits
    Ehalf = E(rmax+rhalf)
    Ehalf2=E(rmax-rhalf)
    Ehalfs.append(Ehalf)
    rhalfs.append(rhalf)


Emaxes=np.array(Emaxes)
rmaxes=np.array(rmaxes)
Ehalfs=np.array(Ehalfs)
rhalfs=np.array(rhalfs)

pl.figure()
pl.scatter(R,vlos,c=Emaxes,edgecolors='none',s=100,cmap='hot')
pl.colorbar()
pl.xlabel('R')
pl.ylabel('$v_{los}$')
pl.title('coloured by energy')

pl.figure()
pl.hist(Emaxes,10)

pl.figure()
pl.scatter(rmaxes,Emaxes)


ra,dec,cz,dist2 = np.loadtxt('/data/ljo31/CFHT/Catalogues/EDD_EVCC_clipped_20.cat',unpack=True,usecols=(0,1,2,3))
vlos = cz - 1284. 
# calculate 3D radial separations from M87
ra_m87, dec_m87 =  187.7059304, 12.3911231
ra2,dec2 = (ra-ra_m87)*np.cos(0.5*(dec+dec_m87)), dec-dec_m87
R = np.sqrt(ra2**2 + dec2**2)*0.083*3600 # in kpc
# R<800
ii = np.where((R<800) & (R>1))
R,vlos = R[ii],vlos[ii]
ra,dec=ra[ii],dec[ii]
ii = np.where(R>200)
R,vlos=R[ii],vlos[ii]
ra,dec=ra[ii],dec[ii]

pl.figure()
pl.scatter(ra,dec,c=Emaxes,edgecolors='none',s=150)
pl.colorbar()
grid = np.column_stack((ra,dec,R,vlos,rmaxes,rhalfs,Emaxes,Ehalfs))
np.save('/data/ljo31/NWE/energies',grid)

# we need to also characterise them by some uncertainty

pl.figure()
pl.scatter(rmaxes,Emaxes,marker='o')
pl.errorbar(rmaxes,Emaxes,yerr=-1.*(Emaxes-Ehalfs),fmt='o')
