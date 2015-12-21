import numpy as np, pylab as pl, pyfits as py
from scipy.special import gamma as Gamma
from scipy.interpolate import splrep, splev, splint
from scipy.optimize import brentq
from scipy import integrate
import GetSigmaML
import cPickle

### COMPARING CIRCULAR VELOCITIES PREDICTED FROM THE POWER-LAW POTENTIAL MODEL WITH THOSE FROM THE JEANS ANALYSIS - USING THE ISOTROPIC MODELS AS THE RADIAL-ORBIT MODEL IS ALSO ISOTROPIC

r = np.linspace(100,1000,500)
#########
### LOG MODEL
rho0 = 10**7.33
r0 = 78
G = 4.3e-6
vs2 = rho0 * r0**2. * G
vs = vs2**0.5
v_log = vs*r/(r**2. + r0**2.)**0.5

pl.figure()
pl.plot(r,v_log,color='CornflowerBlue',label='LOG')



######
### POWER LAW MODEL
#psi0 = 10**9.24
#alpha = 1.05
psi0 = 10**9.49
alpha = 1.13

r0 = 1.
v_pow = np.sqrt(alpha*psi0) * (r/r0)**(-alpha/2.)

#pl.figure()
#pl.plot(r,v_pow)
#pl.figure()
pl.plot(r,v_pow/2,color='Crimson',label='power-law')

######
## gNFW MODEL - use mass profile
## save mass profile so I don't have to KEEP generating it though!!!
#modchain = cPickle.load(open('/data/ljo31/CFHT/Dynamics/ModelsPSF/FINAL_CCC/fchain_MODEL2_anisotropy_AUGUST_andGCs.dat'))
#modlnprob = cPickle.load(open('/data/ljo31/CFHT/Dynamics/ModelsPSF/FINAL_CCC/flnprob_MODEL2_anisotropy_AUGUST_andGCs.dat'))
#ii = np.argmax(modlnprob)
#Mstar, log_rhodm, r0,gamma,noise,betastar,betared,betablue = modchain[ii]
log_rhodm, r0, gamma = 7.93, 121, 0
lp_args = [gamma,r0]
mod = GetSigmaML.BrokenPL
rhodm = 10**log_rhodm
Mdm = rhodm * mod(r,lp_args)
v_gnfw = (G*Mdm/r)**0.5

#pl.figure()
pl.plot(r,v_gnfw,color='Orange',label='gNFW')
pl.legend(loc='upper right')
pl.xlabel('r / kpc')
pl.ylabel('$v_{c}$ / kms$^{-1}$')


'''
######
## CHECK - LOG MODEL from mass
modchain = cPickle.load(open('/data/ljo31/CFHT/Dynamics/ModelsPSF/FINAL_CCC/fchain_MODEL3_anisotropy_AUGUST_andGCs.dat'))
modlnprob = cPickle.load(open('/data/ljo31/CFHT/Dynamics/ModelsPSF/FINAL_CCC/flnprob_MODEL3_anisotropy_AUGUST_andGCs.dat'))
ii = np.argmax(modlnprob)
Mstar, log_rhodm, rs,noise,betastar, betared, betablue = modchain[ii]
lp_args = rs
mod = GetSigmaML.Logged
rhodm = 10**log_rhodm
Mdm = rhodm * mod(r,lp_args)
v_log2 = (G*Mdm/r)**0.5

pl.figure()
pl.plot(r,v_log2)
'''



'''
psi0 = 10**9.24
alpha = 1.05
r0 = 1.
r = np.logspace(0,2.6,2000)
potential = psi0 * (r0/r)**alpha

### gNFW
G = 4.3e-6
a = 121
rho0 = 10**7.93
#gnfw = 1./(a+r) - np.log(r)/a + np.log(a+r)/a
#gnfw *= 0.5*rho0*a**2. * G

gnfw = 2.*np.log(1.+(r/a))/r - 1/(a+r) #+ 1
gnfw *= 0.5*G*rho0*a**3

nfw = np.log(1.+(r/a))/r
nfw*= G*rho0*a**3

### also plot drho by dr
gamma = 0
r0 = 121.
dpdr = (3.*r + gamma*r0) / (r+r0)

## well, the normalisation is wonky, but one thing we can do is to work out the slopes.
## slope of radial orbit model = alpha = 1.05

mod = splrep(r,gnfw,k=1)
slope = splev(r,mod,der=1)
# slope ~ -1 from measuring off the plot. So slopes match but normalisation is still CRAZY

## log model
rho0 = 10**7.33
r0 = 78
G = 4.3e-6
vs2 = rho0 * r0**2. * G
phi = 0.5*vs2 * np.log(r0**2. + r**2.)

pl.figure()
pl.loglog(r,phi)
'''
