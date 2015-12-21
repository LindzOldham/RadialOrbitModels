import numpy as np, pylab as pl, pyfits as py
from scipy.special import gamma as Gamma
from scipy.interpolate import splrep, splev, splint
from scipy.optimize import brentq
from scipy import integrate
import GetSigmaML
import cPickle

### COMPARING CIRCULAR VELOCITIES PREDICTED FROM THE POWER-LAW POTENTIAL MODEL WITH THOSE FROM THE JEANS ANALYSIS - USING THE ISOTROPIC MODELS AS THE RADIAL-ORBIT MODEL IS ALSO ISOTROPIC

r = np.linspace(1,1000,500)
#########
### LOG MODEL
rho0 = 10**7.77
r0 = 48
G = 4.3e-6
vs2 = rho0 * r0**2. * G
vs = vs2**0.5
v_log = vs*r/(r**2. + r0**2.)**0.5

vs_lo, vs_hi = (10**7.67)*43**2. * G, (10**7.87)*53**2. *G
v_log_lo, v_log_hi = vs_lo**0.5 *r/(r**2. + 43**2.)**0.5, vs_hi**0.5 *r/(r**2. + 53**2.)**0.5

pl.figure()
pl.plot(r,v_log,color='CornflowerBlue',label='coreISO')
pl.fill_between(r,v_log,v_log_lo,color='LightCyan')#,alpha=0.5)
pl.fill_between(r,v_log_hi,v_log,color='LightCyan')#,alpha=0.5)
pl.ylim([0,1400])

######
### POWER LAW MODEL
#psi0 = 10**9.24
#alpha = 1.05
psi0 = 10**8.89
alpha = 1.13

r0 = 1.
v_pow = np.sqrt(alpha*psi0) * (r/r0)**(-alpha/2.)
v_lo, v_hi = np.sqrt(1.05*10**8.86) * (r/r0)**(-1.05/2.), np.sqrt(1.2*10**8.95) * (r/r0)**(-1.2/2.)

#pl.figure()
#pl.plot(r,v_pow)
#pl.figure()
pl.plot(r,v_pow,color='Crimson',label='power-law')
pl.fill_between(r,v_pow,v_lo, color='LavenderBlush')
pl.fill_between(r,v_pow,v_hi, color='LavenderBlush')

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
pl.legend(loc='lower right')
pl.xlabel('r / kpc')
pl.ylabel('$v_{c}$ / kms$^{-1}$')
