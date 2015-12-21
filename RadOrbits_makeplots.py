import numpy as np, pylab as pl, pyfits as py
from scipy.special import gamma as Gamma
from scipy.interpolate import splrep, splev, splint
from scipy.optimize import brentq
from scipy import integrate

# model
alpha,gamma,psi0,F0 = 1.13,2.41, 10**8.89, 10**-2.59


# data
ra,dec,cz,dist2 = np.loadtxt('/data/ljo31/CFHT/Catalogues/EDD_EVCC_clipped_20.cat',unpack=True,usecols=(0,1,2,3))
vlos = cz - 1284. 
# calculate 3D radial separations from M87
ra_m87, dec_m87 =  187.7059304, 12.3911231
ra,dec = (ra-ra_m87)*np.cos(0.5*(dec+dec_m87)), dec-dec_m87
R = np.sqrt(ra**2 + dec**2)*0.083*3600 # in kpc
# R<800
ii = np.where((R<800) & (R>1))
R,vlos = R[ii],vlos[ii]
ii = np.where(R>200)
R,vlos=R[ii],vlos[ii]

from DensityCalc import SurfaceDensity
#Rbins,n,nerr = SurfaceDensity(R,11)
Rbins,n,nerr = SurfaceDensity(R,8)

### model
p = (2.-gamma)/alpha + 0.5 ############## corrected
r0 = 1.
r = np.logspace(1,8,1000)
S0_num = np.pi**2. * Gamma(0.5*(gamma-1.)) * Gamma((2.*gamma + alpha - 4.)/(2.*alpha)) * psi0**((gamma-2.)/alpha)
S0_den = np.sqrt(2) * r0 * Gamma(0.5*gamma) * Gamma((gamma+alpha-2)/alpha)
S0 = S0_num / S0_den
SB = F0 * S0 * (r0/r)**gamma
SBmodel = splrep(r,SB)
SB = splev(Rbins,SBmodel)

def CompSB(alpha,gamma,psi0):
    S0_num = np.pi**2. * Gamma(0.5*(gamma-1.)) * Gamma((2.*gamma + alpha - 4.)/(2.*alpha)) * psi0**((gamma-2.)/alpha)
    S0_den = np.sqrt(2) * r0 * Gamma(0.5*gamma) * Gamma((gamma+alpha-2)/alpha)
    S0 = S0_num / S0_den
    SB = F0 * S0 * (r0/r)**gamma
    #SBmodel = splrep(r,SB)
    #SB = splev(Rbins,SBmodel)
    return SB

SB_low, SB_hi = CompSB(1.05,2.3,10**8.86), CompSB(1.2,2.5,10**8.95)
SB = CompSB(1.13,2.41,10**8.89)
print r.size, SB.size

pl.figure()
pl.plot(r,1e4*SB,color='Black',label='power-law fit')
pl.fill_between(r,1e4*SB_low,1e4*SB,color='Gainsboro')
pl.fill_between(r,SB_hi*1e4,SB*1e4,color='Gainsboro')
pl.plot(Rbins,n*1e4,'o',color='CornflowerBlue',label='data')
pl.errorbar(Rbins,n*1e4,yerr=1e4*nerr,color='CornflowerBlue',fmt='o')
pl.xlabel('R/kpc')
pl.ylabel('SB /kpc$^{-2}$')
pl.legend(loc='upper right')
pl.xlim([200,850])
pl.ylim([0,1])

# also velocity dispersion.
def func(x):
    if len(x)==0:
        return 0.
    return np.std(x)#, np.std(x)/np.sqrt(len(x))

sr = np.linspace(100,1000,100)
v02 = psi0 * Gamma(gamma/2.) * Gamma(0.5*(alpha + gamma - 1.))
v02 /= 2. * Gamma(0.5*(gamma-1.)) * Gamma(0.5*(alpha + gamma + 2.))
vlos2 = v02 * (r0/sr)**alpha

def CompVlos(alpha,gamma,psi0):
    v02 = psi0 * Gamma(gamma/2.) * Gamma(0.5*(alpha + gamma - 1.))
    v02 /= 2. * Gamma(0.5*(gamma-1.)) * Gamma(0.5*(alpha + gamma + 2.))
    vlos2 = v02 * (r0/sr)**alpha
    return vlos2

vlos2_lo, vlos2_hi = CompVlos(1.05,2.3,10**8.95), CompVlos(1.2,2.5,10**8.86)

pl.figure()
pl.plot(sr,vlos2**0.5,color='Black',label='power-law fit')
pl.fill_between(sr,vlos2_lo**0.5, vlos2**0.5,color='Gainsboro')
pl.fill_between(sr,vlos2_hi**0.5, vlos2**0.5,color='Gainsboro')



from scipy.stats import binned_statistic
#means,bins,ns = binned_statistic(R,vlos,statistic=func,bins=10)
means,bins,ns = binned_statistic(R,vlos,statistic=func,bins=7)

# assign bin numbers
numbers = np.zeros(len(means))
for i in range(len(numbers)):
    ii = np.where((R>bins[i]) & (R<bins[i+1]))
    if len(R[ii])==0:
        numbers[i] = 1.
    else:
        numbers[i] = len(R[ii])
print numbers
#pl.figure()
pl.plot(bins[1:],means,'o',color='CornflowerBlue',label='data')
pl.errorbar(bins[1:],means,yerr=means/2/np.sqrt(numbers),fmt='o',color='CornflowerBlue')
pl.xlabel('R/kpc')
pl.ylabel(r'$\sigma_{los}$ / kms$^{-2}$')
pl.legend(loc='upper right')
pl.ylim([0,900])
pl.xlim([200,850])

# switch between 7 and 10 bins (8 and 11)
