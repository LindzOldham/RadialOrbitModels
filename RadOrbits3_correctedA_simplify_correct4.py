import numpy as np, pylab as pl, pyfits as py
from scipy.special import gamma as Gamma
from scipy.interpolate import splrep, splev, splint
from scipy.optimize import brentq
from scipy import integrate

## a second calculation of the surface density - Wyn's way. This should be correct now.
#N0_num = np.pi**1.5 * Gamma((2.*gamma + alpha - 4.)/(2.*alpha)) * psi0**((gamma-2.)/alpha)
#N0_den = np.sqrt(2) * r0**2 * Gamma((gamma + alpha - 2.)/alpha)
#N0 = N0_num / N0_den

#S0_num = np.pi**2. * Gamma(0.5*(gamma-1.)) * Gamma((2.*gamma + alpha - 4.)/(2.*alpha)) * psi0**((gamma-2.)/alpha)
#S0_den = np.sqrt(2) * r0 * Gamma(0.5*gamma) * Gamma((gamma+alpha-2)/alpha)
#S0 = S0_num / S0_den

#SB = S0 * (r0/r)**gamma


def Proper(alpha,gamma,psi0):
    p = (2.-gamma)/alpha + 0.5 ############## corrected
    r0 = 1.
    r = np.logspace(1,8,1000)
    S0_num = np.pi**2. * Gamma(0.5*(gamma-1.)) * Gamma((2.*gamma + alpha - 4.)/(2.*alpha)) * psi0**((gamma-2.)/alpha)
    S0_den = np.sqrt(2) * r0 * Gamma(0.5*gamma) * Gamma((gamma+alpha-2)/alpha)
    S0 = S0_num / S0_den

    SB = S0 * (r0/r)**gamma
    SBmodel = splrep(r,SB)
    SB = splev(R,SBmodel)
    Ps=np.zeros(R.size)
    for l in range(R.size):
        r = np.logspace(1,8,1500)
        tosolve = lambda r: (r**2. - R[l]**2.) * psi0 * r0**alpha - 0.5 * vlos[l]**2. * r**(alpha+2.)
        if np.all(tosolve(r)<0):
            return -np.inf, p
        ## trial
        tosolve1 = (r**2. - R[l]**2.) * psi0 * r0**alpha - 0.5 * vlos[l]**2. * r**(alpha+2.)
        tosolve1,r1 = tosolve1[r<1000],r[r<1000]
        sort = tosolve1.argsort()
        solvermod = splrep(tosolve1[sort],r1[sort])
        rmintry = splev(0,solvermod)
        if tosolve(rmintry+150)<0:
            ### do it again!
            tosolve1,r1 = tosolve1[r<700],r[r<700]
            sort = tosolve1.argsort()
            solvermod = splrep(tosolve1[sort],r1[sort])
            rmintry = splev(0,solvermod)
        for xx in [500,400]:
            if tosolve(rmintry+100)<0:
                ### do it again!
                tosolve1,r1 = tosolve1[r<xx],r[r<xx]
                sort = tosolve1.argsort()
                solvermod = splrep(tosolve1[sort],r1[sort])
                rmintry = splev(0,solvermod)
                break
        if rmintry<0:
            tosolve1,r1 = tosolve1[np.where((r>100) & (r<700))],r[np.where((r>100) & (r<700))]
            sort = tosolve1.argsort()
            solvermod = splrep(tosolve1[sort],r1[sort])
            rmintry = splev(0,solvermod)
        if rmintry < 0:
            rmintry = 500
        if tosolve(rmintry+150)<0 and tosolve(rmintry+100)<0:
            #print rmintry, tosolve(rmintry+10),tosolve(rmintry+100)
            rmin = brentq(tosolve,10,rmintry+10)
        elif tosolve(rmintry+150)<0:
            rmin = brentq(tosolve,10,rmintry+100)
        else:
            rmin = brentq(tosolve,10,rmintry+150)
        
        if tosolve(1e8)<0:
            rmax = brentq(tosolve,rmin+10,1e8)
        
        else:
            for xx in range(10,80,2):
                if tosolve(10**xx)<0:
                    rmax = brentq(tosolve,10**(xx-2), 10**(xx))
                    break
        integral = r**(2.-gamma+0.5*alpha) / (tosolve(r)**p * (r**2- R[l]**2.)**(1.-p))   ################## corrected here
        if len(r[integral>0])<2:
            #print 'integral is full of nans: abort!!!'
            return -np.inf,p
        if rmin>1e6 or rmin<0:
            # bad model!
            return -np.inf,p
        ### proper integration
        theta = np.linspace(0,np.pi/2.,1000)
        theta[-1] = np.pi/2.
        ct = np.cos(theta)
        ct[-1] = 0
        c2,s2 = ct**2., np.sin(theta)**2.
        rint = rmin*c2 + rmax*s2
        tosolve = (rint**2. - R[l]**2.) * psi0 * r0**alpha - 0.5 * vlos[l]**2. * rint**(alpha+2.)
        integral =  2.*(rmax-rmin)*ct * np.sin(theta) * rint**(2.-gamma + 0.5*alpha) / (tosolve**p * (rint**2- R[l]**2.)**(1.-p)) ################ corrected here
        integral = np.nan_to_num(integral)
        import scipy
        P = integrate.cumtrapz(integral,theta)[-1]
        #print P
        Ps[l] = P * 2. * np.pi / SB[l]
    lnL = np.sum(np.log(Ps))
    return lnL, p


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
print R.size

# now let's run Proper for a grid of parameters!
#alphas = np.linspace(0.2,2.0,20)
#gammas = np.linspace(1.5,3.0,20)
#psi0s = np.logspace(7.0,10,20)
r0 = 1.
alphas = np.linspace(0.2,2.0,20)
gammas = np.linspace(1.5,3.0,20)
psi0s = np.logspace(7.0,13.5,30)

# now we have done a broad search ,we are gonna fo it again over a narrower range.
# we find alpha = 1.05, gamma = 2.37, logPsi0 = 9.24, sooo...
alphas = np.linspace(0.8,1.5,20)
gammas = np.linspace(2.0,2.6,20)
psi0s = np.logspace(8.8,9.8,30)

# now we find alpha = 1.13, gamma = 2.41, psi0 = 9.48
#alphas = np.linspace(0.9,1.2,20)
#gammas = np.linspace(2.35,2.5,20)
#psi0s = np.logspace(9.2,9.6,30)

#print np.log10(psi0s)
#alphas = np.linspace(1.0,1.6,25)
#gammas = np.linspace(2.0,3.0,25)
#psi0s = np.logspace(9.0,10.5,25)
#print np.log10(psi0s)


lnLs = np.zeros((alphas.size,gammas.size,psi0s.size))
ps = np.zeros((alphas.size,gammas.size,psi0s.size))
from itertools import product
for i,j,k in product(range(alphas.size),range(gammas.size),range(psi0s.size)):
    alpha,gamma,psi0 = alphas[i],gammas[j],psi0s[k]
    lnL,p = Proper(alpha,gamma,psi0)
    lnLs[i,j,k] = lnL
    ps[i,j,k] = p

''' plot results'''
#for i in range(10):
#    pl.figure()
#    pl.imshow(lnLs[i],interpolation='nearest')
#    pl.colorbar()
''' not entirely clear why this should be the case that there remain some nans - but there are only twoentries atm where this is the case. '''
#lnLs[np.where(np.isnan(lnLs)==True)]=-np.inf

lnLf = lnLs.copy()
lnLs[np.where(np.isnan(lnLs)==True)]=-np.inf

ai,gi,pi = np.where(lnLs == np.amax(lnLs))
alpha,gamma,psi0 = alphas[ai],gammas[gi],psi0s[pi]
print alpha,gamma,np.log10(psi0)

''' further restrict to where p<1 '''
ii = np.where(ps<1)
ai,gi,pi = np.where(lnLs == np.amax(lnLs[ii]))
alpha,gamma,psi0 = alphas[ai],gammas[gi],psi0s[pi]
print alpha,gamma,np.log10(psi0)

## make an interpolation object to get alpha, gamma, psi0
## don't know how to do this really

# plot for Wynevere
#print lnLs.shape
#for i in range(lnLs.shape[2]):
#    pl.figure()
#    pl.imshow(lnLs[:,:,i],interpolation='nearest',origin='lower',extent=[alphas[0],alphas[-1],gammas[0],gammas[-1]],aspect='auto')
#    pl.colorbar()
#    pl.xlabel(r'$\alpha$')
#    pl.ylabel(r'$\gamma$')
#for i in range(lnLs.shape[0]):
#    pl.figure()
#    pl.imshow(lnLs[i,:,:],interpolation='nearest',origin='lower',extent=[gammas[0],gammas[-1],psi0s[0],psi0s[-1]],aspect='auto')
#    pl.colorbar()
#    pl.xlabel(r'$\gamma$')
#    pl.ylabel(r'$\Psi_0$')

# compare to Sigma(R) and fit for F0
from DensityCalc import SurfaceDensity
Rbins,n,nerr = SurfaceDensity(R,9)

### CORRECT PSI0 FOR THE FACTOR OF FOUR
psi0 = psi0/4.


### model
p = (2.-gamma)/alpha + 0.5 ############## corrected
r0 = 1.
r = np.logspace(1,8,1000)
S0_num = np.pi**2. * Gamma(0.5*(gamma-1.)) * Gamma((2.*gamma + alpha - 4.)/(2.*alpha)) * psi0**((gamma-2.)/alpha)
S0_den = np.sqrt(2) * r0 * Gamma(0.5*gamma) * Gamma((gamma+alpha-2)/alpha)
S0 = S0_num / S0_den
SB = S0 * (r0/r)**gamma
SBmodel = splrep(r,SB)
SB = splev(Rbins,SBmodel)


# now fit F0
Rbins,n,nerr = SurfaceDensity(R,9)
pl.figure()
pl.plot(Rbins,n,'b-')
pl.errorbar(Rbins,n,yerr=nerr,color='CornflowerBlue',marker='o')
SB = splev(Rbins,SBmodel)
F0s = np.linspace(-4,-2,50)
lnL = np.zeros(50)
for i in range(50):
    F0 = 10**F0s[i]
    model = SB*F0
    chi2 = -0.5*(model-n)**2.
    lnL[i] = np.sum(chi2)

# interpolate: 
lnLmod = splrep(F0s,lnL)
F0s2 = np.linspace(-4,-1.5,1000)
lnLs = splev(F0s2,lnLmod)
nn = np.argmax(lnLs)
print F0s2[nn]

pl.figure()
pl.plot(Rbins,n,color='CornflowerBlue',label='data')
pl.errorbar(Rbins,n,yerr=nerr,color='CornflowerBlue',marker='o')
pl.plot(Rbins,SB*10**F0s2[nn],color='Black',label='power-law fit')
pl.text(600,0.00008,r'$\alpha = $'+'%.2f'%alpha.item())
pl.text(600,0.000075,r'$\gamma = $'+'%.2f'%gamma.item())
pl.text(600,0.000085,r'$\Psi_0 = $'+'%.2f'%np.log10(psi0.item()))
pl.text(600,0.00009,'$F_0 = $'+'%.2f'%F0s2[nn])
pl.xlabel('R/kpc')
pl.ylabel('SB /kpc$^{-2}$')
pl.legend(loc='lower left')
pl.title('surface density')

'''
m = np.argmax(lnL)
print F0s[m]
pl.figure()
pl.plot(Rbins,n,color='CornflowerBlue',label='data')
pl.errorbar(Rbins,n,yerr=nerr,color='CornflowerBlue',marker='o')
pl.plot(Rbins,SB*10**F0s[m],color='Black',label='power-law fit')
pl.text(600,0.00008,r'$\alpha = $'+str(alpha.item()))
pl.text(600,0.000075,r'$\gamma = $'+str(gamma.item()))
pl.text(600,0.000085,r'$\Psi_0 = $'+str(np.log10(psi0.item())))
pl.text(600,0.00009,'$F_0 = $'+str(F0s[m]))
pl.xlabel('R/kpc')
pl.ylabel('SB /kpc$^{-2}$')
pl.legend(loc='lower left')
pl.title('surface density')
'''
# also velocity dispersion.
def func(x):
    return np.std(x)

sr = np.linspace(100,1000,100)
v02 = psi0 * Gamma(gamma/2.) * Gamma(0.5*(alpha + gamma - 1.))
v02 /= 2. * Gamma(0.5*(gamma-1.)) * Gamma(0.5*(alpha + gamma + 2.))
vlos2 = v02 * (r0/sr)**alpha
pl.figure()
pl.plot(sr,vlos2**0.5,color='Black')

from scipy.stats import binned_statistic
means,bins,ns = binned_statistic(R,vlos,statistic=func,bins=10)
#pl.figure()
pl.plot(bins[1:],means,'o',color='CornflowerBlue')
pl.xlabel('R/kpc')
pl.ylabel(r'$\sigma_{los}^2$')
pl.title('velocity dispersion')
pl.legend(loc='upper right')

## also plot likelihood surface
