import numpy as np, pyfits as py, cPickle
from scipy.interpolate import splrep, splev, splint
from itertools import product
import pylab as pl

ra,dec,R,vlos,rmaxes,rhalfs,Emaxes,Ehalfs = np.load('/data/ljo31/NWE/energies.npy').T

# for each object, construct D_KL
# start with object 0
#x = np.linspace(-1*10**6.1,2e6,15000)
x1 = np.logspace(4.,6.5,8000)
x = np.concatenate((x1[::-1]*-1., x1))
#E = Emaxes[0]
#dE2 = ((Emaxes[0]-Ehalfs[0])*50.)**2. # giving it a big uncertainty for now
#lp = -0.5*(x-E)**2./dE2 - 0.5*np.log(2.*np.pi*dE2)
#lq = -0.5*(x-Emaxes[1])**2./((Emaxes[1]-Ehalfs[1])*500)**2. - 0.5*np.log(2.*np.pi*(Emaxes[1]-Ehalfs[1])**2.*500**2.)
#pl.figure()
#pl.plot(x,np.exp(lp))
#pl.plot(x,np.exp(lq))

lpfunc = lambda x,E,dE2: -0.5*(x-E)**2./dE2 - 0.5*np.log(2.*np.pi*dE2)
lps = np.zeros((Emaxes.size,x.size))
for i in range(Emaxes.size):
    E,dE2 = Emaxes[i], 0.25*rmaxes[i]**2. * vlos[i]**2. / (rmaxes[i]**2. - R[i]**2.)
    dE2 = dE2**2.
    lps[i] = lpfunc(x,E,dE2)
    if E<0:
        j = np.where(E<x)
    elif E>0:
        j = np.where(E>x)
    lps[i,j] = -1e6


#pl.figure()
#pl.plot(x,np.exp(lps).T)

#np.save('/data/ljo31/NWE/lps_1D_v2',lps) # just energies



#lps = np.load('/data/ljo31/NWE/lps_1D_v2.npy')
    
Dkls = np.zeros((Emaxes.size,Emaxes.size))
for i,j in product(range(Emaxes.size),range(Emaxes.size)):
    lp,lq=lps[i],lps[j]
    func = np.exp(lp) * (lp-lq)
    model = splrep(x,func)
    Dkl = splint(x[0],x[-1],model)
    Dkls[i,j] = Dkl
    #print Dkl

pl.figure()
pl.imshow(Dkls,interpolation='nearest',origin='lower',vmax=0.5)#,cmap='Blues_r')
cbar = pl.colorbar()
cbar.set_label('$D_{kl}$', rotation=270,labelpad=30)

np.save('/data/ljo31/NWE/Dkls_1D_v2',Dkls)
pl.xlabel('galaxy index 1')
pl.ylabel('galaxy index 2')
#pl.title('coloured by $D_{kl}$ - Gaussian model')

positives = []
pairs = []
for i,j in product(range(45),range(45)):
    if i!=j:
        if Dkls[i,j] < 0.5:
            #print  i, j
            positives.append([Dkls[i,j],Emaxes[i],Emaxes[j],i,j])
            pairs.append([i,j])

pairs=np.array(pairs)

#np.save('/data/ljo31/NWE/matches_gaussian',positives)

# need to think about how to recognise groups.
# we could make groups for each one, then compare and see if they're the same.
u = np.unique(pairs[:,0])
groups = []
for i in u:                 
    ii = np.where(pairs[:,0] == i)
    groups.append(np.concatenate((np.array([i]),pairs[ii,1][0])))


clumps = []
for i in range(len(u)):
    cal = u[i]
    cals = []
    for j in range(len(groups)):
        g = groups[j]
        if cal in g:
            [cals.append(f) for f in g]
    clumps.append([cals])

clump2s = []
for i in range(len(clumps)):
    c = np.array(clumps[i])
    us = np.unique(c)
    for j in range(len(clump2s)):
        if us[0] in clump2s[j] or us[1] in clump2s[j]:
            #print us, clump2s[j]
            clump2s[j] = np.concatenate((clump2s[j],us))
            continue
        continue
    clump2s.append(us)

for p in range(len(clump2s)):
    clump2s[p] = np.unique(clump2s[p])

from collections import defaultdict

arrs = []
# but this seems a bit messy...
g1 = np.array([0,3,4,7,9,13,14,17,23,24,39,40])
g2 = np.array([1,6,12,18,19,21,27,29,30,32,33,36,38])
g3 = np.array([15,20])
#g4 = np.array([12,19,33])
g5 = np.array([2,8,28,31,35,37,42,43])
gl = np.array([5,10,11,16,22,25,26,34,41,44])
print len(g1),len(g2),len(g3),len(g5),len(gl)
pl.figure()
pl.scatter(rmaxes[gl],Emaxes[gl]*1e-6,s=100,color='k',label='ungrouped')
pl.scatter(rmaxes[g1],Emaxes[g1]*1e-6,s=100,color='CornflowerBlue',label='group 1')
pl.scatter(rmaxes[g2],Emaxes[g2]*1e-6,s=100,color='Crimson',label='group 2')
pl.scatter(rmaxes[g3],Emaxes[g3]*1e-6,s=100,color='Purple',label='group 3')
#pl.scatter(rmaxes[g4],Emaxes[g4]*1e-6,s=100,color='Orange',label='group 4')
pl.scatter(rmaxes[g5],Emaxes[g5]*1e-6,s=100,color='SeaGreen',label='group 4')
pl.legend(loc='lower right')
pl.xlabel('r / kpc')
#pl.ylabel('E/m / kms$^{-1}$')
pl.ylabel('E/m x $10^{-6}$/ km$^{2}$s$^{-2}$')
pl.xlim([100,2700])
pl.ylim([-1.2,0.6])


## also plot Gaussian distributions
lpfunc = lambda x,E,dE2: -0.5*(x-E)**2./dE2 - 0.5*np.log(2.*np.pi*dE2)
lps = np.zeros((Emaxes.size,x.size))
pl.figure()
#for i in gl:
#    E,dE2 = Emaxes[i], 0.25*rmaxes[i]**2. * vlos[i]**2. / (rmaxes[i]**2. - R[i]**2.)
#    dE2 = dE2**2.
#    gauss = np.exp(lpfunc(x,E,dE2))
#    if E<0:
#        j = np.where(E<x)
#    elif E>0:
#        j = np.where(E>x)
#    gauss[j] = 0
#    gauss /=np.sum(gauss)
#    pl.plot(x*1e-6,gauss,color='k')
for i in g1:
    E,dE2 = Emaxes[i], 0.25*rmaxes[i]**2. * vlos[i]**2. / (rmaxes[i]**2. - R[i]**2.)
    dE2 = dE2**2.
    gauss = np.exp(lpfunc(x,E,dE2))
    if E<0:
        j = np.where(E<x)
    elif E>0:
        j = np.where(E>x)
    gauss[j] = 0
    gauss /=np.sum(gauss)
    pl.plot(x*1e-6,gauss,color='CornflowerBlue')
for i in g2:
    E,dE2 = Emaxes[i], 0.25*rmaxes[i]**2. * vlos[i]**2. / (rmaxes[i]**2. - R[i]**2.)
    dE2 = dE2**2.
    gauss = np.exp(lpfunc(x,E,dE2))
    if E<0:
        j = np.where(E<x)
    elif E>0:
        j = np.where(E>x)
    gauss[j] = 0
    gauss /=np.sum(gauss)
    pl.plot(x*1e-6,gauss,color='Crimson')
for i in g3:
    E,dE2 = Emaxes[i], 0.25*rmaxes[i]**2. * vlos[i]**2. / (rmaxes[i]**2. - R[i]**2.)
    dE2 = dE2**2.
    gauss = np.exp(lpfunc(x,E,dE2))
    if E<0:
        j = np.where(E<x)
    elif E>0:
        j = np.where(E>x)
    gauss[j] = 0
    gauss /=np.sum(gauss)
    pl.plot(x*1e-6,gauss,color='Purple')
#for i in g4:
#    E,dE2 = Emaxes[i], 0.25*rmaxes[i]**2. * vlos[i]**2. / (rmaxes[i]**2. - R[i]**2.)
#    dE2 = dE2**2.
#    gauss = np.exp(lpfunc(x,E,dE2))
#    if E<0:
#        j = np.where(E<x)
#    elif E>0:
#        j = np.where(E>x)
#    gauss[j] = 0
 ##   gauss /=np.sum(gauss)
#    pl.plot(x*1e-6,gauss,color='Orange')    
for i in g5:
    E,dE2 = Emaxes[i], 0.25*rmaxes[i]**2. * vlos[i]**2. / (rmaxes[i]**2. - R[i]**2.)
    dE2 = dE2**2.
    gauss = np.exp(lpfunc(x,E,dE2))
    if E<0:
        j = np.where(E<x)
    elif E>0:
        j = np.where(E>x)
    gauss[j] = 0
    gauss /=np.sum(gauss)
    pl.plot(x*1e-6,gauss,color='SeaGreen') 
pl.xlabel('E/m x $10^{-6}$/ km$^2$s$^{-2}$')
pl.ylabel('p(E)')
pl.xlim([-1.5,1.5])
pl.ylim([0,0.012])


'''
for i in range(len(clump2s)):
    pl.figure()
    pl.scatter(rmaxes[clump2s[i]],Emaxes[clump2s[i]],s=100)

pl.figure()
pl.scatter(rmaxes[clump2s[0]],Emaxes[clump2s[0]],s=100,color='CornflowerBlue')
pl.scatter(rmaxes[clump2s[1]],Emaxes[clump2s[1]],s=100,color='Crimson')
pl.scatter(rmaxes[clump2s[2]],Emaxes[clump2s[2]],s=100,color='Purple')
pl.scatter(rmaxes[clump2s[3]],Emaxes[clump2s[3]],s=100,color='Orange')
pl.scatter(rmaxes[clump2s[4]],Emaxes[clump2s[4]],s=100,color='Gray')
pl.scatter(rmaxes[clump2s[5]],Emaxes[clump2s[5]],s=100,color='ForestGreen')
'''
