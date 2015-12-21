import numpy as np, pylab as pl, pyfits as py

# data
ra,dec,cz,dist2 = np.loadtxt('/data/ljo31/CFHT/Catalogues/EDD_EVCC_clipped_20.cat',unpack=True,usecols=(0,1,2,3))
vlos = cz - 1284. 
# calculate 3D radial separations from M87
ra_m87, dec_m87 =  187.7059304, 12.3911231
ra1,dec1 = (ra-ra_m87)*np.cos(0.5*(dec+dec_m87)), dec-dec_m87
R = np.sqrt(ra1**2 + dec1**2)*0.083*3600 # in kpc
# R<800
ii = np.where((R<800) & (R>1))
R,vlos = R[ii],vlos[ii]
ra,dec[ii]=ra[ii],dec[ii]
ii = np.where(R>200)
R,vlos=R[ii],vlos[ii]
ra,dec=ra[ii],dec[ii]

pl.figure()
pl.scatter(ra,dec,c=vlos,s=100,vmin=-1000,vmax=1000,cmap='Blues_r')
cbar = pl.colorbar()
cbar.set_label('$v_{los}$ / kms$^{-1}$', rotation=270)
pl.axis([185.3,191,9.5,16.5])
pl.scatter(ra_m87,dec_m87,s=200,marker='*',color='k',label='M87')
#pl.legend(loc='upper left')
pl.xlabel('ra / deg')
pl.ylabel('dec / deg')
