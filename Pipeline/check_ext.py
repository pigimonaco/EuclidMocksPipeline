import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
from astropy.io import fits
import sys

# check of the application of extinction to a mock

Nside=256

z1 = 1.1
z2 = 1.2

print("Reading catalog...")
cat = fits.getdata('../Products/GalaxyCatalogs/hodcat_8614_100sqdeg_m3_cMdiemer19.fits')

mysel = fits.getdata('../Products/Selections/sel_MWext_hodcat_8614_100sqdeg_m3_cMdiemer19.fits')['SELECTION']


print("Filtering catalog...")
filter = (cat['true_redshift_gal']>=z1) & (cat['true_redshift_gal']<z2)

print("Finding pixels for {} galaxies...".format(filter.sum()))
conv = np.pi/180.
pix  = hp.ang2pix(Nside, np.pi/2. - cat['dec_gal'][filter]*conv, cat['ra_gal'][filter]*conv)
map1 = np.zeros((hp.nside2npix(Nside)), dtype=np.float)
for p in pix:
    map1[p]+=1

pix  = hp.ang2pix(Nside, np.pi/2. - cat['dec_gal'][filter & mysel]*conv, cat['ra_gal'][filter & mysel]*conv)
map2 = np.zeros((hp.nside2npix(Nside)), dtype=np.float)
for p in pix:
    map2[p]+=1

fp = fits.getdata('../Products/Footprints/100sqdeg.fits')['FOOTPRINT_G']
fp = hp.ud_grade(fp,Nside)>0.5

mw = hp.read_map('../Products/ExtinctionMaps/HFI_CompMap_ThermalDustModel_2048_R1.20.fits',field=2)

mw = hp.ud_grade(mw,Nside)
mw[~fp]=hp.UNSEEN
map1[~fp]=hp.UNSEEN
map2[~fp]=hp.UNSEEN

ff=map1==0
rat = map2/map1
rat[ff]=0
rat[~fp] = hp.UNSEEN

hp.mollview(map1,rot=[0,0,0],title='no mask')
hp.mollview(map2,rot=[0,0,0],title='mask')
hp.mollview(rat,rot=[0,0,0],title='ratio')
hp.mollview(mw,rot=[0,0,0],title='mw',max=0.2)

plt.figure()
plt.scatter(map1[fp], map2[fp], marker='.', s=1)

plt.figure()
plt.scatter(rat[fp], mw[fp], marker='.', s=1)

plt.show()
