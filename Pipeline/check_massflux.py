import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

flag=fits.getdata('../Products/GalaxyCatalogs/flagship_8614_100sqdeg_m3.fits')

ffc=(flag['observed_redshift_gal']>=1.0) & (flag['observed_redshift_gal']<1.05) & (flag['kind']==0)
ffs=(flag['observed_redshift_gal']>=1.0) & (flag['observed_redshift_gal']<1.05) & (flag['kind']==1)

plt.figure()
plt.scatter(flag['halo_lm'][ffc],flag['logf_halpha_model3_ext'][ffc],marker='.',s=1)
plt.title('flagship centrals')
plt.ylim([-15.8,-14.6])
plt.xlabel('Log halo mass')
plt.ylabel('Log Halpha flux')

plt.figure()
plt.scatter(flag['halo_lm'][ffs],flag['logf_halpha_model3_ext'][ffs],marker='.',s=1)
plt.title('flagship satellites')
plt.ylim([-15.8,-14.6])
plt.xlabel('Log halo mass')
plt.ylabel('Log Halpha flux')

hod=fits.getdata('../Products/GalaxyCatalogs/hodcat_8614_100sqdeg_m3_cMdiemer19.fits')
hhc=(hod['observed_redshift_gal']>=1.0) & (hod['observed_redshift_gal']<1.05) & (hod['kind']==0)
hhs=(hod['observed_redshift_gal']>=1.0) & (hod['observed_redshift_gal']<1.05) & (hod['kind']==1)

plt.figure()
plt.scatter(hod['halo_lm'][hhc],hod['logf_halpha_model3_ext'][hhc],marker='.',s=1)
plt.title('sdhod centrals')
plt.ylim([-15.8,-14.6])
plt.xlabel('Log halo mass')
plt.ylabel('Log Halpha flux')

plt.figure()
plt.scatter(hod['halo_lm'][hhs],hod['logf_halpha_model3_ext'][hhs],marker='.',s=1)
plt.title('sdhod satellites')
plt.ylim([-15.8,-14.6])
plt.xlabel('Log halo mass')
plt.ylabel('Log Halpha flux')

plt.show()
