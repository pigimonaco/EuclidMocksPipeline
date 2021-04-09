################################################################################
### Author: Tiago Castro, Pierluigi Monaco                                   ###
###                                                                          ###
################################################################################
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
import healpy as hp
import sys

if len(sys.argv)<2:
    print("Usage: pyton {} [my input file]".format(sys.argv[0]))
    sys.exit(0)
try: 
    input = __import__(sys.argv[1],  globals(), locals(), [], 0)
except ModuleNotFoundError:
    print("input file not found")
    print("Usage: pyton {} [my input file]".format(sys.argv[0]))
    sys.exit(0)

print("# Running plot_dndz.py with {}".format(sys.argv[1]))

if not input.PLOT:
    print("Plots are not required, exiting...")
    sys.exit(0)

# special behaviour
second_dndz_fname=None
#second_dndz_fname=input.outdir+'NumberCounts/dndz_MWext_hodcat_8614_100sqdeg_m3_cMdiemer19_smL5_obsz.fits'
CHECK_WITH_RANDOM=True

# survey footprint in equatorial coordinate
footprint_res, footprint_zrange, sky_fraction, footprint = input.read_footprint()
sky_coverage=input.sqdegonthesky * sky_fraction
print("This survey covers {} sq deg".format(sky_coverage))
del footprint

print("Reading dndz from file {}".format(input.dndz_fname()))
dndz=fits.getdata(input.dndz_fname())

fig=plt.figure(figsize=(8, 8))
plt.suptitle(input.exclude_dir(input.dndz_fname()))

gs = gridspec.GridSpec(2, 1, height_ratios=[2.5, 1], hspace=0)

panel1 = plt.subplot(gs[0])
panel2 = plt.subplot(gs[1])

panel1.set_xlim([footprint_zrange[0]-0.1,footprint_zrange[1]+0.1])
panel2.set_xlim([footprint_zrange[0]-0.1,footprint_zrange[1]+0.1])
if input.lf_model=='1':
    panel1.set_ylim([800, 1e4])
elif input.lf_model=='3':
    panel1.set_ylim([0, 0.55e4])
panel2.set_yscale('linear')
panel2.set_ylim([0.7,1.3])
panel2.set_xlabel(r'redshift')
panel1.set_ylabel(r'$dn/dz$, deg$^{-2}$ $(\Delta z)^{-1}$')
panel2.set_ylabel(r'residuals vs model')
#uuu=panel1.set_xticklabels('',visible=False)
plt.tight_layout(pad=1.5)

# redshift binning for computing dn/dz
ztab = np.append(dndz['z_lower'],dndz['z_upper'][-1])

Nmodel=[]
for z1,z2 in zip(ztab[:-1],ztab[1:]):
    if z1<0.1:
        Nmodel.append(0.)
    else:
        Nmodel.append(input.Pozzetti_dndz(z1,z2)/input.deltazbin)
Nmodel=np.asarray(Nmodel)
pos=Nmodel>0
panel1.plot(dndz['z_center'],Nmodel,label='model',c='k')
panel2.plot(dndz['z_center'],np.ones_like(Nmodel),c='k')

panel1.plot(dndz['z_center'],dndz['N_gal']/sky_coverage/input.deltazbin,label='data',c='red')
panel2.plot(dndz['z_center'][pos],dndz['N_gal'][pos]/sky_coverage/input.deltazbin/Nmodel[pos],c='red')
panel1.plot(dndz['z_center'],dndz['N_gal_gaus']/sky_coverage/input.deltazbin,label='smoothed data',c='orange')
panel2.plot(dndz['z_center'][pos],dndz['N_gal_gaus'][pos]/sky_coverage/input.deltazbin/Nmodel[pos],c='orange')

panel1.plot(dndz['z_center'],dndz['N_cen']/sky_coverage/input.deltazbin,'--',label='data, centrals',c='red')
panel1.plot(dndz['z_center'],(dndz['N_gal']-dndz['N_cen'])/sky_coverage/input.deltazbin,':',label='data, satellites',c='red')


if CHECK_WITH_RANDOM:
    print("Reading random catalog {}...".format(input.random_fname()))
    random = fits.getdata(input.random_fname())

    Ng=np.histogram(random['observed_redshift_gal'], bins=ztab)[0]/sky_coverage/input.deltazbin/input.alpha
    panel1.plot(dndz['z_center'],Ng,':',label='from random',c='cyan')
    panel2.plot(dndz['z_center'][pos],Ng[pos]/Nmodel[pos],':',c='cyan')


if second_dndz_fname is not None:
    print("Reading second dndz from file {}".format(second_dndz_fname))
    dndz2  = fits.getdata(second_dndz_fname)
    panel1.plot(dndz2['z_center'],dndz2['N_gal']/sky_coverage/input.deltazbin,label='second set',c='green')
    panel1.plot(dndz2['z_center'],dndz2['N_cen']/sky_coverage/input.deltazbin,'--',label='second set, centrals',c='green')
    panel1.plot(dndz2['z_center'],(dndz2['N_gal']-dndz2['N_cen'])/sky_coverage/input.deltazbin,':',label='second set, satellites',c='green')
    panel2.plot(dndz2['z_center'][pos],dndz2['N_gal'][pos]/sky_coverage/input.deltazbin/Nmodel[pos],c='green')

panel1.legend()

if input.SHOW:
    plt.show()

plt.savefig(input.plot_dndz_fname())
print("## written image in file {}".format(input.plot_dndz_fname()))


