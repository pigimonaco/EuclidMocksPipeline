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
import filenames

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
labd='lookup'
labs='complete'
second_dndz_fname=None
#second_dndz_fname=input.project+'NumberCounts/dndz_flagship_8614_BenSC8_m3_smL5_truez.fits'
CHECK_WITH_RANDOM=True
COMPLETENESS=True

# survey footprint in equatorial coordinate
footprint_res, footprint_zrange, sky_fraction, footprint = input.read_footprint()
sky_coverage=input.sqdegonthesky * sky_fraction
print("This survey covers {} sq deg".format(sky_coverage))
del footprint

fname = filenames.dndz(input)
print("Reading dndz from file {}".format(fname))
dndz=fits.getdata(fname)

fig=plt.figure(figsize=(8, 8))
plt.suptitle(filenames.exclude_dir(fname))

gs = gridspec.GridSpec(2, 1, height_ratios=[2.5, 1], hspace=0)

panel1 = plt.subplot(gs[0])
panel2 = plt.subplot(gs[1])

panel1.set_xlim([footprint_zrange[0]-0.1,footprint_zrange[1]+0.1])
panel2.set_xlim([footprint_zrange[0]-0.1,footprint_zrange[1]+0.1])
if input.lf_model=='1':
    panel1.set_ylim([0, 1e4])
elif input.lf_model=='3':
    panel1.set_ylim([0, 0.55e4])
panel2.set_yscale('linear')
panel2.set_ylim([0.85,1.15])
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

panel1.plot(dndz['z_center'],dndz['N_gal']/sky_coverage/input.deltazbin,label=labd,c='red')
panel2.plot(dndz['z_center'][pos],dndz['N_gal'][pos]/sky_coverage/input.deltazbin/Nmodel[pos],c='red')
panel1.plot(dndz['z_center'],dndz['N_gal_gaus']/sky_coverage/input.deltazbin,label='smoothed '+labd,c='orange')
panel2.plot(dndz['z_center'][pos],dndz['N_gal_gaus'][pos]/sky_coverage/input.deltazbin/Nmodel[pos],c='orange')

panel1.plot(dndz['z_center'],dndz['N_cen']/sky_coverage/input.deltazbin,'--',label=labd+', centrals',c='red')
panel1.plot(dndz['z_center'],(dndz['N_gal']-dndz['N_cen'])/sky_coverage/input.deltazbin,':',label=labd+', satellites',c='red')


if second_dndz_fname is not None:
    print("Reading second dndz from file {}".format(second_dndz_fname))
    dndz2  = fits.getdata(second_dndz_fname)
    panel1.plot(dndz2['z_center'],dndz2['N_gal']/sky_coverage/input.deltazbin,label=labs,c='b')
    panel1.plot(dndz2['z_center'],dndz2['N_cen']/sky_coverage/input.deltazbin,'--',label=labs+', centrals',c='b')
    panel1.plot(dndz2['z_center'],(dndz2['N_gal']-dndz2['N_cen'])/sky_coverage/input.deltazbin,':',label=labs+', satellites',c='b')
    panel2.plot(dndz2['z_center'][pos],dndz2['N_gal'][pos]/sky_coverage/input.deltazbin/Nmodel[pos],c='b')


if CHECK_WITH_RANDOM:
    print("Reading random catalog {}...".format(filenames.random(input)))
    random = fits.getdata(filenames.random(input))

    Ng=np.histogram(random[input.redshift_key], bins=ztab)[0]/sky_coverage/input.deltazbin/input.alpha
    panel1.plot(dndz['z_center'],Ng,'-.',label='from random',c='cyan')
    panel2.plot(dndz['z_center'][pos],Ng[pos]/Nmodel[pos],'-.',c='cyan')

    if (not input.apply_dataselection_to_random) & (input.selection_random_tag is not None):
        print("Reading random selection {}...".format(filenames.selection_random(input)))        
        sel = fits.getdata(filenames.selection_random(input))['SELECTION']
        Ng=np.histogram(random[input.redshift_key][sel], bins=ztab)[0]/sky_coverage/input.deltazbin/input.alpha
        panel1.plot(dndz['z_center'],Ng,'-.',label='random with selection',c='yellow')
        panel2.plot(dndz['z_center'][pos],Ng[pos]/Nmodel[pos],'-.',c='yellow')
        



panel1.legend()

if input.SHOW:
    plt.show()

plt.savefig(filenames.plot_dndz(input))
print("## written image in file {}".format(filenames.plot_dndz(input)))

if COMPLETENESS and second_dndz_fname is not None:

    plt.figure()
    plt.plot(dndz['z_center'],dndz['N_gal']/dndz2['N_gal'])
    plt.xlabel('redshift')
    plt.ylabel('completeness')
    plt.ylim([0,1.1])
    plt.plot(dndz['z_center'],np.ones_like(dndz['z_center']),c='k')
    plt.show()
