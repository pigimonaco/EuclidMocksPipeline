################################################################################
### Authors: Tiago Castro, Pierluigi Monaco                                  ###
###                                                                          ###
################################################################################
import numpy as np
from astropy.io import fits
import healpy as hp
from scipy.ndimage import gaussian_filter1d
import sys
from os import path


if len(sys.argv)<2:
    print("Usage: pyton {} [my input file]".format(sys.argv[0]))
    sys.exit(0)
try: 
    input = __import__(sys.argv[1],  globals(), locals(), [], 0)
except ModuleNotFoundError:
    print("input file not found")
    print("Usage: pyton {} [my input file]".format(sys.argv[0]))
    sys.exit(0)


print("# Running dndz.py with {}".format(sys.argv[1]))

fname=input.galcat_fname()
print("# loading catalog {}...".format(fname))

if not path.exists(fname):
    print("ERROR: galaxy catalog {} does not exist".format(fname))
    sys.exit(-1)

cat = fits.getdata(fname)

# loads the survey footprint in equatorial coordinates
footprint_res, footprint_zrange, sky_fraction, footprint = input.read_footprint()
        
# redshift binning for computing dn/dz
print("# setting binnings...")

ztab = np.linspace(footprint_zrange[0], footprint_zrange[1], int( np.round( (footprint_zrange[1]-footprint_zrange[0])/input.deltazbin) + 1) )

if np.abs((ztab[1]-ztab[0])/input.deltazbin-1) > 1e-5:
    print("!The chosen redshift binning and redshift range are not compatible!")
    sys.exit(-1)

# comoving distances for bins
bin_edges = np.asarray([ input.LCDMmodel.comoving_distance(z).to(input.l_unit).value for z in ztab ])
# bin volume and center
bin_volume = 4. * np.pi/3. * (bin_edges[1:]**3 - bin_edges[:-1]**3) * sky_fraction
bin_center = 0.5 * (bin_edges[:-1] + bin_edges[1:])

# selection
if input.selection_data_tag is not None:
    print("# loading selection {}...".format(input.selection_data_fname()))
    mysel = fits.getdata(input.selection_data_fname())['SELECTION']
else:
    mysel = np.ones(len(cat),dtype=bool)


print("# Processing catalog...")

# Histogram

Ngal     = np.histogram(cat[input.redshift_key][mysel], bins=ztab)[0]
isCen    = cat['kind'][mysel]==0
Ncen     = np.histogram(cat[input.redshift_key][mysel][isCen], bins=ztab)[0]


## Writes on file
fname = input.dndz_fname()
print("# Writing results on file {}...".format(fname))

dndz = np.empty(Ngal.size,
                         dtype=[('N_gal', int), ('N_gal_gaus', np.float), 
                                ('N_cen', int), ('N_cen_gaus', np.float),
                                ('z_center', np.float), ('z_lower', np.float), ('z_upper', np.float),
                                ('bin_center', np.float), ('bin_lower', np.float), 
                                ('bin_upper', np.float), ('bin_volume', np.float)])

dndz['N_gal'] = Ngal;
dndz['N_cen'] = Ncen;
dndz['bin_center'] = bin_center;
dndz['bin_lower'] = bin_edges[:-1]; dndz['bin_upper'] = bin_edges[1:];
dndz['bin_volume'] = bin_volume;
dndz['z_center'] = (ztab[1:] + ztab[:-1])/2; dndz['z_lower'] = ztab[:-1]; dndz['z_upper'] = ztab[1:];

dndz['N_gal_gaus'] = gaussian_filter1d(Ngal, input.smoothing_length);
dndz['N_cen_gaus'] = gaussian_filter1d(Ncen, input.smoothing_length);

fits.writeto(fname, dndz, overwrite=True)

print("# done!")