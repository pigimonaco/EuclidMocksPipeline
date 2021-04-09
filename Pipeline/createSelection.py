################################################################################
### Authors: Tiago Castro, Pierluigi Monaco                                  ###
###                                                                          ###
################################################################################
import numpy as np
from astropy.io import fits
import healpy as hp
import sys
from os import path


def extinction_curve(z):

    if sel_input.extinction_curve is 'standard':
        wavel = 6562.8 * (1.0+z)
        mwl = np.array([3600.,4500.,5500.,6600.,8000.,12500.,16500.,22000.,35000.,48000.])
        mwA = np.array([1.531,1.324,1.000,0.748,0.482,0.282,0.175,0.112,0.058,0.023])
        return 3.086 * np.interp(wavel,mwl,mwA)
    else:
        return None


if len(sys.argv)<2:
    print("Usage: pyton {} [my input file]".format(sys.argv[0]))
    sys.exit(0)
try: 
    input = __import__(sys.argv[1],  globals(), locals(), [], 0)
except ModuleNotFoundError:
    print("input file not found")
    print("Usage: pyton {} [my input file]".format(sys.argv[0]))
    sys.exit(0)

print("# Running selection.py with {}".format(sys.argv[1]))


use_data=True
if len(sys.argv)>=3:
    if (sys.argv[2]=='1') | (sys.argv[2]=='random'):
        print("# Selection will be applied to the random catalog")
        use_data=False
    else:
        print("# WARNING: unrecognised command-line option {}".format(sys.argv[2]))

if use_data:
    print("# Selection will be applied to the data catalog")

if use_data:
    if input.selection_data_tag is None:
        print("No selection specified for galaxy catalog, exiting")
        sys.exit(0)
    else:
        sel_input_fname = 'sel_input_{}'.format(input.selection_data_tag)

else:
    if input.selection_random_tag is None:
        print("No selection specified for random catalog, exiting")
        sys.exit(0)
    else:
        sel_input_fname = 'sel_input_{}'.format(input.selection_random_tag)

try: 
    sel_input = __import__(sel_input_fname,  globals(), locals(), [], 0)
except ModuleNotFoundError:
    print("ERROR: {}.py file not found".format(sel_input_fname))
    sys.exit(0)


my_selections=['extinction','visibilitymask','fluxcut','central','satellite']

for key in sel_input.selection_keys:
    if key not in my_selections:
        print("Error: selection key {} not recognised".format(key))
        sys.exit(1)

if use_data:
    print("Opening galaxy catalog {}...".format(input.galcat_fname()))
    cat  = fits.getdata(input.galcat_fname())
else:
    print("Opening random catalog {}...".format(input.random_fname()))
    cat  = fits.getdata(input.random_fname())

Ngal = len(cat)

selection=np.ones((Ngal),dtype=bool)

print("The catalog contains {} galaxies".format(Ngal))

if 'flux_limit' in sel_input.selection_keys:
    my_flux_limit = sel_input.selection_logflux_limit
else:
    my_flux_limit = input.logflux_limit

if 'extinction' in sel_input.selection_keys:

    print("# applying extinction...")

    conv      = np.pi/180.

    # THIS PROCESS MAY BE HEAVY FOR THE RANDOM CATALOG, WE MAY SPLIT IT INTO SEVERAL SECTIONS

    fname = input.outdir + sel_input.extinctionmap_fname
    print("# loading reddening map {}...".format(fname))
    reddening = hp.read_map(fname, field=sel_input.extinctionmap_field)
    if sel_input.extinctionmap_res != hp.npix2nside(reddening.size):
        print("# resampling extinction map...")
        reddening = hp.ud_grade(reddening, sel_input.extinctionmap_res)
    print("# finding pixels for galaxies...")
    pix       = hp.ang2pix( sel_input.extinctionmap_res, np.pi/2. - cat['dec_gal']*conv, cat['ra_gal']*conv )
    print("# computing extinction...")
    ext_mag   = reddening[pix] * extinction_curve(cat[input.redshift_key])

    print("# constructing selection...")
    selection &= cat[input.my_flux_key] >= my_flux_limit + ext_mag / 2.5


elif 'flux_limit' in sel_input.selection_keys:

    print("# applying flux limit {}...".format(my_flux))

    selection &= cat[input.my_flux_key] >= sel_input.selection_logflux_limit


elif 'visibilitymask' in sel_input.selection_keys:
    
    print("# applying visibility mask {}...".format(sel_input.selection_VM_fname))

    VM = fits.getdata(sel_input.selection_VM_fname)
    VM = hp.ud_grade(VM, sel_input.selection_VM_res)

    conv      = np.pi/180.
    pix       = hp.ang2pix( sel_input.selection_VM_res, np.pi/2. - cat['dec_gal']*conv, cat['ra_gal']*conv )
 
    selection &= cat[input.my_flux_key] >= VM[pix]
    

# selection of centrals and satellites applies only to data catalogs
if use_data:
    if 'central' in sel_input.selection_keys:

        print("# applying selection of centrals...")

        selection &= cat['kind'] == 0

    elif 'satellite' in sel_input.selection_keys:

        print("# applying selection of satellites...")

        selection &= cat['kind'] == 1


if use_data:
    fname = input.selection_data_fname()
else:
    fname = input.selection_random_fname()

print("# Writing file {}...".format(fname))

tofits = np.empty(Ngal,dtype=[('SELECTION', bool)])
tofits['SELECTION'] = selection

fits.writeto(fname, tofits, overwrite=True)

print("# Done!")

