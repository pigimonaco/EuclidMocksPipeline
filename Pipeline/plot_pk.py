################################################################################
### Authors: Tiago Castro, Pierluigi Monaco                                  ###
###                                                                          ###
################################################################################
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from glob import glob
import matplotlib.gridspec as gridspec
import sys
import importlib

if len(sys.argv)<2:
    print("Usage: python {} [my input file]".format(sys.argv[0]))
    sys.exit(0)
try: 
    input = __import__(sys.argv[1],  globals(), locals(), [], 0)
except ModuleNotFoundError:
    print("input file {} not found".format(sys.argv[1]))
    print("Usage: python {} [my input file]".format(sys.argv[0]))
    sys.exit(0)

print("# Running plot_pk.py with {}".format(sys.argv[1]))

type1='sdhod'
type2=None #'flagship'

for z1,z2 in input.finalCatZShell:

    importlib.reload(input)
    #input.cat_type=type1
    fname = glob(input.pk_fname(z1,z2)+"/EUC_LE3_GCL_PK_*0Z_1D.fits")[0]
    print("Reading PK from file {}".format(fname))
    data1 = fits.getdata( fname )

    if type2 is not None:
        input.cat_type=type2
        #input.otherlab=None
        fname = glob(input.pk_fname(z1,z2)+"/EUC_LE3_GCL_PK_*0Z_1D.fits")[0]
        print("Reading PK from file {}".format(fname))
        data2 = fits.getdata( fname )
    

    for order in [0,2,4]:

        title = r"z=[{},{}], PK{}".format(z1,z2,order)
        fig  = plt.figure(figsize=(6, 8))
        plt.suptitle(title)
        gs   = gridspec.GridSpec(2,1,height_ratios=[1.5,1], hspace=0)
        plot = fig.add_subplot(gs[0])
        res  = fig.add_subplot(gs[1], sharex=plot)
        plt.xscale('log')
        #plt.yscale('log')
        plot.set_yscale('log')
        string=r'P_{}'.format(order)
        plot.set_ylabel(r'$'+string+'(k)$ ($h^{-3}$ Mpc)')
        res.set_xlabel(r'$k$ ($h$/Mpc)')
        res.set_ylabel('data2/data1 -1')
        res.set_ylim([-0.15,0.15])

        plot.plot(data1['K'], data1['PK{}'.format(order)], label=type1,c='r')
        plot.plot(data1['K'],-data1['PK{}'.format(order)], ':', c='r')

        if type2 is not None:

            plot.plot(data2['K'], data2['PK{}'.format(order)], label=type2,c='b')
            plot.plot(data2['K'],-data2['PK{}'.format(order)], ':', c='b')

            res.plot(data1['K'], data2['PK{}'.format(order)]/data1['PK{}'.format(order)]-1., c='r')

        res.plot(data1['K'],np.zeros_like(data1['K']),c='k')
        plot.legend()

        fname=input.plot_pk_fname(z1,z2)
        fig.savefig(fname)
        print("## written image in file "+fname)



if input.SHOW:
    plt.show()

print("# DONE!")
