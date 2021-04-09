################################################################################
### Authors: Tiago Castro, Pierluigi Monaco                                  ###
###                                                                          ###
################################################################################
from string import Template
from shutil import copyfile
import numpy as np

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

print("# Running createPkScript.py with {}".format(sys.argv[1]))

# Open a file: params
f = open('parameters_pk_slave_fits.ini',mode='r')
params = Template( f.read() )
f.close()

print("# starting...")

count=0
nscript=0

for z1, z2 in input.finalCatZShell:

    print("writing %s"%input.params_fname(z1,z2))

    if input.Lbox is None:
        clbox='true'
        lbox=0.
    else:
        clbox='false'
        lbox=input.Lbox
    f = open(input.params_fname(z1,z2), "w")
    f.write(params.substitute(GRID   = input.ngrid,
                              DATA   = input.exclude_dir(input.LE3_data_fname(z1, z2)), 
                              RANDOM = input.exclude_dir(input.LE3_random_fname(z1, z2)),
                              OUTPUT = 'Pks/'+input.exclude_dir(input.pk_fname(z1,z2)),
                              CLBOX  = clbox, 
                              LBOX   = lbox))
    f.close()

    if count==0:

        copyfile("eden_slave.sh",input.script_fname(nscript))
        print("Writing {}".format(input.script_fname(nscript)))
        script=open(input.script_fname(nscript),"a")

    script.write("\n")
    script.write("mkdir -p "+input.pk_fname(z1,z2)+"\n")
    script.write("E-Run LE3_GC_PowerSpectrum  LE3_GC_ComputePowerSpectrum --log-level=DEBUG --parfile=PkParams/"+input.exclude_dir(input.params_fname(z1,z2))+" --workdir=/euclid_data/pmonaco/EuclidMocks/Products"+"\n\n")

    count += 1
    if count==input.max_PKs_in_script:
        count=0
        script.write("echo '### {} DONE! ###'".format(nscript))
        script.close()
        nscript+=1

script.write("echo '### {} DONE! ###'".format(nscript))
script.close()

print("# DONE!")
