################################################################################
### Authors: Tiago Castro, Pierluigi Monaco                                  ###
###                                                                          ###
################################################################################
from astropy.io import fits as astrofits
from fitsio import FITS, FITSHDR
import numpy as np
import sys
import healpy as hp

if len(sys.argv)<2:
    print("Usage: pyton {} [my input file]".format(sys.argv[0]))
    sys.exit(0)
try: 
    input = __import__(sys.argv[1],  globals(), locals(), [], 0)
except ModuleNotFoundError:
    print("input file not found")
    print("Usage: pyton {} [my input file]".format(sys.argv[0]))
    sys.exit(0)


def write_catalog(fname, x, y, z, w, type="DATA", format="fits", coord="PSEUDO_EQUATORIAL"):
    if format=="txt":
        
        np.savetxt(fname, np.transpose( [x, y, z, w] ) )

    elif format=="fits":
        # current DM supported
        DM_VERSION = "2.5.0"

        # Extracted from AsciiToFits.py written by Daniele Tavagnacco
        
        if coord=="PSEUDO_EQUATORIAL":
            columns = [['SOURCE_ID', -1, 'f8'],
                       ['RA',         x, 'f8'],
                       ['DEC',        y, 'f8'],
                       ['REDSHIFT',   z, 'f8'],
                       ['WEIGHT',    -1, 'f8'],
                       ['DENSITY',    w, 'f8']]
        elif coord=="CARTESIAN":
            columns = [['SOURCE_ID', -1, 'f8'],
                       ['COMOV_X',    x, 'f8'],
                       ['COMOV_Y',    y, 'f8'],
                       ['COMOV_Z',    z, 'f8'],
                       ['WEIGHT',    -1, 'f8'],
                       ['DENSITY',    w, 'f8']]
 
        header_keywords = {"TELESCOP": "EUCLID",
                           "INSTRUME": "LE3GC-MOCKS",
                           "FILENAME": fname,
                           "CAT_ID"  : "MOCK",
                           "COORD"   : coord,
                           "ANGLE"   : "DEGREES"}

        extension="CATALOG"

        xmlKeys = {"pf"    : "PK_LE3_GC_WindowMultipoles",
                   "instr" : "LE3_GC_MOCKS",
                   "id"    : "MOCK",
                   "coord" : coord}

        print("Preparing FITS structure")
        types = []
        # keep just wanted columns ...bad but..
        for c in columns :
            types.append((c[0], c[2]))

        hdr = FITSHDR()
        print("+ Add keywords")
        for k in header_keywords :
            hdr[k] = header_keywords[k]

        keep_table = {}
        for c in columns :
            # if required but not existing (-1) fill with ones
            if (c[1] is -1) :
                # add tmp column with correct name and position, bu only ones
                print(str("+  -Col '%s' filled" % c[0]))
                keep_table[c[0]] = np.ones_like(columns[1][1],dtype=np.float64)
            else :
                # keep column with  requested position in the input file
                print(str("+  -Col '%s' from '%s'" % (c[0], c[1])))
                keep_table[c[0]] = c[1].astype(np.float64)

        fullsize = len(keep_table)*len(c[1])* 8 / 1024 / 1024.
        print(str("+   ~%.2f MB in memory" % fullsize))

        # now write some data
        print(str("+ Write FITS: %s" % fname))
        fits = FITS(fname,'rw', clobber=True)
        fits.write_table(data=keep_table, header=hdr, extname=extension)
        fits.close()

        print("+ Preparing XML product")

        with open (fname.replace(input.cat4le3_format,"xml"), "w+") as f :

            f.write('''<?xml version="1.0" encoding="UTF-8" standalone="no" ?>\n''')
            if (type == "RANDOM") :
                f.write('''<p1:DpdLE3gcInputRandCat xmlns:p1="http://euclid.esa.org/schema/dpd/le3/gc/inp/catrandin">\n''')
            else:
                f.write('''<p1:DpdLE3gcInputDataCat xmlns:p1="http://euclid.esa.org/schema/dpd/le3/gc/inp/catdatain">\n''')
            f.write('''  <Header>\n''')
            f.write(str('''    <ProductId>%s</ProductId>\n''' % header_keywords["FILENAME"].split('.')[0]))
            f.write('''    <ProductType>dpdLE3gcInputRandCat</ProductType>\n''')
            f.write('''    <SoftwareName>LE3_GC_test</SoftwareName>\n''')
            f.write('''    <SoftwareRelease>1.0</SoftwareRelease>\n''')
            f.write('''    <ManualValidationStatus>UNKNOWN</ManualValidationStatus>\n''')
            f.write('''    <PipelineRun>LE3_GC_Test_Inputs</PipelineRun>\n''')
            f.write('''    <ExitStatusCode>OK</ExitStatusCode>\n''')
            f.write(str('''    <DataModelVersion>%s</DataModelVersion>\n''' % DM_VERSION))
            f.write(str('''    <MinDataModelVersion>%s</MinDataModelVersion>\n''' % DM_VERSION))
            f.write('''    <ScientificCustodian>LE3</ScientificCustodian>\n''')
            f.write('''    <AccessRights>\n''')
            f.write('''      <EuclidConsortiumRead>true</EuclidConsortiumRead>\n''')
            f.write('''      <EuclidConsortiumWrite>true</EuclidConsortiumWrite>\n''')
            f.write('''      <ScientificGroupRead>true</ScientificGroupRead>\n''')
            f.write('''      <ScientificGroupWrite>true</ScientificGroupWrite>\n''')
            f.write('''    </AccessRights>\n''')
            f.write('''    <Curator>\n''')
            f.write('''      <Name>SDC-IT</Name>\n''')
            f.write('''    </Curator>\n''')
            f.write(str('''    <Creator>%s</Creator>\n''' % xmlKeys["pf"]))
            f.write('''    <CreationDate>2019-10-31T12:12:12Z</CreationDate>\n''')
            f.write('''  </Header>\n''')
            f.write('''  <Data>\n''')
            f.write(str('''  <Instrument>%s</Instrument>\n''' % xmlKeys["instr"]))
            f.write(str('''  <Catalog_ID>%s</Catalog_ID>\n'''% xmlKeys["id"]))
            f.write(str('''  <CoordType>%s</CoordType>\n''' % xmlKeys["coord"]))
            f.write('''  <Catalog format="le3.gc.cat.test" version="0.2">\n''')
            f.write('''    <DataContainer filestatus="PROPOSED">\n''')
            f.write(str('''      <FileName>%s.fits</FileName>\n''' % fname))
            f.write('''    </DataContainer>\n''')
            f.write('''  </Catalog>\n''')
            f.write('''  </Data>\n''')
            if (type == "RANDOM") :
                f.write('''</p1:DpdLE3gcInputRandCat>\n''')
            else :
                f.write('''</p1:DpdLE3gcInputDataCat>\n''')
            f.close()

        print("files %s and %s written"%(fname,fname.replace(input.cat4le3_format,"xml")))

    else:
        
        print("ERROR: unrecognized format in write_catalog")
        sys.exit(-1)

###################################################################################


print("# Running writeCatalog4LE3.py with {}".format(sys.argv[1]))

dndz = astrofits.getdata(input.dndz_fname())

if input.WriteLE3Random:
    print("# reading random catalog {}...".format(input.random_fname()))
    randomcat = astrofits.getdata(input.random_fname())

    if input.selection_random_tag is not None:
        fname=input.selection_random_fname()
        print("# reading selection of random from file {}...".format(fname))
        selection = astrofits.getdata(fname)['SELECTION']
    else:
        selection = np.ones(len(randomcat), dtype=bool)

    for zshell in input.finalCatZShell:

        zmin = zshell[0]
        zmax = zshell[1]

        print("# selecting the shell at z=%3.1f-%3.1f..."%(zmin,zmax))
        sel = selection & (randomcat[input.redshift_key] >= zmin) & (randomcat[input.redshift_key] < zmax)
        print("# computing weights...")
        density = np.interp(randomcat[input.redshift_key][sel], dndz['z_center'], dndz['N_gal']/dndz['bin_volume'])

        print("# writing random file {}".format(input.LE3_random_fname(zmin,zmax)))
        write_catalog(input.LE3_random_fname(zmin,zmax),
                      randomcat['ra_gal'][sel], randomcat['dec_gal'][sel], 
                      randomcat[input.redshift_key][sel],
                      density, type = "RANDOM", format = 'fits')

    del randomcat


print("# reading data catalog {}...".format(input.galcat_fname()))
cat = astrofits.getdata(input.galcat_fname())
zused = cat[input.redshift_key]
if input.selection_data_tag is not None:
    fname = input.selection_data_fname()
    print("# loading selection {}...".format(fname))
    selection = astrofits.getdata(fname)['SELECTION']
else:
    selection = np.ones(zused.size,dtype=bool)

for zshell in input.finalCatZShell:

    zmin = zshell[0]
    zmax = zshell[1]

    print("# selecting the shell at z=%3.1f-%3.1f..."%(zmin,zmax))
    mysel = (zused >= zmin) & (zused < zmax) & selection
    print("# computing weights...")
    density = np.interp(zused[mysel], dndz['z_center'], dndz['N_gal']/dndz['bin_volume'])

    print("# writing data file {}".format(input.LE3_data_fname(zmin,zmax)))
    write_catalog(input.LE3_data_fname(zmin,zmax),
                  cat['ra_gal'][mysel], cat['dec_gal'][mysel], zused[mysel], 
                  density, type = "DATA", format = 'fits')


print("# DONE!")
