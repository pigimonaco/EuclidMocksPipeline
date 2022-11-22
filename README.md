# Pipeline for processing galaxy or halo mocks on the lightcone

## Installation

Setup a typical CONDA installation with python3. Activate your environment first. You may install the requirements using:
```
python3 -m pip install -r requirements.txt
```

Then you may run the installation of the euclid observational systematics package (`euclid_obssys`) using
```
python3 setup.py [install, develop]
```

If you want only use the pipeline, use "install", but if you want your changes to the scripts to be immediately active, use "develop"

Two basic tool are made available through `$PYTHON -m euclid_obssys.tool` and `$PYTHON -m euclid_obssys.view`. Each command comes with an help message describing the input/output and the optional arguments. Please consult `quick.sh` for a complete example.

## Starting point:

Build and install python package:
```
python2 setup.py install
```


### Download Flagship 2.0 catalog

You can obtain the selection of 2.0 Flagship mock needed by the pipeline as follows:
```
SELECT `ra_gal`, `dec_gal`, `x_gal`, `y_gal`, `z_gal`, `vx_gal`, `vy_gal`, `vz_gal`, `true_redshift_gal`, `observed_redshift_gal`, `logf_halpha_model3_ext`, `logf_hbeta_model3_ext`, `logf_o2_model3_ext`, `logf_n2_model3_ext`, `logf_o3_model3_ext`, `logf_s2_model3_ext`, `euclid_vis`, `euclid_nisp_h`, `euclid_nisp_j`, `euclid_nisp_y`, `abs_mag_r01`, `lm_halo`, `kind`, `halo_id`, `galaxy_id`, `bulge_fraction`, `bulge_r50`, `disk_r50`, `disk_axis_ratio` 
       FROM fs2_mock_std_spv3_epic_20220422 
       WHERE `true_redshift_halo`>=0.8 
       AND `true_redshift_halo`<2.0 
       AND `lm_halo` >=11 
```

The catalog will be provided as a fits file like "12345.fits", the query number is a relevant information that will be needed.

### Project and Repository directories

You will be asked to specify paths for a Repo directory, that contains the Flagship catalog and other important files, and a Project directory, where the product of your analysis will be contained. More details later.


### Example pipeline

A default script that runs an entire example pipeline is available in `example.sh`. You need to have downloaded the flagship halos in `Repo/RawCatalogs/catalogs_8614.fits`.


## Repo Directory Tree:
[repo]/RawCatalogs            # Cosmohub queries, from the database or with footprint applied
      /Footprints             # survey footprints (see below)
      /ExtinctionMaps         # extinction and reddening maps from various sources
      /SDHOD                  # Stripped-down HOD versions
      /LookUpTable	      # look-up tables to apply systematics
      /SelectionInputs	      # config files to obtain specific selections


## Project Directory Tree:
[project]/GalaxyCatalogs         # complete galaxy catalogs
         /RandomCatalogs         # complete random catalogs
         /Selections             # selections to be applied to galaxy catalogs
         /Catalogs4LE3           # catalogs readable by LE3 PFs, in redshift slices
         /NumberCounts           # galaxy counts
         /Cls                    # angular maps and angular clustering measurements
	 /PK                     # directory for OU-LE3 estimators, each directory
	 /2PCF                   #   contains: 
	 /BK			 #   /Scripts   where you find the scripts to launch the estimator
	 /3PCF			 #   /Params    where parameter files are stored
	 /CM-PK			 #   /Measures  where LE3 output is written
	 /CM-2PCF		 #


## Definitions:

* config file: most scripts require a config file that contains all
  needed information for ALL the scripts. The name of the config file
  is provided to the scripts as an argument.

* Survey footprint: a fits file with the fields FOOTPRINT_G,
  containing a boolean healpy map that is True/False in the sky pixels
  that are/are not observed, in galactic coordinates, and a redshift
  range of validity. It is defined by a footprint_tag in the config
  file. One can define in the config file a default footprint for a
  catalog, that would be called using footprint_tag=None.
  The footprint is read by the read_footprint() function defined in
  the config file, that returns the healpy resolution of the map, the
  redshift range in which the footprint is defined, the fraction of
  the sky covered by the survey and the footprint itself.

  NOTE: it is assumed below that sky coordinates are ALWAYS galactic.

* Master catalog: a catalog containing information on all dark matter
  halos and its central/satellite galaxies, coming from a CosmoHub
  query, or from applying a footprint on a Cosmohub query.

* Galaxy catalog: a catalog that contains only galaxies subject to a
  basic selection, like Halpha flux > 2e-16 erg/cm^2/s. The basic
  selection is based on a flux_key, that may be logf_halpha_model1_ext
  (model '1'), logf_halpha_model3_ext (model '3'), euclid_nisp_H
  (model '9') or halo_lm (model '0'). The threshold criterion is based
  on logflux_limit, specified in the config file. A galaxy catalog
  contains the following fields (with the same meaning as the Flagship
  catalog): x_gal, y_gal, z_gal, ra_gal, dec_gal, kind,
  true_redshift_gal, observed_redshift_gal, halo_lm, id, halo_id,
  [flux_key], sh_[flux_key]. This last column contains the same
  fluxes, shuffled among galaxies at fixed redshift to suppress
  luminosity-dependent bias.

* Selection: a boolean array of the same size of its reference galaxy
  or random catalog, that selects galaxies according to:
  - fluxcut: a simple cut in a reference variable (may be Halpha flux,
    H band magnitude, halo mass);
  - central / satellite: central or satellite galaxies;
  - extinction: galaxies are selected accorting to their extincted
    flux;
  - visibilitymask: a boolean visibility mask is applied to the catalog.
  - lookup: a lookup table is used to determine what galaxies are 
    detected.
  A selection is defined by a tag specified in the config file, and
  requires a separate sel_input_[tag].py file to be created. It can be
  applied to a galaxy or a random catalog. The config file specifies
  two selection tags, selection_data_tag for the galaxy catalog and
  selection_random_tag for the random.

* Random catalog: it is created by (randomly) replicating the data
  catalog alpha times in redshift and flux, and by distributing
  galaxies uniformly on the sky (within the footprint). It assumed to
  be unmodulated on the sky. It can follow a specific number density
  if required **TO BE IMPLEMENTED**.

* Catalogs for LE3: these are extracted from galaxy and random
  catalogs, give galaxies in specified redshift bins, and are readable
  by LE3 codes.

* Pinocchio light cones: the pipeline can process a set of pinocchio
  light cones. In this case it will be possible to compute the galaxy
  number density averaged over the mocks, and a random that follows
  this averaged number density.


## Preliminary scripts, to be run once:

* createIndicesForSats.py: must be run on a master catalog to create
  pointers of satellite galaxies to their main halo in that catalog.
  They are needed to create the SDHOD. This script must be edited to
  type in the name of the master catalog.

* createFullOctantFootprint.py: creates a survey footprint for the
  Flagship octant.

* createFootprint.py: this is an example script for creating a survey
  footprint and writing its file. It creates either a square or a circle
  field around a given position of the sky.

* applyFootprintToMasterCatalog.py [config file]: it applies a
  footprint to the master catalog, extracting another (smaller) master
  catalog.


## Galaxy catalogs from master:

* extractGalaxyCatalogFromMaster.py [config file]: extracts a galaxy
  catalog from a master catalog by applying the standard flux limit,
  according to the model. It must be run separately for each model.


## Stripped-Down HOD from master:

* createSmoothHOD.py [config file]: based on a master file, it measures
  the HOD curves and smooths them in redshift. It should be run on the
  full octant, there is no reason to run it many times for a specified
  model.


## Galaxy catalog created from master DM halos using SDHOD

* createSDHODCatalog.py [config file]: it creates a galaxy catalog by
  populating the halos of the master catalog with galaxies, using the
  SDHOD, including standard recipes to distribute satellites in halos.


## Galaxy catalog created from pinocchio light cones using SDHOD

* createSDHODfromPinocchio.py [config file] [first run] [last run]: it
  creates a set of galaxy catalog by populating the halos of pinocchio
  light cones with galaxies, using the SDHOD, including standard
  recipes to distribute satellites in halos. A calibration of halo
  masses is applied, to match them to the Flagship catalog either in
  abundance or in clustering amplitude.


## Random catalog

* dndz.py [config file]: given a galaxy catalog, it measures the number
  density of its galaxies (number of galaxies brighter than the flux
  limit in each redshift interval as a function of redshift). This is
  needed to create the random catalog. It also provides number
  densities smoothed in redshift. It accepts application of a
  selection as an option, though the random catalog will typically be
  created with the dndz of the unselected data catalog.

* createRandom.py [config file]: it creates a random for the galaxy
  catalog, with abundace alpha times the data catalog, alpha as
  specified in the config catalog. TO BE IMPLEMENTED: on request it can
  follow a specified galaxy number density.


## Selection

* createSelection.py [config file] [random]: creates a selection for a
  galaxy or a random catalog; the second option is chosen if "random"
  or "1" is provided when calling the script. The parameters of the
  selection are specified in the `sel_input_[input.selection_tag].py`,
  where `selection_tag` is given in the config file.


## Measurements of galaxy catalogs:

* writeCatalogs4LE3.py [config file]: given a selection, it writes
  galaxy and random catalogs in redshift shells, in files that are
  readable by LE3 processing functions.

* createPKScripts.py [config file]: creates the parameter files and the
  scripts needed to run the LE3 PK code.

* numbercounts.py [config file]: given a data catalog and a selection,
  it measures the differential number counts (number of galaxies per
  sq deg per redshift interval per magnitude as a function of flux) in
  redshift intervals.

* `maps_and_cls.py` [config file]: it builds the density contrast maps
  data and random galaxy catalogs, read from LE3 files, and measures
  their angular clustering with anafast.
  THIS IS STILL WORK IN PROGRESS


## Plots

THIS PART IS STILL UNDER DEVELOPMENT

* visualizeHOD.py [config file]: plots the Ncen(Mh,z) and Nsat(Mh,z)
  curves for the SDHOD.

* plot_dndz.py [config file]: plots the number density of galaxies as a
  function of redshift, also separated into centrals and satellites,
  and compares it to the relative Pozzetti model if relevant. It is
  easy to edit the script to plot two sets of curves.

* plot_numbercounts.py [config file]: plots the differential number
  counts as a function of flux for specified redshift intervals, also
  separated into centrals and satellites, and compares it to Pozzetti
  model if relevant.

* plot_pk.py [config file]: plot the measured PKs (monopole, quadrupole,
  exadecapole) in redshift intervals. It is easy to edit the script to
  plot two sets of curves.

* plot_cls.py [config file]: plot the measured angular clustering in redshift
  intervals. It is easy to edit the script to plot two sets of curves.
  THIS IS STILL WORK IN PROGRESS

* angular_map.py [filename] [selection] [zmin] [zmax]: reads the
  galaxy or random catalog specified in the command line argument, its
  selection ('None' to skip this) and plots its angular map of density
  contrast (with Nside=256) in the specified redshift interval. Very
  useful for quick checks of what a catalog is like.


