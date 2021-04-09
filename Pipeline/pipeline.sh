#!/bin/bash
#SBATCH -J Mocks
#SBATCH -p euclid
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=240000
#SBATCH -o batch-%j.out
#SBATCH -e batch-%j.err
#SBATCH --time=360

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/u/pmonaco/anaconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/u/pmonaco/anaconda3/etc/profile.d/conda.sh" ]; then
        . "/u/pmonaco/anaconda3/etc/profile.d/conda.sh"
    else
        export PATH="/u/pmonaco/anaconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<

conda activate py3
cd /euclid_data/pmonaco/EuclidMocks/Pipeline

#python -u applyFootprintToMasterCatalog.py input

#python -u createSmoothHOD.py input_m1_FO
#python -u createSmoothHOD.py input_m3_FO

#python -u extractGalaxyCatalogFromMaster.py input_m1_FO
#python -u extractGalaxyCatalogFromMaster.py input_m3_FO

# SDHOD, model 1, real space
cat input_template.py | sed -e 's\FOOTTAG\None\' | sed -e 's\LFMODEL\1\' | sed -e 's\CATTYPE\sdhod\' | sed -e 's\SHUFFLE\False\' | sed -e 's\SELDATA\None\' | sed -e 's\SELRAND\None\' | sed -e 's\RSDFLAG\False\' > my_input.py

python -u createSDHOD_Catalog.py my_input
python -u createRandom.py my_input
# python -u createSelection.py my_input
# python -u createSelection.py my_input random
python -u dndz.py my_input
python -u numbercounts.py my_input
python -u writeCatalogs4LE3.py my_input
python -u createPkScripts.py my_input

# SDHOD, model 1, redshift space
cat input_template.py | sed -e 's\FOOTTAG\None\' | sed -e 's\LFMODEL\1\' | sed -e 's\CATTYPE\sdhod\' | sed -e 's\SHUFFLE\False\' | sed -e 's\SELDATA\None\' | sed -e 's\SELRAND\None\' | sed -e 's\RSDFLAG\True\' > my_input.py

python -u createSDHOD_Catalog.py my_input
python -u createRandom.py my_input
# python -u createSelection.py my_input
# python -u createSelection.py my_input random
python -u dndz.py my_input
python -u numbercounts.py my_input
python -u writeCatalogs4LE3.py my_input
python -u createPkScripts.py my_input

# SDHOD, model 3, real space
cat input_template.py | sed -e 's\FOOTTAG\None\' | sed -e 's\LFMODEL\3\' | sed -e 's\CATTYPE\sdhod\' | sed -e 's\SHUFFLE\False\' | sed -e 's\SELDATA\None\' | sed -e 's\SELRAND\None\' | sed -e 's\RSDFLAG\False\' > my_input.py

python -u createSDHOD_Catalog.py my_input
python -u createRandom.py my_input
# python -u createSelection.py my_input
# python -u createSelection.py my_input random
python -u dndz.py my_input
python -u numbercounts.py my_input
python -u writeCatalogs4LE3.py my_input
python -u createPkScripts.py my_input

# SDHOD, model 3, redshift space
cat input_template.py | sed -e 's\FOOTTAG\None\' | sed -e 's\LFMODEL\3\' | sed -e 's\CATTYPE\sdhod\' | sed -e 's\SHUFFLE\False\' | sed -e 's\SELDATA\None\' | sed -e 's\SELRAND\None\' | sed -e 's\RSDFLAG\True\' > my_input.py

python -u createSDHOD_Catalog.py my_input
python -u createRandom.py my_input
# python -u createSelection.py my_input
# python -u createSelection.py my_input random
python -u dndz.py my_input
python -u numbercounts.py my_input
python -u writeCatalogs4LE3.py my_input
python -u createPkScripts.py my_input


