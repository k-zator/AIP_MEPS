# AIP_MEPS
Code for generating custom isosurface MEPS cube files for the calculation of AIPs
Author: Katarzyna Zator, kz265

#### Environment setup.
The script can be installed in python >3.6

Clone the repository:

    git clone https://github.com/k-zator/AIP_MEPS.git
    cd AIP_MEPS

Then install this environment instead (based on python 3.6):

        conda env create -f environment.yml
        conda activate aip_map

#### Calculation.

The code is written as an executable script with one variable which specifies the molecule file 
for which to calculate the MEPS isosurfaces:

        python calculate_MEPS water.xyz

the script runs a DFT (B3LYP/6-31G*) single point calculation, then extracts electron density, MEPS 
and the specific isosurfaces (0.0020, 0.0300, 0.0104) as required for AIP_map module as presented at:
https://github.com/k-zator/AIP_map

The script produces 3 key files, in the format {name}_{isosurface}.cube, e.g. (water_0.0020.cube).
Do note, though the extension specifies "cube", those are not equivalent to the Gaussian cube files.


Who do I talk to?
Any queries please contact Katarzyna Zator, kz265@cam.ac.uk

License
© Katarzyna Zator, Maria Chiara Storer, Christopher Hunter at the University of Cambridge
This is released under an AGPLv3 license for academic use.

