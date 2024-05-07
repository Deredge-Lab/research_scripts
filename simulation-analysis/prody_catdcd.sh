#!/bin/bash

### File paths
dcd_folder='/home/dderedge/Desktop/'
psf='/home/dderedge/Desktop/step3_input.psf'      # Must have same number of atoms as DCD
pdb='/home/dderedge/Desktop/step3_input.pdb'      # Must have same number of atoms as DCD
outfile='/home/dderedge/Desktop/system_prot.dcd'

### DCD files
prefix='prod_'       

### User selections
selection='protein'   # Default is all
align='protein'       # Default is protein
stride='1'            # Default is 1
first='0'             # Default is 0 (first frame)
last='-1'             # Default is -1 (last frame)

### Sort dcd files numerically
cd $dcd_folder
files=()
for file in $(ls -v $prefix*.dcd); do files=(${files[*]} "$file"); done
echo "${files[*]}"

### Activate conda environment where ProDy is installed
source activate openMM

### Concatenate traj files in order
prody catdcd -o $outfile -s $selection --psf $psf --pdb $pdb --first $first --last $last --stride $stride --align $align dcd ${files[*]}
