#!/bin/bash

### Bash script for calculating hydrogen-deuterium exchange (HDX) rates from a simulated protein ensemble using the CalcHDX module of HDXer
### Link to install HDXer: https://github.com/Lucy-Forrest-Lab/HDXer

### To display help text, type: !python $HDXER_PATH/HDXer/calc_hdx.py -h
 

### Arguments are space-delimited
### Topology and trajectory files must have the same number of atoms
### If selection includes non-protein atoms, be sure to add the following argument to mopt: "{ 'protonly' : False }"
### Double check that all necessary arguments are listed in the command

### Mandatory inputs
top='system.pdb'                                             # path to topology/parameter file
traj='system.dcd'                                            # path to trajectory file(s) 
seg='cropped_seg.list'                                       # path to peptide segment list (used to calculate segment-averaged deuterated fractions)
method='BestVendruscolo'                                     # protection factor model (default is BestVendruscolo)
calchdx='/home/username/Desktop/HDXer/HDXer/calc_hdx.py'     # path to calc_hdx.py script

### Optional inputs
expt='system.txt'                                            # path to experimental HDX-MS data file
times='0.167 1.0 10.0 60.0 120.0'                            # labeling times used to calculate HDX-MS deuterated fractions (defaults to [ 0.167, 1.0, 10.0, 120.0 ])
out='system_'                                                # prefix for output files
log='system'                                                 # prefix for logfile (defaults to 'HDX_analysis.log')
sel='protein'                                                # atom selection string in MDTraj format (default = all)
start='1'                                                    # starting frame (default = 1)
end='-1'                                                     # ending frame (defaults to final frame)
stride='1'                                                   # stride trajectory every "n" frames (default = 1)                
chunk='1000'                                                 # number of frames to read in traj (default = 1000)
mopt='"{ 'save_detailed' : True }"'                          # method options (see HDXer.methods.BV.params or HDXer.methods.PH.params for additional info)
#apot=''                                                     # analysis options (see HDXer.analysis.Analyze.params attribute for details)

### Activate conda env with HDXer installed
conda activate HDXER_ENV

### Make subdirectory for CalcHDX outputs 
cd CalcHDX
echo "CalcHDX outputs will be written to $PWD"

### Run CalcHDX cmd (uncomment aopt if needed)
python $calchdx -p $top -t $traj -c $chunk -s $start -e $end -str $stride -seg $seg -m $method -exp $expt -dt $times -out $out -log $log -sel $sel -mopt $mopt # -aopt $aopt
