import pandas as pd
import MDAnalysis as mda
import matplotlib.pyplot as plt

from MDAnalysis.analysis import align, rms
from MDAnalysis.tests.datafiles import PSF, DCD

### Load topology/trajectory files 
top = 'system-name.psf'
traj = 'system-name.dcd'

### Select the atoms of interest 
sel = 'resid 1:100 and name CA'

### Align atomic coords to first frame and save avg positions over traj as ref structure
u = mda.Universe(top, traj)
average = align.AverageStructure(u, u, select=sel,
                                 ref_frame=0).run()
ref = average.results.universe

### Align the trajectory to the first frame of the simulation
alignment = align.AlignTraj(u, ref,
                          select=sel,
                          filename='system-name_aligned_traj.dcd',
                          in_memory=False).run()
u = mda.Universe(top, 'system-name_aligned_traj.dcd')

### Select protein atoms of interest 
sel = u.select_atoms(sel)
print('atom selection range: ', sel[0], sel[-1])

### Append data to df
rmsf_df = pd.DataFrame(columns = ['#Res', 'AtomicFlx'])
for i in range(len(sel.resids)):
    rmsf_df = rmsf_df.append({
    "#Res" : sel.resids[i],
    "AtomicFlx" : rmsf.results.rmsf[i]
    }, ignore_index=True)

### Round decimal to 4 digits
rmsf_df = rmsf_df.round(decimals=4)

### Write out data to .csv file
rmsf_df.to_csv("system-name_rmsf.dat", sep='\t', index=False, encoding='utf8', header=True)
