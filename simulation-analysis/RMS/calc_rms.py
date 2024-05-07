import pandas as pd
import MDAnalysis as mda
import matplotlib.pyplot as plt

from MDAnalysis.analysis import align, rms
from MDAnalysis.tests.datafiles import PSF, DCD

### Load topology/trajectory files 
top = 'system-name.psf'
traj = 'system-name.dcd'

### Select protein atoms of interest (mdanalysis notation, endpoint inclusive)
sel = 'resid 1:100 and name CA'

### Align atomic coords to first frame and save avg positions over traj as ref structure
u = mda.Universe(top, traj)
average = align.AverageStructure(u, u, select=sel,
                                 ref_frame=0).run()
ref = average.results.universe

### Align traj to initial simulation coords
alignment = align.AlignTraj(u, ref,
                          select=sel,
                          filename='system-name_prot_aligned.dcd',
                          in_memory=False).run()
u = mda.Universe(top, 'system-name_prot_aligned.dcd')

### Select protein atoms of interest 
sel = u.select_atoms(sel)
print('atom selection range: ', sel[0], sel[-1])

### Calculate RMSD and RMSF
rmsd = rms.RMSD(sel).run()
rmsf = rms.RMSF(sel).run()

### Append data to df
rmsd_df = pd.DataFrame(columns = ['#Frame', 'RMSD'])
for i in range(len(rmsd.rmsd)):
    rmsd_df = rmsd_df.append({
    "#Frame" : rmsd.rmsd[i][0],
    "RMSD" : rmsd.results.rmsd[i]
    }, ignore_index=True)

rmsf_df = pd.DataFrame(columns = ['#Res', 'AtomicFlx'])
for i in range(len(sel.resids)):
    rmsf_df = rmsf_df.append({
    "#Res" : sel.resids[i],
    "AtomicFlx" : rmsf.results.rmsf[i]
    }, ignore_index=True)

### Round decimal to 4 digits
rmsd_df = rmsd_df.round(decimals=4)
rmsf_df = rmsf_df.round(decimals=4)

### Write out data to .csv files
rmsd_df.to_csv("system-name_rmsd.dat", sep='\t', index=False, encoding='utf8', header=True)
rmsf_df.to_csv("system-name_rmsf.dat", sep='\t', index=False, encoding='utf8', header=True)

### Plotting functions
plt.plot(rmsd.rmsd, rmsd.results.rmsd)
#plt.xlim(1, 50000000)
#plt.ylim(0, 10)
plt.xlabel('Frame')
plt.ylabel('RMSD ($\AA$)')
plt.savefig('system-name_rmsd.png', dpi=300)

plt.plot(sel.resids, rmsf.results.rmsf)
#plt.xlim(1, 100)
#plt.ylim(0, 10)
plt.xlabel('Residue')
plt.ylabel('RMSF ($\AA$)')
plt.savefig('system-name_rmsf.png', dpi=300)
