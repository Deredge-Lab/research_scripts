import pandas as pd
import MDAnalysis as mda

from MDAnalysis.analysis import align, rms
from MDAnalysis.tests.datafiles import PSF, DCD

### Load topology/trajectory files 
top = 'system-name.psf'
traj = 'system-name.dcd'

### Select protein residues/atoms (mdanalysis notation, endpoint inclusive)
sel = 'resid 1:100 and name CA'

### Align traj to initial simulation coords
u = mda.Universe(top, traj)
ref = mda.Universe(top, traj)
u.trajectory[-1]	# Set mobile trajectory to last frame
ref.trajectory[0]   # Set reference trajectory to first frame
alignment = align.AlignTraj(u, ref, select=sel, filename='system-name_aligned_frame1.dcd').run()
u = mda.Universe(top, 'system-name_aligned_frame1.dcd')

u.trajectory[-1]	
ref.trajectory[0]   

### Select protein atoms of interest 
u_res = ref.select_atoms(sel)
ref_res = ref.select_atoms(sel)
print('atom selection range: ', u_res[0], u_res[-1])

### Calculate RMSD
rmsd = rms.RMSD(u_res.positions, ref_res.positions)
rmsd.run()

### Write out .csv file with cols: 'Frame #' , 'Time' , 'RMSD' 
rmsd_df =  pd.DataFrame(rmsd.rmsd, columns = ['Frame', 'Time (ns)', 'RMSD (CA)'])
for i in range(len(rmsd.rmsd)):
    rmsd_df = rmsd_df.append({
    "#Frame" : rmsd.rmsd[i][0],
    "Time" : rmsd.rmsd[i][1],
    "RMSD" : rmsd.rmsd[i][2]
    }, ignore_index=True)
rmsd_df = rmsd_df.round(4)
rmsd_df.to_csv("system-name_rmsd.dat", sep='\t', index=False, encoding='utf8', header=True)

