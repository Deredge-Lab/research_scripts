import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as PathEffects
import pyemma
import pandas as pd
import seaborn as sns

### Plotting defaults
plt.rcParams["axes.labelweight"] = "bold"
plt.rcParams["axes.titleweight"] = "bold"

### Load .csv files containing Rg values and inter-residue distances
Rg = pd.read_csv('system-name_rg.dat')
resdist = pd.read_csv('system-name_res1-res2_dist.dat')

### Append data to lists
Rg_data = []
resdist_data = []

x = 0 
for row in Rg.iterrows():
    data = Rg.at[x,'Rg'] 
    Rg_data.append(data)      
    x = x + 1

z = 0
for row in resdist.iterrows():
    data = resdist.at[z,'Distance'] 
    resdist_data.append(data)    
    z = z + 1

### Plot KDE of Rg vs. resid distance
### Optionally, plot points corresponding to experimental SAXS/FRET observations and x-ray crystal/NMR structure
fig, ax = plt.subplots()
sns.kdeplot(Rg_data, resdist_data, shade=True, cbar=False, cmap='viridis', ax=ax)
ax.set(xlim=(0,50), ylim=(0,50), title='System Name Probability Density')
ax.set_xlabel('Radius of Gyration (\u212B)', labelpad=10)
ax.set_ylabel('Res1 - Res2 Distance (\u212B)', labelpad=10)
#expt = ax.scatter(21.2, 38.9, color='#ed0dd9', marker='*', s=125, edgecolors='white', lw=0.75)
#pdb = ax.scatter(20.1, 35.5, color='#ee3e1e', marker='p', s=75, edgecolors='white', lw=0.75)
#legend = ax.legend((expt, pdb), ('Expt (SAXS/FRET)', 'PDB structure'), scatterpoints=1, ncol=1, fontsize=10, loc='upper left')
#for ha in ax.legend_.legendHandles:
#    ha.set_edgecolor("black")
#    ha.set_linewidth(1)
plt.savefig('kde_rgyr_resdist_system.png', dpi=300)
plt.close()

