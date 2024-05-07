import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob
import itertools
import seaborn as sns
import csv

from scipy.stats import linregress
from cycler import cycler

### Matplotlib settings for plotting
plt.rc('lines', linewidth=3, markersize=4)
plt.rc('axes', labelweight='heavy', labelsize=22, titlesize=22) 
plt.rc('axes.spines', top=False, right=False) 
plt.rc('legend', fontsize=16) 
plt.rc('xtick', labelsize=12) 
plt.rc('ytick', labelsize=16) 
plt.rc('figure', titlesize=22, titleweight='heavy') 
custom_cycler = (cycler(color=['#ff7f00', 'lightseagreen', 'mediumvioletred']))

### Create a NumPy array with time points in minutes
times = np.array([0.167, 1.0, 10.0, 60.0, 120.0]) 


def read_to_df(file):
    """Read and create a pandas DataFrame for the given argument.
    
    Args:
        file: There are four possible options:
            'segs' - peptide segments
            'expt' - experimental HDX deuterated fractions
            'pred' - calculated HDX deuterated fractions
            'reweighted' - reweighted HDX deuterated fractions
    
    Returns:
        df: A pandas DataFrame containing data for the given argument.
    """
    if file == 'segs':
        # Read and create a pandas DataFrame using a residue segments file
        df = pd.read_csv(os.path.expandvars('/file/path/to/HDXMS_expt_data/system-name_residue_segs.txt'),
                         sep='\s+', header=None, names=['ResStr', 'ResEnd'])
    elif file == 'expt':
        # Read and create a pandas DataFrame using an experimental deuterated fractions file
        df = pd.read_csv(os.path.expandvars('/file/path/to/HDXMS_expt_data/system-name_expt_dfracs.dat'), 
                         sep='\s+', skiprows=[0], header=None, usecols=[2, 3, 4, 5, 6], names=times)
    elif file == 'pred':
        # Read and create a pandas DataFrame using a computed deuterated fractions file
        df = pd.read_csv(os.path.expandvars('/file/path/to/CalcHDX/system-name_SUMMARY_segment_average_fractions.dat'), 
                         sep='\s+', skiprows=[0], header=None, usecols=[2, 3, 4, 5, 6], names=times)
    elif file == 'weight':
        df = pd.read_csv(os.path.expandvars('/file/path/to/CalcHDX/system-name_weight.csv'), 
                    header=None, usecols=[1, 2, 3, 4, 5], names=times)
    elif file == 'norm_weight':
        df = pd.read_csv(os.path.expandvars('/file/path/to/CalcHDX/system-name_norm_weight.csv'), 
                    header=None, usecols=[1, 2, 3, 4, 5], names=times)
    else:
        print("Incorrect argument given. Please choose one of the following: 'segs' 'expt' 'pred', 'weight', 'norm_weight'")
    return df


def plot_dfracs(*args):
    """Plot HDX deuterated fractions for each time point.
    
    Args:
        *args: 
            'expt' - experimental HDX deuterated fractions
            'pred' - computed HDX deuterated fractions
            'weight' - weighted HDX deuterated fractions
    """
    fig, axs = plt.subplots(len(times), 1, figsize=(18, 21), sharex=True)
    for i, (ax, t) in enumerate(zip(axs, times)):
        ax.set_prop_cycle(custom_cycler)
        for arg in args:
            if arg in ['expt', 'pred', 'weight', 'norm_weight']:
                xs = np.arange(0, read_to_df(arg).iloc[:, 0].shape[0])
                ax.plot(xs, read_to_df(arg).iloc[:, i], label=arg)
                ax.legend(bbox_to_anchor=(1.04,1), loc='upper left', borderaxespad=0)
                ax.set_xticks(xs)
                ax.set_xlim(xs[0], xs[-1])
                ax.set_xticklabels(PEPTIDES, rotation=90)
                ax.xaxis.set_tick_params(labelsize=18)
                ax.yaxis.set_tick_params(labelsize=18)
                ax.set_ylim(0, 1.01)
                ax.set_ylabel(str(t) + ' min', fontsize='xx-large')
            else:
                print("Incorrect argument given. Please choose one or more of the following: 'expt' 'pred', 'weight', 'norm_weight'")
    fig.text(0.5, 0.001, 'Peptide', ha='center', fontsize=24)
    fig.text(0.001, 0.5, 'HDX Deuterated Fractions', va='center', rotation='vertical', fontsize=24)
    fig.tight_layout(pad=3.5)
    fig.savefig('system-name_AS_dfracs_wRMSE.pdf')
    

### Read in the experimental and calculated HDX-MS data
expt_data = read_to_df('expt')
pred_data = read_to_df('pred')
res_segs = read_to_df('segs')

### Flatten the data arrays to combine all the labeling times together in a single dimension
expt_data_alltimes = expt_data.values.flatten()
pred_data_alltimes = pred_data.values.flatten()
res_segs_alltimes = res_segs.values.flatten()

### Make a list of peptides and timepoints
peptide_list = list(zip(res_segs_alltimes[::2], res_segs_alltimes[1::2]))
peptides_iter = itertools.repeat(peptide_list, len(times))
peptides_iterlist = (list(itertools.chain.from_iterable(peptides_iter)))
peptides_iterlist.sort(key=lambda x: x[:])
peptides_space = ','.join(map(lambda x: str(x[0]) + '-' + str(x[1]), peptides_iterlist))
peptides = peptides_space.split(",")
timepoint_list = list(times)
timepoints = timepoint_list*len(peptide_list)

### Calculate the non-normalized and normalized weights for each peptide and timepoint based on the absolute value of their residuals 
expt_pred_diff = expt_data_alltimes - pred_data_alltimes
expt_pred_absdiff = np.absolute(expt_pred_diff)
sum_absdiff = sum(expt_pred_absdiff)
absdiff_norm = expt_pred_absdiff/sum_absdiff

### Generate a dataframe with HDX data, residuals, and normalized weights for each peptide and timepoint
zipped = list(zip(peptides, timepoints, expt_data_alltimes, pred_data_alltimes, expt_pred_absdiff, absdiff_norm))
data = pd.DataFrame(zipped, columns=['Peptide', 'Timepoint', 'Experimental', 'Predicted', '|Expt-Pred|', 'Norm Weight'])

### Write out a .csv file containing the non-normalized weights for each peptide across all timepoints
norm_weight = data['|Expt-Pred|'].values.reshape(36,5)
nw = pd.DataFrame(weight)
nw.to_csv('system-name_weight.csv', header=None)

### Write out a .csv file containing the normalized weights for each peptide across all timepoints
norm_weight = data['Norm Weight'].values.reshape(36,5)
nw = pd.DataFrame(norm_weight)
nw.to_csv('system-name_norm_weight.csv', header=None)

### Verify that sum of normalized peptide weights equals 1 
pep_weight = data.groupby('Peptide', as_index=False, sort=False)['|Expt-Pred|'].sum()
weight_sum = np.sum(pep_weight)
pep_weight['|Expt-Pred| norm'] = pep_weight['|Expt-Pred|'].div(sum_absdiff)
pep_weight_norm = pep_weight['|Expt-Pred| norm']
print('Sum of normalized weights: ', np.sum(pep_weight_norm))

### Calculate RMSE and wRMSE
rmse = np.sqrt(np.mean( (pred_data_alltimes - expt_data_alltimes)**2 ))
wrmse = np.sqrt(np.sum (absdiff_norm*(expt_pred_absdiff**2)) )
print("The RMSE between computed and experimental HDX-MS is: %3.4f" % rmse)
print("The wRMSE between computed and experimental HDX-MS is: %3.4f" % wrmse)

### Write out RMSE and wRMSE values
with open("rmse.txt", "w") as file1:
    file1.writelines("The RMSE between computed and experimental HDX-MS is: %f" % rmse "\n", "The wRMSE between computed and experimental HDX-MS is: %f" % wrmse "\n")

### Read in the weights file and combine labeling times in a single dimension
weight_data = read_to_df('weight')
weight_data_alltimes = weight_data.values.flatten()

### Read in the normalized weights file and combine labeling times in a single dimension
norm_weight_data = read_to_df('norm_weight')
norm_weight_data_alltimes = norm_weight_data.values.flatten()    

### Plot peptide deuterated fractions for each timepoint w/ weights
plot_dfracs('expt', 'pred', 'weight', 'norm_weight')
