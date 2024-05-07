import pandas as pd
import numpy as np

data = pd.read_csv("/file/path/to/HDXMS_expt_data/system-name_expt_dfracs.csv")  # Path to experimental dfracs file

### List of column headers in HDX-MS data file
data.columns = [
    "Protein", 
    "Start", 
    "End", 
    "Sequence", 
    "Modification", 
    "Fragment", 
    "MaxUptake", 
    "MHP", 
    "State", 
    "Exposure", 
    "File",
    "z",
    "RT",
    "Inten",
    "Center"
    ]
data.head()

### Ignore column headers
data = data.drop(columns=["Protein", "Modification", "Fragment", "MaxUptake", "MHP", "State", "Exposure", "File", "RT", "Inten", "Center"])

### Ignore duplicate peptide sequences
data.drop_duplicates(subset="Sequence",
                     keep='first', inplace=True)

### Verify that correct file is being read
data = data.loc[data['File'].str.startswith('system-name_seq_0s')]
print(data)

### Calculate avg protein charge state 
avg_charge = data["z"].mean()
d = {"Avg_z": [avg_charge]}
d2 = pd.DataFrame(data=d)
data = data.assign(Avg_z=d2)

### Write out file with avg charge state
data.to_csv('system-name_hdxms-data_clean.csv', sep=',', index=False)
