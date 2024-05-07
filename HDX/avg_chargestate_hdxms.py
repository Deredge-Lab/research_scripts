import pandas as pd
import numpy as np

data = pd.read_csv("/file/path/to/system-name_hdxms-data.csv")

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

data = data.drop(columns=["Protein", "Modification", "Fragment", "MaxUptake", "MHP", "State", "Exposure", "File", "RT", "Inten", "Center"])
  
data.drop_duplicates(subset="Sequence",
                     keep='first', inplace=True)

data = data.loc[data['File'].str.startswith('system-name_seq_0s')]
print(data)

avg_charge = data["z"].mean()
d = {"Avg_z": [avg_charge]}
d2 = pd.DataFrame(data=d)
data = data.assign(Avg_z=d2)

data.to_csv('system-name_hdxms-data_clean.csv', sep=',', index=False)
