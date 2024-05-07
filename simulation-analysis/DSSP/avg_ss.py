import math
import pandas as pd
import os
import numpy as np
import statistics

file_path = '/file/path/to/system-name_dssp.csv'    # Path to DSSP output file
save_loc = '/file/path/to/outdir'                   # Location to write out file w/ avg sec struct values (alpha helix, beta sheet, random coil)

class AvgFile:
    def __init__(self, *args):
        ...


def find_files(file_path):
    file_list = []
    for fname in os.listdir(file_path):
        if fname.endswith("system-name_dssp.csv"):
            avg_file = AvgFile()
            fp = str(file_path) + '/' + str(fname)
            df = pd.read_csv(fp, sep='\s+', header=7)
            df.columns = ['Residue', 'Resname', 'ChainID', 'H', 'E', 'C']
            setattr(avg_file, 'df', df)
            file_list.append(avg_file)
    print(len(file_list))
    return file_list

def get_avg(file_list):
    hvalues = []
    evalues = []
    cvalues = []
    for f in file_list:
        z=0
        for row in f.df.iterrows():
            hvalues.append(f.df.at[z, 'H'])
            evalues.append(f.df.at[z, 'E'])
            cvalues.append(f.df.at[z, 'C'])
            z = z+1 
    print(len(hvalues))
    avg_h = sum(hvalues)/len(hvalues)
    avg_e = sum(evalues)/len(evalues)
    avg_c = sum(cvalues)/len(cvalues)
    with open('system-name_dsspAVG.txt', 'w+') as doc:
        doc.write("AVG SS (Helix, Beta Strand, Coil):" + '\t' + str(avg_h) + '\t' + str(avg_e) + '\t' + str(avg_c) + '\n')
        doc.close()

file_list = find_files(file_path)
get_avg(file_list)
          

