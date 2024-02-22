import napari
import shared.display as dis
import matplotlib.pyplot as plt
import seaborn as sns
import tifffile as tif
from skimage.measure import regionprops
import numpy as np
import seaborn as sns
import shared.segmentation as seg
import pandas as pd
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20231104_analysis_ColoDMandHSR_Wee1_6hr/"

sample = 'ColoHSR_Wee1_degrader_gH2AX_pH3'
conc_lst = ['XY37', 'XY38', 'XY41', 'XY42']

red_cut = 6.5

for conc in conc_lst:
    df = pd.read_csv('%s/%s_%s.txt' % (master_folder, sample, conc), na_values=['.'], sep='\t')
    cellcycle_lst = []
    df['ln_mean_int_gH2AX'] = np.log(df['gH2AX'])
    df['ln_mean_int_green'] = np.log(df['green'])
    df['ln_mean_int_red'] = np.log(df['red'])
    df['ln_mean_int_pH3'] = np.log(df['pH3'])
    for i in range(len(df)):
        if df['ln_mean_int_pH3'][i] < 6:
            cellcycle_lst.append('NA')
        elif (df['ln_mean_int_red'][i] < red_cut) & (df['ln_mean_int_green'][i] < 6.5) & (df['ln_mean_int_pH3'][i] <= 8):
            cellcycle_lst.append('G1S')
        elif (df['ln_mean_int_red'][i] < red_cut) & (df['ln_mean_int_green'][i] >= 6.5) & (df['ln_mean_int_pH3'][i] <= 8):
            cellcycle_lst.append('S')
        elif (df['ln_mean_int_red'][i] >= red_cut) & (df['ln_mean_int_green'][i] < 7) & (df['ln_mean_int_pH3'][i] <= 8):
            cellcycle_lst.append('G1')
        elif (df['ln_mean_int_red'][i] >= red_cut) & (df['ln_mean_int_green'][i] >= 7) & (df['ln_mean_int_pH3'][i] <= 8):
            cellcycle_lst.append('G2')
        elif df['ln_mean_int_pH3'][i] > 8:
            cellcycle_lst.append('M')
        else:
            print(df[i])

    df['cellcycle'] = cellcycle_lst
    df.to_csv('%s/%s_%s_cellcycle.txt' % (master_folder, sample, conc), index=False, sep='\t')


