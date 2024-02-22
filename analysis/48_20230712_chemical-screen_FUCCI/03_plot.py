import napari
import shared.display as dis
import matplotlib.pyplot as plt
import tifffile as tif
from skimage.measure import regionprops
import numpy as np
import seaborn as sns
import shared.segmentation as seg
import pandas as pd
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230712_analysis_chemical-screen_FUCCI/"
data_dir = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

hue_order = ['G1', 'S', 'G2', 'M']
line_colors = [(220/255, 20/255, 60/255), (154/255, 205/255, 50/255), (255/255, 127/255, 80/255), (0.85, 0.35, 0.25)]

plate = 'HSR_FUCCI_2hr'
sample_lst = ['B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B9', 'B10', 'B11', 'C11', 'C10', 'C9', 'C8', 'C7', 'C6',
              'C5', 'C4', 'C3', 'C2', 'D2', 'D3', 'D4', 'D5', 'D6', 'D7', 'D8', 'D9', 'D10', 'D11', 'E11', 'E10',
              'E9', 'E8', 'E7', 'E6', 'E5', 'E4', 'E3', 'E2', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8', 'F9', 'F10',
              'F11', 'G11', 'G10', 'G9', 'G8', 'G7', 'G6', 'G5', 'G4', 'G3', 'G2']

data = pd.DataFrame(columns=['plate', 'sample', 'G1', 'S/G2', 'G2/M', 'G1/S'])

for sample in sample_lst:
    print(sample)
    df = pd.read_csv('%s%s/txt/%s.txt' % (data_dir, plate, sample), na_values=['.'], sep='\t')

    df['ln_mean_int_green'] = np.log(df['mean_int_green'])
    df['ln_mean_int_red'] = np.log(df['mean_int_red'])

    cellcycle = []
    for i in range(len(df)):
        if (df['ln_mean_int_red'][i] > 1.5) & (df['ln_mean_int_green'][i] < 2):
            cellcycle.append('G1')
        elif (df['ln_mean_int_red'][i] < 1.5) & (df['ln_mean_int_green'][i] > 2):
            cellcycle.append('S/G2')
        elif (df['ln_mean_int_red'][i] > 1.5) & (df['ln_mean_int_green'][i] > 2):
            cellcycle.append('G2/M')
        elif (df['ln_mean_int_red'][i] < 1.5) & (df['ln_mean_int_green'][i] < 2):
            cellcycle.append('G1/S')
        else:
            cellcycle.append('NA')

    df['cellcycle'] = cellcycle
    print(len(df[df['cellcycle'] == 'NA']))
    data.loc[len(data.index)] = [plate, sample, len(df[df['cellcycle']=='G1']), len(df[df['cellcycle'] == 'S/G2']),
                                 len(df[df['cellcycle'] == 'G2/M']), len(df[df['cellcycle'] == 'G1/S'])]
    df.to_csv('%s%s/txt/%s_cellcycle.txt' % (output_dir, plate, sample), index=False, sep='\t')

    plt.subplots(figsize=(9, 6))
    sns.scatterplot(data=df, y='ln_mean_int_green', x='ln_mean_int_red', color=(0.30, 0.30, 0.30), s=5, alpha=0.5)
    sns.scatterplot(data=df[df['cellcycle']=='G1'], y='ln_mean_int_green', x='ln_mean_int_red', color=line_colors[0], s=5, alpha=0.5)
    sns.scatterplot(data=df[df['cellcycle'] == 'S/G2'], y='ln_mean_int_green', x='ln_mean_int_red', color=line_colors[1],
                    s=5, alpha=0.5)
    sns.scatterplot(data=df[df['cellcycle'] == 'G2/M'], y='ln_mean_int_green', x='ln_mean_int_red', color=line_colors[2],
                    s=5, alpha=0.5)
    sns.scatterplot(data=df[df['cellcycle'] == 'G1/S'], y='ln_mean_int_green', x='ln_mean_int_red', color=line_colors[3],
                    s=5, alpha=0.5)
    plt.xlim([0, 5.5])
    plt.ylim([1, 5.5])
    if not os.path.exists("%s%s/cellcycle/" % (output_dir, plate)):
        os.makedirs("%s%s/cellcycle/" % (output_dir, plate))
    plt.savefig('%s/%s/cellcycle/%s.pdf' % (output_dir, plate, sample))
    plt.close()
data.to_csv('%s%s/summary.txt' % (output_dir, plate), index=False, sep='\t')
print("DONE!")
