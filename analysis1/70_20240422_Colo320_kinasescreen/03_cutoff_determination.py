import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
import skimage.io as skio
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240422_analysis_Colo320_kinaseScreen/"
data_dir = "%sprocessed/" % master_folder
data_dir1 = "%sdata/" % master_folder
output_dir = "%sprocessed/" % master_folder

ref = pd.read_excel('%s/kinase_screen.xlsx' % master_folder, na_values=['.'])
print(ref.head())
hc = [5000, 40000]
cutoff = 2.95
n_fov = 9

batch = 'point5uM_24hr'
plate = 3
if plate == 1:
    DM_wells = ['D4', 'D6', 'D8', 'D10', 'D12']
    HSR_wells = ['D14', 'D16', 'D18', 'D20', 'D22']
    DM_samples = ['XY02', 'XY04', 'XY06', 'XY08', 'XY10']
    HSR_samples = ['XY12', 'XY14', 'XY16', 'XY18', 'XY20']
elif plate == 2:
    DM_wells = ['I4', 'I6', 'I8', 'I10', 'I12']
    HSR_wells = ['I14', 'I16', 'I18', 'I20', 'I22']
    DM_samples = ['XY119', 'XY117', 'XY115', 'XY113', 'XY111']
    HSR_samples = ['XY109', 'XY107', 'XY105', 'XY103', 'XY101']
else:
    DM_wells = ['N3', 'N5', 'N7', 'N9', 'N11']
    HSR_wells = ['N4', 'N6', 'N8', 'N10', 'N12']
    DM_samples = ['XY201', 'XY203', 'XY205', 'XY207', 'XY209']
    HSR_samples = ['XY202', 'XY204', 'XY206', 'XY208', 'XY210']

samples = DM_samples
wells = DM_wells

for i in range(len(samples)):
    print(i)
    sample = samples[i]
    well = wells[i]
    data = pd.read_csv('%s/%s/%s_%s/txt/%s.txt' % (output_dir, batch, batch, plate, sample), na_values=['.'], sep='\t')
    data['log10_emiRFP670'] = np.log10(data['emiRFP670'])

    fig, ax = plt.subplots(figsize=(9, 6))
    fig.subplots_adjust(right=0.8)
    sns.histplot(data=data, x='hoechst')
    plt.axvline(hc[0], 0, 1000, c='r')
    plt.axvline(hc[1], 0, 1000, c='r')
    if not os.path.exists("%s%s/%s_%s/hoechst_hist/" % (output_dir, batch, batch, plate)):
        os.makedirs("%s%s/%s_%s/hoechst_hist/" % (output_dir, batch, batch, plate))
    plt.savefig('%s%s/%s_%s/hoechst_hist/%s.pdf' % (output_dir, batch, batch, plate, sample))
    plt.show()

    fig, ax = plt.subplots(figsize=(9, 6))
    sns.histplot(data=data[(data['hoechst'] > hc[0]) & (data['hoechst'] < hc[1])], x='log10_emiRFP670', bins=20)  # annot=True
    plt.xlim([2, 6])
    plt.axvline(cutoff, 0, 1000, c='r')
    if not os.path.exists("%s%s/%s_%s/emiRFP670_hist/" % (output_dir, batch, batch, plate)):
        os.makedirs("%s%s/%s_%s/emiRFP670_hist/" % (output_dir, batch, batch, plate))
    plt.savefig("%s%s/%s_%s/emiRFP670_hist/%s.pdf" % (output_dir, batch, batch, plate, sample))
    plt.show()


