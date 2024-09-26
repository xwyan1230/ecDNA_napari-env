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

df = pd.DataFrame(columns=['sample', 'well', 'identity', 'n_filtered', 'percentage'])

batch = '5uM_24hr'
plates = [1, 2, 3]
for plate in plates:
    if plate == 1:
        neg_wells = ['D4', 'D6', 'D8', 'D10', 'D12']
        pos_wells = ['D14', 'D16', 'D18', 'D20', 'D22']
        neg_samples = ['XY02', 'XY04', 'XY06', 'XY08', 'XY10']
        pos_samples = ['XY12', 'XY14', 'XY16', 'XY18', 'XY20']
    elif plate == 2:
        neg_wells = ['I4', 'I6', 'I8', 'I10', 'I12']
        pos_wells = ['I14', 'I16', 'I18', 'I20', 'I22']
        neg_samples = ['XY119', 'XY117', 'XY115', 'XY113', 'XY111']
        pos_samples = ['XY109', 'XY107', 'XY105', 'XY103', 'XY101']
    else:
        neg_wells = ['N3', 'N5', 'N7', 'N9', 'N11']
        pos_wells = ['N4', 'N6', 'N8', 'N10', 'N12']
        neg_samples = ['XY201', 'XY203', 'XY205', 'XY207', 'XY209']
        pos_samples = ['XY202', 'XY204', 'XY206', 'XY208', 'XY210']

    samples = neg_samples
    wells = neg_wells

    for i in range(len(samples)):
        sample = samples[i]
        well = wells[i]
        data = pd.read_csv('%s/%s/%s_%s/txt/%s.txt' % (output_dir, batch, batch, plate, sample), na_values=['.'],
                           sep='\t')
        data['log10_emiRFP670'] = np.log10(data['emiRFP670'])
        n_filtered = len(data[(data['hoechst'] > hc[0]) & (data['hoechst'] < hc[1])])
        n_neg = len(data[(data['hoechst'] > hc[0]) & (data['hoechst'] < hc[1]) & (data['log10_emiRFP670'] < cutoff)])
        percent = n_neg/n_filtered
        df.loc[len(df.index)] = [sample, well, 'neg', n_filtered, percent]

    samples = pos_samples
    wells = pos_wells

    for i in range(len(samples)):
        sample = samples[i]
        well = wells[i]
        data = pd.read_csv('%s/%s/%s_%s/txt/%s.txt' % (output_dir, batch, batch, plate, sample), na_values=['.'],
                           sep='\t')
        data['log10_emiRFP670'] = np.log10(data['emiRFP670'])
        n_filtered = len(data[(data['hoechst'] > hc[0]) & (data['hoechst'] < hc[1])])
        n_neg = len(data[(data['hoechst'] > hc[0]) & (data['hoechst'] < hc[1]) & (data['log10_emiRFP670'] < cutoff)])
        percent = (n_filtered - n_neg) / n_filtered
        df.loc[len(df.index)] = [sample, well, 'pos', n_filtered, percent]

df.to_csv('%s/%s/%s_cutoff_qc.txt' % (output_dir, batch, batch), index=False, sep='\t')

plt.subplots(figsize=(9, 9))
sns.violinplot(data=df[df['n_filtered'] > 500], x='identity', y='percentage')
sns.swarmplot(data=df[df['n_filtered'] > 500], x='identity', y='percentage', color="white")
plt.ylim([0, 1.05])
if not os.path.exists('%s/%s/' % (output_dir, batch)):
    os.makedirs('%s/%s/' % (output_dir, batch))
plt.savefig('%s/%s/%s_cutoff_qc.pdf' % (output_dir, batch, batch))
plt.show()

