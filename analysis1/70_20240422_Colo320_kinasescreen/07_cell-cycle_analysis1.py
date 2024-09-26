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
green_cutoff = 10.8
red_cutoff = 9

batch = '5uM_48hr'
plate = 3
if plate != 3:
    samples = ['XY0%s' % (x+1) for x in range(9)] + ['XY%s' % (x+10) for x in range(191)]
    wells = ['D%s' % (x+3) for x in range(20)] + ['E%s' % (x+3) for x in range(20)][::-1] + \
            ['F%s' % (x+3) for x in range(20)] + ['G%s' % (x+3) for x in range(20)][::-1] + \
            ['H%s' % (x+3) for x in range(20)] + ['I%s' % (x+3) for x in range(20)][::-1] + \
            ['J%s' % (x+3) for x in range(20)] + ['K%s' % (x+3) for x in range(20)][::-1] + \
            ['L%s' % (x+3) for x in range(20)] + ['M%s' % (x+3) for x in range(20)][::-1]
else:
    samples = ['XY0%s' % (x + 1) for x in range(9)] + ['XY%s' % (x + 10) for x in range(201)]
    wells = ['D%s' % (x + 3) for x in range(20)] + ['E%s' % (x + 3) for x in range(20)][::-1] + \
            ['F%s' % (x + 3) for x in range(20)] + ['G%s' % (x + 3) for x in range(20)][::-1] + \
            ['H%s' % (x + 3) for x in range(20)] + ['I%s' % (x + 3) for x in range(20)][::-1] + \
            ['J%s' % (x + 3) for x in range(20)] + ['K%s' % (x + 3) for x in range(20)][::-1] + \
            ['L%s' % (x + 3) for x in range(20)] + ['M%s' % (x + 3) for x in range(20)][::-1] + \
            ['N%s' % (x + 3) for x in range(10)]

df = pd.DataFrame(columns=['screen', 'rep', 'group', 'plate', 'sample', 'well', 'cell', 'treatment', 'hoechst_cutoff',
                           'log10_emiRFP670_cutoff', 'log2_H2-3_cutoff', 'log2_AzaleaB5_cutoff',
                           'n_total', 'n_filtered', 'n_neg', 'n_pos',
                           'n_G1', 'per_G1', 'n_G1S', 'per_G1S', 'n_S', 'per_S', 'n_G2M', 'per_G2M',
                           'n_neg_G1', 'per_neg_G1', 'n_neg_G1S', 'per_neg_G1S', 'n_neg_S', 'per_neg_S', 'n_neg_G2M', 'per_neg_G2M',
                           'n_pos_G1', 'per_pos_G1', 'n_pos_G1S', 'per_pos_G1S', 'n_pos_S', 'per_pos_S', 'n_pos_G2M', 'per_pos_G2M'])

for i in range(len(samples)):
    print(i)
    sample = samples[i]
    well = wells[i]
    cell = ref[(ref['screen_plate'] == plate) & (ref['screen_well'] == well)]['cell'].tolist()[0]
    rep = ref[(ref['screen_plate'] == plate) & (ref['screen_well'] == well)]['rep'].tolist()[0]
    group = ref[(ref['screen_plate'] == plate) & (ref['screen_well'] == well)]['group'].tolist()[0]
    treatment = ref[(ref['screen_plate'] == plate) & (ref['screen_well'] == well)]['treatment'].tolist()[0]
    data = pd.read_csv('%s/%s/%s_%s/txt/%s.txt' % (output_dir, batch, batch, plate, sample), na_values=['.'], sep='\t')
    data['log10_emiRFP670'] = np.log10(data['emiRFP670'])
    data['log2_H2-3'] = np.log2(data['H2-3'])
    data['log2_AzaleaB5'] = np.log2(data['AzaleaB5'])

    data_flt = data[(data['hoechst'] > hc[0]) & (data['hoechst'] < hc[1])].copy().reset_index(drop=True)
    data_neg = data_flt[data_flt['log10_emiRFP670'] < cutoff].copy().reset_index(drop=True)
    data_pos = data_flt[data_flt['log10_emiRFP670'] >= cutoff].copy().reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(8, 8))
    sns.scatterplot(data=data_flt, x='log2_AzaleaB5', y='log2_H2-3', s=4, alpha=0.5)
    plt.axvline(x=red_cutoff, linestyle='--', color='r')
    plt.axhline(y=green_cutoff, linestyle='--', color='r')
    plt.xlim([7, 16.5])
    plt.ylim([9, 16.5])
    if not os.path.exists("%s%s/%s_%s/cellcycle/" % (output_dir, batch, batch, plate)):
        os.makedirs("%s%s/%s_%s/cellcycle/" % (output_dir, batch, batch, plate))
    plt.savefig('%s%s/%s_%s/cellcycle/%s.pdf' % (output_dir, batch, batch, plate, sample))
    plt.close()

    fig, ax = plt.subplots(figsize=(8, 8))
    sns.scatterplot(data=data_neg, x='log2_AzaleaB5', y='log2_H2-3', s=4, alpha=0.5)
    plt.axvline(x=red_cutoff, linestyle='--', color='r')
    plt.axhline(y=green_cutoff, linestyle='--', color='r')
    plt.xlim([7, 16.5])
    plt.ylim([9, 16.5])
    if not os.path.exists("%s%s/%s_%s/cellcycle/" % (output_dir, batch, batch, plate)):
        os.makedirs("%s%s/%s_%s/cellcycle/" % (output_dir, batch, batch, plate))
    plt.savefig('%s%s/%s_%s/cellcycle/%s_neg.pdf' % (output_dir, batch, batch, plate, sample))
    plt.close()

    fig, ax = plt.subplots(figsize=(8, 8))
    sns.scatterplot(data=data_pos, x='log2_AzaleaB5', y='log2_H2-3', s=4, alpha=0.5)
    plt.axvline(x=red_cutoff, linestyle='--', color='r')
    plt.axhline(y=green_cutoff, linestyle='--', color='r')
    plt.xlim([7, 16.5])
    plt.ylim([9, 16.5])
    if not os.path.exists("%s%s/%s_%s/cellcycle/" % (output_dir, batch, batch, plate)):
        os.makedirs("%s%s/%s_%s/cellcycle/" % (output_dir, batch, batch, plate))
    plt.savefig('%s%s/%s_%s/cellcycle/%s_pos.pdf' % (output_dir, batch, batch, plate, sample))
    plt.close()

    n_total = len(data)
    n_filtered = len(data_flt)
    n_neg = len(data_neg)
    n_pos = len(data_pos)
    cc = []
    for cc_df in [data_flt, data_neg, data_pos]:
        if len(cc_df) != 0:
            G2M = cc_df[(cc_df['log2_H2-3'] >= green_cutoff) & (cc_df['log2_AzaleaB5'] >= red_cutoff)]
            G1 = cc_df[(cc_df['log2_H2-3'] < green_cutoff) & (cc_df['log2_AzaleaB5'] >= red_cutoff)]
            G1S = cc_df[(cc_df['log2_H2-3'] < green_cutoff) & (cc_df['log2_AzaleaB5'] < red_cutoff)]
            S = cc_df[(cc_df['log2_H2-3'] >= green_cutoff) & (cc_df['log2_AzaleaB5'] < red_cutoff)]
            cc = cc + [len(G1), len(G1)/len(cc_df), len(G1S), len(G1S)/len(cc_df), len(S), len(S)/len(cc_df), len(G2M), len(G2M)/len(cc_df)]
        else:
            cc = cc + [0]*8

    df.loc[len(df.index)] = [data['screen'][0], rep, group, plate, sample, well, cell, treatment,
                             hc, cutoff, green_cutoff, red_cutoff,
                             n_total, n_filtered, n_neg, n_pos] + cc

df.to_csv('%s/%s/%s_%s_cc.txt' % (output_dir, batch, batch, plate), index=False, sep='\t')




