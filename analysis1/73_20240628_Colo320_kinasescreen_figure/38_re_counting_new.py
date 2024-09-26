import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
import skimage.io as skio
import math
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
# cutoff = 3.2
n_fov = 9

batch = 'point5uM_48hr'

skip = pd.read_csv('%s/skip.txt' % data_dir1, na_values=['.'], sep='\t')

plate = 3
ori_data = pd.read_csv('%s/%s/%s_%s.txt' % (output_dir, batch, batch, plate), na_values=['.'], sep='\t')
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


def get_ellipse(x, y, a, b, A, data, data_x, data_y):
    data_flt = data[(((data[data_x] - x) * math.cos(A) + (data[data_y] - y) * math.sin(A)) ** 2) / (a ** 2) + (
                ((data[data_x] - x) * math.sin(A) + (data[data_y] - y) * math.cos(A)) ** 2) / (
                                b ** 2) <= 1].copy()
    return data_flt


def get_ellipse1(x, y, a, b, A, data, data_x, data_y):
    c = math.sqrt(a ** 2 - b ** 2)
    data_flt = data[
        np.sqrt((data[data_x] - x + c * math.cos(A)) ** 2 + (data[data_y] - y + c * math.sin(A)) ** 2) + np.sqrt(
            (data[data_x] - x - c * math.cos(A)) ** 2 + (
                        data[data_y] - y - c * math.sin(A)) ** 2) - 2 * a <= 0].copy()
    return data_flt


df_fov = pd.DataFrame(columns=['screen', 'rep', 'group', 'plate', 'sample', 'well', 'cell', 'treatment', 'hoechst_cutoff', 'log10_emiRFP670_cutoff',
                               'fov', 'hoechst_mean', 'n_total', 'n_filtered', 'n_neg', 'n_pos', 'per_neg', 'per_pos', 'fov_hoechst'])
df = pd.DataFrame(columns=['screen', 'rep', 'group', 'plate', 'sample', 'well', 'cell', 'treatment', 'hoechst_cutoff', 'log10_emiRFP670_cutoff',
                   'hoechst_mean', 'n_total', 'n_filtered', 'n_neg', 'n_pos', 'per_neg', 'per_pos', 'fov_hoechst',
                           'n_G1', 'n_G1S', 'n_S', 'n_G2M', 'n_G2MG1', 'n_neg_G1', 'n_neg_G1S', 'n_neg_S', 'n_neg_G2M', 'n_neg_G2MG1',
                           'n_pos_G1', 'n_pos_G1S', 'n_pos_S', 'n_pos_G2M', 'n_pos_G2MG1', 'per_green', 'per_red'])


def get_cc(data):
    data_G1 = data[(data['log2_AzaleaB5'] > green_cutoff) & (data['log2_H2-3'] <= red_cutoff)]
    data_G1S = data[(data['log2_AzaleaB5'] <= green_cutoff) & (data['log2_H2-3'] < red_cutoff)]
    data_S = data[(data['log2_AzaleaB5'] < green_cutoff) & (data['log2_H2-3'] >= red_cutoff)]
    data_G2M = data[(data['log2_AzaleaB5'] >= green_cutoff) & (data['log2_H2-3'] - k * data['log2_AzaleaB5'] - b > 0)]
    data_G2MG1 = data[(data['log2_H2-3'] > red_cutoff) & (data['log2_H2-3'] - k * data['log2_AzaleaB5'] - b <= 0)]
    return data_G1, data_G1S, data_S, data_G2M, data_G2MG1


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
    plt.subplots(figsize=(9, 7))
    m = plt.hist(data[(data['hoechst'] > hc[0]) & (data['hoechst'] < hc[1])]['log10_emiRFP670'],
                 weights=np.ones(len(data[(data['hoechst'] > hc[0]) & (data['hoechst'] < hc[1])])) / len(data[(data['hoechst'] > hc[0]) & (data['hoechst'] < hc[1])]),
                 range=[2.5, 5], bins=60, color='w', edgecolor='black')
    test_index = m[0].tolist().index(np.max(m[0][:15]))
    test_val = m[1][test_index]
    cutoff = test_val + 0.3
    plt.close()

    for j in range(n_fov):
        fov_hoechst = ori_data[(ori_data['plate'] == plate)&(ori_data['sample'] == sample)]['fov_hoechst'].tolist()[0]
        temp = data[data['fov'] == j].copy().reset_index(drop=True)
        hoechst_mean = np.mean(temp['hoechst'])
        n_total = len(temp)
        n_filtered = len(temp[(temp['hoechst'] > hc[0]) & (temp['hoechst'] < hc[1])])
        n_neg = len(temp[(temp['hoechst'] > hc[0]) & (temp['hoechst'] < hc[1]) & (temp['log10_emiRFP670'] < cutoff)])
        n_pos = len(temp[(temp['hoechst'] > hc[0]) & (temp['hoechst'] < hc[1]) & (temp['log10_emiRFP670'] >= cutoff)])
        per_neg = n_neg / (n_filtered+0.01)
        per_pos = n_pos / (n_filtered+0.01)
        df_fov.loc[len(df_fov.index)] = [data['screen'][0], rep, group, plate, sample, well, cell, treatment, hc, cutoff,
                                         j + 1, hoechst_mean, n_total, n_filtered, n_neg, n_pos, per_neg, per_pos, fov_hoechst]

    hoechst_mean = np.mean(data['hoechst'])
    n_total = len(data)
    n_filtered = len(data[(data['hoechst'] > hc[0]) & (data['hoechst'] < hc[1])])
    n_neg = len(data[(data['hoechst'] > hc[0]) & (data['hoechst'] < hc[1]) & (data['log10_emiRFP670'] < cutoff)])
    n_pos = len(data[(data['hoechst'] > hc[0]) & (data['hoechst'] < hc[1]) & (data['log10_emiRFP670'] >= cutoff)])
    per_neg = n_neg/(n_filtered+0.01)
    per_pos = n_pos/(n_filtered+0.01)
    fov_hoechst_mean = np.mean(df_fov[df_fov['sample'] == sample]['fov_hoechst'])

    data['log2_H2-3'] = np.log2(data['H2-3'])
    data['log2_AzaleaB5'] = np.log2(data['AzaleaB5'])
    # print(len(data))
    data = data[(data['hoechst'] > hc[0]) & (data['hoechst'] < hc[1])].reset_index(drop=True)
    plt.subplots(figsize=(9, 7))
    m = plt.hist(data['log2_H2-3'], bins=50)
    red_index = m[0].tolist().index(np.max(m[0][:-2]))
    red_axis = m[1][red_index]
    if red_index >= 15:
        red_axis = 10.1
    # print(red_index)
    # print(red_axis)
    plt.close()

    plt.subplots(figsize=(9, 7))
    m = plt.hist(data['log2_AzaleaB5'], bins=50)
    green_index = m[0].tolist().index(np.max(m[0][:-2]))
    green_axis = m[1][green_index]
    if green_index >= 15:
        green_axis = 8.3
    # print(green_index)
    # print(green_axis)
    plt.close()

    green_cutoff = green_axis + 0.65
    red_cutoff = red_axis + 0.65
    k = (16 - red_cutoff) / (16 - green_cutoff)
    b = (red_cutoff - green_cutoff) * 16 / (16 - green_cutoff)

    per_green = len(data[data['log2_AzaleaB5'] < green_cutoff]) / len(data)
    per_red = len(data[data['log2_H2-3'] < red_cutoff]) / len(data)

    # print(len(data))
    data_neg = data[(data['hoechst'] > hc[0]) & (data['hoechst'] < hc[1]) & (data['log10_emiRFP670'] < cutoff)].copy().reset_index(drop=True)
    data_pos = data[(data['hoechst'] > hc[0]) & (data['hoechst'] < hc[1]) & (data['log10_emiRFP670'] >= cutoff)].copy().reset_index(drop=True)
    data_G1, data_G1S, data_S, data_G2M, data_G2MG1 = get_cc(data)
    data_neg_G1, data_neg_G1S, data_neg_S, data_neg_G2M, data_neg_G2MG1 = get_cc(data_neg)
    data_pos_G1, data_pos_G1S, data_pos_S, data_pos_G2M, data_pos_G2MG1 = get_cc(data_pos)

    df.loc[len(df.index)] = ['Colo320_GrayKinase', rep, group, plate, sample, well, cell, treatment, hc, cutoff, hoechst_mean,
                             n_total, n_filtered, n_neg, n_pos, per_neg, per_pos, fov_hoechst_mean,
                             len(data_G1), len(data_G1S), len(data_S), len(data_G2M), len(data_G2MG1), len(data_neg_G1), len(data_neg_G1S), len(data_neg_S),
                             len(data_neg_G2M), len(data_neg_G2MG1), len(data_pos_G1), len(data_pos_G1S), len(data_pos_S), len(data_pos_G2M), len(data_pos_G2MG1), per_green, per_red]

df.to_csv('%s/%s/%s_%s_update1.txt' % (output_dir, batch, batch, plate), index=False, sep='\t')
df_fov.to_csv('%s/%s/%s_%s_fov_update1.txt' % (output_dir, batch, batch, plate), index=False, sep='\t')



