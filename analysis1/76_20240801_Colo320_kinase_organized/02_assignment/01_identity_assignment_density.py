import numpy as np
import pandas as pd
import utilities as uti
import os

# INPUT PARAMETERS
# file info
data_dir = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240801_analysis/processed/"
data_dir1 = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240801_analysis/notes/"
output_dir = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240801_analysis/summary/"
batch = '48hr_density'
plates = [1]

# PARAMETERS
exp = 'Colo320_GrayKinase_%s' % batch
hc = [5000, 40000]
n_fov = 9
pd_plates = pd.read_csv('%s/plates.txt' % data_dir1, na_values=['.'], sep='\t')

df = pd.DataFrame(columns=['screen', 'plate', 'sample', 'well',
                           'hoechst_cutoff_low', 'hoechst_cutoff_high', 'log10_emiRFP670_cutoff', 'red_cutoff', 'green_cutoff', 'check_red', 'check_green',
                           'hoechst_mean', 'n_total', 'n_filtered', 'n_neg', 'n_pos', 'fov_hoechst',
                           'n_G1', 'n_G1S', 'n_S', 'n_G2M', 'n_G2MG1',
                           'n_neg_G1', 'n_neg_G1S', 'n_neg_S', 'n_neg_G2M', 'n_neg_G2MG1',
                           'n_pos_G1', 'n_pos_G1S', 'n_pos_S', 'n_pos_G2M', 'n_pos_G2MG1'])

for plate in plates:
    total_wells = pd_plates[(pd_plates['batch'] == batch) & (pd_plates['plate'] == plate)]['total_wells'].tolist()[0]
    samples_lst = ['XY0%s' % (x + 1) for x in range(9)] + ['XY%s' % (x + 10) for x in range(201)]
    wells_lst = ['D%s' % (x + 3) for x in range(20)] + ['E%s' % (x + 3) for x in range(20)][::-1] + \
                ['F%s' % (x + 3) for x in range(20)] + ['G%s' % (x + 3) for x in range(20)][::-1] + \
                ['H%s' % (x + 3) for x in range(20)] + ['I%s' % (x + 3) for x in range(20)][::-1] + \
                ['J%s' % (x + 3) for x in range(20)] + ['K%s' % (x + 3) for x in range(20)][::-1] + \
                ['L%s' % (x + 3) for x in range(20)] + ['M%s' % (x + 3) for x in range(20)][::-1] + \
                ['N%s' % (x + 3) for x in range(10)]
    samples = samples_lst[:total_wells]
    wells = wells_lst[:total_wells]

    for i in range(len(samples)):
        sample = samples[i]
        print('plate: %s; sample: %s (%s/%s)' % (plate, sample, i+1, total_wells))
        # general info
        well = wells[i]
        # load data
        data = pd.read_csv('%s/%s/%s_%s/txt/%s.txt' % (data_dir, batch, batch, plate, sample), na_values=['.'], sep='\t')
        # calculate log10_emiRFP670_cutoff
        cutoff = uti.cell_cutoff(data, hc)
        # calculate cell cycle cutoff
        red_cutoff, green_cutoff, check_red, check_green = uti.cc_cutoff(data, hc)
        # survival
        hoechst_mean = np.mean(data['hoechst'])
        n_total = len(data)
        n_filtered = len(data[(data['hoechst'] > hc[0]) & (data['hoechst'] < hc[1])])
        n_neg = len(data[(data['hoechst'] > hc[0]) & (data['hoechst'] < hc[1]) & (data['log10_emiRFP670'] < cutoff)])
        n_pos = len(data[(data['hoechst'] > hc[0]) & (data['hoechst'] < hc[1]) & (data['log10_emiRFP670'] >= cutoff)])
        fov_hoechst_lst = []
        for j in range(n_fov):
            temp = data[data['fov'] == j]
            if len(temp) != 0:
                fov_hoechst = data[data['fov'] == j]['fov_hoechst'].tolist()[0]
                fov_hoechst_lst.append(fov_hoechst)
        fov_hoechst_mean = np.mean(fov_hoechst_lst)
        # cell cycle
        data_neg = data[(data['hoechst'] > hc[0]) & (data['hoechst'] < hc[1]) & (data['log10_emiRFP670'] < cutoff)].copy().reset_index(drop=True)
        data_pos = data[(data['hoechst'] > hc[0]) & (data['hoechst'] < hc[1]) & (data['log10_emiRFP670'] >= cutoff)].copy().reset_index(drop=True)
        data_G1, data_G1S, data_S, data_G2M, data_G2MG1 = uti.get_cc(data, red_cutoff, green_cutoff)
        data_neg_G1, data_neg_G1S, data_neg_S, data_neg_G2M, data_neg_G2MG1 = uti.get_cc(data_neg, red_cutoff, green_cutoff)
        data_pos_G1, data_pos_G1S, data_pos_S, data_pos_G2M, data_pos_G2MG1 = uti.get_cc(data_pos, red_cutoff, green_cutoff)

        df.loc[len(df.index)] = [exp, plate, sample, well,
                                 hc[0], hc[1], cutoff, red_cutoff, green_cutoff, check_red, check_green,
                                 hoechst_mean, n_total, n_filtered, n_neg, n_pos, fov_hoechst_mean,
                                 len(data_G1), len(data_G1S), len(data_S), len(data_G2M), len(data_G2MG1),
                                 len(data_neg_G1), len(data_neg_G1S), len(data_neg_S), len(data_neg_G2M), len(data_neg_G2MG1),
                                 len(data_pos_G1), len(data_pos_G1S), len(data_pos_S), len(data_pos_G2M), len(data_pos_G2MG1)]

if batch == '24hr_density':
    densities = [1] * 20 + [2] * 20 + [3] * 20 + [4] * 20 + [6] * 20 + [8] * 20 + \
            [10] * 20 + [12] * 20 + [14] * 20 + [16] * 20 + [-1] * 10  # -1 is control
    df['density'] = densities
elif batch == '48hr_density':
    densities = [1] * 20 + [2] * 20 + [3] * 20 + [4] * 20 + [5] * 20 + [6] * 20 + \
                [7] * 20 + [8] * 20 + [9] * 20 + [10] * 20 + [-1] * 10  # -1 is control
    df['density'] = densities
else:
    densities = []

if not os.path.exists("%s/01_summary/" % output_dir):
    os.makedirs("%s/01_summary/" % output_dir)
df.to_csv('%s/01_summary/%s.txt' % (output_dir, batch), index=False, sep='\t')
print("DONE!")



