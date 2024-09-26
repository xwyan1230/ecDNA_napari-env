import pandas as pd
import os
import numpy as np

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240517_analysis_KOscreen1/"
data_dir = "%sprocessed/" % master_folder
output_dir = "%ssummary/" % master_folder

version = 2

samples = ['B%s' % (x+2) for x in range(10)] + ['C%s' % (x+2) for x in range(10)] + \
          ['D%s' % (x+2) for x in range(10)] + ['E%s' % (x+2) for x in range(10)] + \
          ['F%s' % (x+2) for x in range(10)] + ['G%s' % (x+2) for x in range(10)]
hueorder = ['mCherry', 'GFP']
reps = [1, 2]

data = pd.DataFrame(columns=['sample', 'group', 'n', 'mean_area_nuclear', 'median_area_nuclear',
                             'mean_copy_number', 'median_copy_number',
                             'n_DNAFISH', 'mean_n_ecDNA', 'median_n_ecDNA',
                             'n_point02', 'mean_per_AUC', 'median_per_AUC',
                             'peak_relativer', 'fwhm_relativer', 'center_relativer', 'peakval_relativer',
                             'peak_absoluter', 'fwhm_absoluter', 'center_absoluter', 'peakval_absoluter', 'gene'])
data_mean = pd.DataFrame(columns=['sample', 'mean_area_nuclear', 'median_area_nuclear',
                                  'mean_copy_number', 'median_copy_number', 'mean_n_ecDNA', 'median_n_ecDNA',
                                  'mean_per_AUC', 'median_per_AUC',
                                  'peak_relativer', 'fwhm_relativer', 'center_relativer', 'peakval_relativer',
                                  'peak_absoluter', 'fwhm_absoluter', 'center_absoluter', 'peakval_absoluter', 'gene'])
data_IF = pd.DataFrame(columns=['sample', 'group', 'rep', 'n', 'IF', 'gene'])
data_RNAFISH = pd.DataFrame(columns=['sample', 'group', 'rep', 'n', 'RNAFISH', 'gene'])
data_IF_mean = pd.DataFrame(columns=['sample', 'rep', 'n', 'IF', 'gene'])
data_RNAFISH_mean = pd.DataFrame(columns=['sample', 'rep', 'n', 'RNAFISH', 'gene'])

pd_gene = pd.read_csv("%s/genelist_koscreen1_1.txt" % master_folder, na_values=['.'], sep='\t')


def get_max_index(lst: list):
    return lst.index(max(lst))


def get_surrounding_points(lst: list, y):
    out = []
    for i in range(len(lst)):
        if (y < lst[i]) & (len(out) == 0):
            out.append([i-1, lst[i-1]])
            out.append([i, lst[i]])
        elif (y > lst[i]) & (len(out) == 2):
            out.append([i - 1, lst[i - 1]])
            out.append([i, lst[i]])
    return out


def find_x_from_y_and_line(point1: list, point2: list, y):
    return (y-point1[1])*(point2[0]-point1[0])/(point2[1]-point1[1]) + point1[0]


def get_fwhm_and_center(lst: list):
    halfpeak_value = 1+0.5*(max(lst)-1)
    surrounding_points = get_surrounding_points(lst, halfpeak_value)
    x1 = find_x_from_y_and_line(surrounding_points[0], surrounding_points[1], halfpeak_value)
    x2 = find_x_from_y_and_line(surrounding_points[2], surrounding_points[3], halfpeak_value)
    fwhm = x2-x1
    center = (x2+x1)/2
    return fwhm, center


for sample in samples:
    print(sample)
    gene = pd_gene[pd_gene['location'] == sample]['gene'].tolist()[0]
    for rep in reps:
        pd_IF = pd.read_csv("%s/%s/25_%s_IF_rep%s_red_green_group.txt" % (data_dir, sample, sample, rep), na_values=['.'], sep='\t')
        pd_RNAFISH = pd.read_csv("%s/%s/28_%s_RNAFISH_rep%s_red_green_group.txt" % (data_dir, sample, sample, rep),
                                 na_values=['.'], sep='\t')

        for hue in hueorder:
            data_IF.loc[len(data_IF.index)] = [sample, hue, rep, len(pd_IF[pd_IF['group'] == hue]), np.mean(pd_IF[pd_IF['group'] == hue]['IF']), gene]
            data_RNAFISH.loc[len(data_RNAFISH.index)] = [sample, hue, rep, len(pd_RNAFISH[pd_RNAFISH['group'] == hue]),
                                                         np.mean(pd_RNAFISH[pd_RNAFISH['group'] == hue]['RNAFISH']), gene]
        data_IF_temp = data_IF[(data_IF['sample'] == sample) & (data_IF['rep'] == rep)].copy().reset_index(drop=True)
        print(len(data_IF_temp))
        data_IF_mean.loc[len(data_IF_mean.index)] = [sample, rep, data_IF_temp['n'][1] / data_IF_temp['n'][0],
                                                     data_IF_temp['IF'][1]/data_IF_temp['IF'][0], gene]
        data_RNAFISH_temp = data_RNAFISH[(data_RNAFISH['sample'] == sample) & (data_RNAFISH['rep'] == rep)].copy().reset_index(drop=True)
        print(len(data_RNAFISH_temp))
        data_RNAFISH_mean.loc[len(data_RNAFISH_mean.index)] = [sample, rep, data_RNAFISH_temp['n'][1] / data_RNAFISH_temp['n'][0],
                                                               data_RNAFISH_temp['RNAFISH'][1]/data_RNAFISH_temp['RNAFISH'][0], gene]

    pd_cp = pd.read_csv("%s/%s/20_%s_copy_number_group.txt" % (data_dir, sample, sample), na_values=['.'], sep='\t')
    pd_cp['log10_DNAFISH_total_int_merge'] = np.log10(pd_cp['DNAFISH_total_int_merge'] + 1)
    if version == 1:
        pd_cluster = pd.read_csv("%s/%s/21_%s_cluster_group.txt" % (data_dir, sample, sample), na_values=['.'], sep='\t')
        df_relativer = pd.read_csv('%s/%s/23_%s_radial_summary_relativer.txt' % (data_dir, sample, sample), na_values=['.'],
                                   sep='\t')
        df_absoluter = pd.read_csv('%s/%s/23_%s_radial_summary_absoluter.txt' % (data_dir, sample, sample), na_values=['.'],
                                   sep='\t')
    else:
        pd_cluster = pd.read_csv("%s/%s/33_%s_cluster_group.txt" % (data_dir, sample, sample), na_values=['.'],
                                 sep='\t')
        df_relativer = pd.read_csv('%s/%s/35_%s_radial_summary_relativer.txt' % (data_dir, sample, sample),
                                   na_values=['.'],
                                   sep='\t')
        df_absoluter = pd.read_csv('%s/%s/35_%s_radial_summary_absoluter.txt' % (data_dir, sample, sample),
                                   na_values=['.'],
                                   sep='\t')

    for h in range(len(hueorder)):
        hue = hueorder[h]
        relativer_lst = df_relativer.iloc[h].tolist()[20:]
        absoluter_lst = df_absoluter.iloc[h].tolist()[20:]

        relativer_peak_value = max(relativer_lst)
        absoluter_peak_value = max(absoluter_lst)

        relativer_peak = get_max_index(relativer_lst) + 20
        absoluter_peak = get_max_index(absoluter_lst) + 20

        relativer_fwhm, relativer_center = get_fwhm_and_center(relativer_lst)
        absoluter_fwhm, absoluter_center = get_fwhm_and_center(absoluter_lst)

        relativer_center = relativer_center + 20
        absoluter_center = absoluter_center + 20

        data.loc[len(data.index)] = [sample, hue, len(pd_cluster[pd_cluster['group'] == hue]),
                                     np.mean(pd_cluster[pd_cluster['group'] == hue]['area_nuclear']),
                                     np.median(pd_cluster[pd_cluster['group'] == hue]['area_nuclear']),
                                     np.mean(pd_cp[pd_cp['group'] == hue]['log10_DNAFISH_total_int_merge']),
                                     np.median(pd_cp[pd_cp['group'] == hue]['log10_DNAFISH_total_int_merge']),
                                     len(pd_cluster[(pd_cluster['group'] == hue) & (pd_cluster['n_ecDNA'] != -1) & (pd_cluster['total_area_ratio_ecDNA'] < 0.5)]),
                                     np.mean(pd_cluster[(pd_cluster['group'] == hue) & (pd_cluster['n_ecDNA'] != -1) & (pd_cluster['total_area_ratio_ecDNA'] < 0.5)]['n_ecDNA']),
                                     np.median(pd_cluster[(pd_cluster['group'] == hue) & (pd_cluster['n_ecDNA'] != -1) & (pd_cluster['total_area_ratio_ecDNA'] < 0.5)]['n_ecDNA']),
                                     len(pd_cluster[(pd_cluster['group'] == hue) & (pd_cluster['total_area_ratio_ecDNA'] > 0.02) & (pd_cluster['n_ecDNA'] != -1) & (pd_cluster['total_area_ratio_ecDNA'] < 0.5)]),
                                     np.mean(pd_cluster[(pd_cluster['group'] == hue) & (pd_cluster['total_area_ratio_ecDNA'] > 0.02) & (pd_cluster['n_ecDNA'] != -1) & (pd_cluster['total_area_ratio_ecDNA'] < 0.5)]['per_AUC']),
                                     np.median(pd_cluster[(pd_cluster['group'] == hue) & (pd_cluster['total_area_ratio_ecDNA'] > 0.02) & (pd_cluster['n_ecDNA'] != -1) & (pd_cluster['total_area_ratio_ecDNA'] < 0.5)]['per_AUC']),
                                     relativer_peak, relativer_fwhm, relativer_center, relativer_peak_value,
                                     absoluter_peak, absoluter_fwhm, absoluter_center, absoluter_peak_value, gene]

    data_temp = data[data['sample'] == sample].copy().reset_index(drop=True)
    print(len(data_temp))
    data_mean.loc[len(data_mean.index)] = [sample, data_temp['mean_area_nuclear'][1] / data_temp['mean_area_nuclear'][0],
                                           data_temp['median_area_nuclear'][1] / data_temp['median_area_nuclear'][0],
                                           data_temp['mean_copy_number'][1]/data_temp['mean_copy_number'][0],
                                           data_temp['median_copy_number'][1] / data_temp['median_copy_number'][0],
                                           data_temp['mean_n_ecDNA'][1] / data_temp['mean_n_ecDNA'][0],
                                           data_temp['median_n_ecDNA'][1] / data_temp['median_n_ecDNA'][0],
                                           data_temp['mean_per_AUC'][1]/data_temp['mean_per_AUC'][0],
                                           data_temp['median_per_AUC'][1] / data_temp['median_per_AUC'][0],
                                           data_temp['peak_relativer'][1]/data_temp['peak_relativer'][0],
                                           data_temp['fwhm_relativer'][1]/data_temp['fwhm_relativer'][0],
                                           data_temp['center_relativer'][1]/data_temp['center_relativer'][0],
                                           data_temp['peakval_relativer'][1]/data_temp['peakval_relativer'][0],
                                           data_temp['peak_absoluter'][1] / data_temp['peak_absoluter'][0],
                                           data_temp['fwhm_absoluter'][1] / data_temp['fwhm_absoluter'][0],
                                           data_temp['center_absoluter'][1] / data_temp['center_absoluter'][0],
                                           data_temp['peakval_absoluter'][1] /data_temp['peakval_absoluter'][0], gene]

if version == 1:
    data.to_csv('%s/02_summary.txt' % output_dir, index=False, sep='\t')
    data_mean.to_csv('%s/02_summary_mean.txt' % output_dir, index=False, sep='\t')
else:
    data.to_csv('%s/02_summary_v2.txt' % output_dir, index=False, sep='\t')
    data_mean.to_csv('%s/02_summary_mean_v2.txt' % output_dir, index=False, sep='\t')
data_IF.to_csv('%s/02_summary_IF.txt' % output_dir, index=False, sep='\t')
data_IF_mean.to_csv('%s/02_summary_IF_mean.txt' % output_dir, index=False, sep='\t')
data_RNAFISH.to_csv('%s/02_summary_RNAFISH.txt' % output_dir, index=False, sep='\t')
data_RNAFISH_mean.to_csv('%s/02_summary_RNAFISH_mean.txt' % output_dir, index=False, sep='\t')
print("DONE!")

