import pandas as pd
import numpy as np

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240517_analysis_KOscreen1/"
data_dir = "%sprocessed/" % master_folder
output_dir = "%ssummary/" % master_folder

version = 2

pd_IF = pd.read_csv("%s/02_summary_IF_mean.txt" % output_dir, na_values=['.'], sep='\t')
pd_RNAFISH = pd.read_csv("%s/02_summary_RNAFISH_mean.txt" % output_dir, na_values=['.'], sep='\t')
if version == 1:
    pd_DNAFISH = pd.read_csv("%s/02_summary_mean.txt" % output_dir, na_values=['.'], sep='\t')
else:
    pd_DNAFISH = pd.read_csv("%s/02_summary_mean_v2.txt" % output_dir, na_values=['.'], sep='\t')
pd_n = pd.DataFrame()
pd_n['sample'] = pd_IF['sample'].tolist() + pd_RNAFISH['sample'].tolist()
pd_n['category'] = ['IF'] * len(pd_IF) + ['RNAFISH'] * len(pd_RNAFISH)
pd_n['rep'] = pd_IF['rep'].tolist() + pd_RNAFISH['rep'].tolist()
pd_n['n'] = pd_IF['n'].tolist() + pd_RNAFISH['n'].tolist()
pd_n['gene'] = pd_IF['gene'].tolist() + pd_RNAFISH['gene'].tolist()

samples = ['B%s' % (x+2) for x in range(10)] + ['C%s' % (x+2) for x in range(10)] + \
          ['D%s' % (x+2) for x in range(10)] + ['E%s' % (x+2) for x in range(10)] + \
          ['F%s' % (x+2) for x in range(10)] + ['G%s' % (x+2) for x in range(10)]

ctrl_AAVS = ['B6', 'C9', 'D8', 'E7', 'F2', 'G3']
ctrl_NTC = ['B10', 'C7', 'D11', 'E10', 'F9', 'G5']
ctrl_water = ['B2', 'C3', 'D4', 'E5', 'F6', 'G7']
ctrl = ctrl_AAVS + ctrl_NTC + ctrl_water

features = ['n', 'IF', 'RNAFISH', 'mean_area_nuclear', 'median_area_nuclear',
            'mean_copy_number', 'median_copy_number', 'mean_n_ecDNA', 'median_n_ecDNA',
            'mean_per_AUC', 'median_per_AUC', 'peak_relativer',
            'fwhm_relativer', 'center_relativer', 'peakval_relativer', 'peak_absoluter', 'fwhm_absoluter',
            'center_absoluter', 'peakval_absoluter']
pds = [pd_n, pd_IF, pd_RNAFISH] + [pd_DNAFISH]*16
feature_means = []
for i in range(len(features)):
    feature = features[i]
    pd_temp = pds[i]
    feature_mean = np.mean(pd_temp[pd_temp['sample'].isin(ctrl)][feature])
    feature_means.append(feature_mean)
    pd_temp['%s_log2fc' % feature] = [np.log2(pd_temp[feature][j] / feature_mean) for j in range(len(pd_temp))]

pd_IF.to_csv('%s/03_summary_IF_log2fc.txt' % output_dir, index=False, sep='\t')
pd_RNAFISH.to_csv('%s/03_summary_RNAFISH_log2fc.txt' % output_dir, index=False, sep='\t')
pd_DNAFISH.to_csv('%s/03_summary_log2fc.txt' % output_dir, index=False, sep='\t')
pd_n.to_csv('%s/03_summary_n_log2fc.txt' % output_dir, index=False, sep='\t')
pds = [pd_n, pd_IF, pd_RNAFISH] + [pd_DNAFISH]*16

pd_final = pd.DataFrame(columns=['sample'] + features)
for i in range(len(samples)):
    s = samples[i]
    features_log2fc = []
    for f in range(len(features)):
        feature = features[f]
        pd_temp = pds[f]
        pd_s = pd_temp[pd_temp['sample'] == s].copy().reset_index(drop=True)
        features_log2fc.append(np.mean(pd_s['%s_log2fc' % feature]))
    pd_final.loc[len(pd_final.index)] = [s] + features_log2fc

if version == 1:
    pd_final.to_csv('%s/03_summary_log2fc_final.txt' % output_dir, index=False, sep='\t')
else:
    pd_final.to_csv('%s/03_summary_log2fc_final_v2.txt' % output_dir, index=False, sep='\t')
print(feature_means)
print("DONE!")


