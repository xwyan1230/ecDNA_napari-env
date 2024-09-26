import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
import seaborn as sns
import collections
from scipy.stats import ks_2samp, ttest_ind
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240422_analysis_Colo320_kinaseScreen/"
data_dir = "%sprocessed/" % master_folder
output_dir = "%sprocessed/" % master_folder

target_df = pd.read_excel('%s/kinase_inhibitor_target.xlsx' % master_folder, na_values=['.'])
targets = [str(target_df['TargetGene'][i]).split(', ') for i in range(len(target_df))]
targets_temp = []
for i in range(len(targets)):
    for j in range(len(targets[i])):
        targets_temp.append(targets[i][j])
targets_set = list(set(targets_temp))
targets_set.remove('CDK')
targets_set.remove('pan-TRK')
targets_set.remove('nan')
targets_duplicated = [x for x, count in collections.Counter(targets_temp).items() if count > 1]
targets_duplicated.remove('nan')
print(len(targets_set))
print(targets_set)
print(len(targets_duplicated))
print(targets_duplicated)

batch = '5uM_24hr'
type = 'pos_vs_neg_reverse'

data = pd.read_csv('%s/%s/%s_summary.txt' % (data_dir, batch, batch), na_values=['.'], sep='\t')
data_cc = pd.read_csv('%s/%s/%s_summary_cc.txt' % (data_dir, batch, batch), na_values=['.'], sep='\t')
data_cc = data_cc.drop(columns=['screen', 'group', 'cell', 'treatment', 'target'])
data = pd.concat([data, data_cc], axis=1)
data['label'] = ['%s(%s)' % (data['treatment'][i], data['target'][i][:15]) for i in range(len(data))]
data['divide'] = data['mean_fov_hoechst']/data['mean_n_filtered']
mean_divide = np.mean(data['divide'])
# data_sort = data[(data['log2fc_fov_hoechst']<-1) & (data['log2fc_n_filtered']<-1)].copy().reset_index(drop=True)
data_sort = data[(data['mean_n_filtered'] >= 200) & (data['divide'] <= 1.6*mean_divide) & (data['divide'] >= 0.4*mean_divide) & (data['log2fc_pos_vs_neg'] < -0.7)].copy().reset_index(drop=True)
print(data_sort['label'])
print(len(data_sort))
genes = [str(data_sort['target'][i]).split(', ') for i in range(len(data_sort))]
genes_temp = []
for i in range(len(genes)):
    for j in range(len(genes[i])):
        genes_temp.append(genes[i][j])
genes_set = list(set(genes_temp))
test = ['CDK', 'pan-TRK', 'nan']
for i in test:
    if i in genes_set:
        genes_set.remove(i)
print(len(genes_set))
print(genes_set)
genes_duplicated = []
for i in genes_set:
    if i in targets_duplicated:
        genes_duplicated.append(i)
print(len(genes_duplicated))
print(genes_duplicated)

output_list = genes_set
output_pd = pd.DataFrame({'gene': output_list})
output_pd.to_csv('%s/%s_%s_gene_set_list_%s.txt' % (output_dir, batch, type, len(genes_set)), index=False, sep='\t', header=None)




