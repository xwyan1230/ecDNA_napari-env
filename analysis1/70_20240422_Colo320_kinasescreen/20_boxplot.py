import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import shared.dataframe as dat
from scipy.optimize import curve_fit
import numpy as np
import seaborn as sns
import os
from venny4py.venny4py import *

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240422_analysis_Colo320_kinaseScreen/"
data_dir = "%sprocessed/" % master_folder
output_dir = "%sfigures/" % master_folder

conA = ['5uM_48hr', 'survival_all']
batch = conA[0]
hits = pd.read_csv('%s/hit.txt' % master_folder, na_values=['.'], sep='\t')
hits['hits'] = [dat.str_to_list(hits['hits'][i]) for i in range(len(hits))]
setA = hits[(hits['batch'] == conA[0]) & (hits['category'] == conA[1])]['hits'].tolist()[0]

data = pd.read_csv('%s/%s/%s_summary.txt' % (data_dir, batch, batch), na_values=['.'], sep='\t')
data_cc = pd.read_csv('%s/%s/%s_summary_cc.txt' % (data_dir, batch, batch), na_values=['.'], sep='\t')
data_cc = data_cc.drop(columns=['screen', 'group', 'cell', 'treatment', 'target'])
data = pd.concat([data, data_cc], axis=1)
data['label'] = ['%s(%s)' % (data['treatment'][i], data['target'][i][:15]) for i in range(len(data))]
data['identity'] = ['hit' if data['label'][i] in setA else 'non-hit' for i in range(len(data))]
print(len(data))
data_flt = data[data['treatment']!='DMSO']
data_flt = data_flt[((data['mean_n_filtered']>=200) & (data['log2fc_cytokinesis']>=-1) & (data['log2fc_cytokinesis']<=1))]
print(len(data_flt))

features = ['cc_score', 'log2fc_per_G1', 'log2fc_per_G1S', 'log2fc_per_S', 'log2fc_per_G2M']
for feature in features:
    fig, ax = plt.subplots(figsize=(3, 5))
    fig.subplots_adjust(left=0.4)
    sns.boxplot(data=data_flt, x='identity', y=feature)
    # sns.swarmplot(data=data_flt, x='identity', y=feature, s=2, color=(255/255, 140/255, 0/255))
    if feature == 'cc_score':
        plt.ylim([-1,15])
    else:
        plt.ylim([-3,3])
    if not os.path.exists('%s/%s/' % (output_dir, batch)):
        os.makedirs('%s/%s/' % (output_dir, batch))
    plt.savefig('%s/%s/%s_%s_%s.pdf' % (output_dir, batch, batch, conA[1], feature))
    plt.show()