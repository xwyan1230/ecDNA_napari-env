import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA as sklearnPCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
import scipy.cluster.hierarchy as shc
from pandas.plotting import parallel_coordinates
import matplotlib.pyplot as plt
import matplotlib.cm as pcm
import matplotlib
import seaborn as sns
import shared.dataframe as dat
import shared.display as dis
import numpy as np
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240422_analysis_Colo320_kinaseScreen/"
data_dir = "%sprocessed/" % master_folder
output_dir = "%sfigures/" % master_folder

batches = ['point5uM_24hr', 'point5uM_48hr', '5uM_24hr', '5uM_48hr']
exclude = ['Alisertib', 'AZD2811', 'Barasertib', 'PF-03814735', 'BDP5290', 'GSK269962A', 'Y-33075', 'Axitinib', 'Midostaurin', 'MRT67307']
rainboo_colors = ['#9e4747', '#bc4d4a', '#ce6a6b', '#ecada3', '#f0d7c7', '#bfd3c3', '#669daa', '#3e7287', '#395b77', '#1e2f54']
target_df = pd.read_excel('%s/kinase_inhibitor_target.xlsx' % master_folder, na_values=['.'])
survival_seq = pd.read_csv('%s/average_survival_hoechst_scale_log2.txt' % (output_dir), na_values=['.'], sep='\t')

features = ['log2_mean_n_neg', 'log2_mean_n_pos', 'log2_ratio']
for feature in features:
    df = pd.DataFrame()
    for batch in batches:
        data = pd.read_csv('%s/%s/%s_average_update1.txt' % (data_dir, batch, batch), na_values=['.'], sep='\t')
        df['index'] = data.index
        df['treatment'] = data['treatment']
        df['%s_%s' % (batch, feature)] = data[feature]

    df['seq'] = survival_seq['seq']
    print(len(df))
    df_drop = df.copy()
    df_drop = df_drop[~df_drop['treatment'].isin(exclude)].copy().reset_index(drop=True)
    print(len(df_drop))
    df_drop = df_drop.sort_values(by='seq').reset_index(drop=True)
    pd_ctrl = df_drop.copy()
    df_drop.index = df_drop['treatment']
    df_drop = df_drop.drop(['index', 'treatment', 'seq'], axis=1)

    limit = 6
    if feature == 'log2_ratio':
        limit = 3
    # heat map
    fig, ax = plt.subplots(figsize=(6, 40))
    fig.subplots_adjust(left=0.4)
    ax1 = sns.heatmap(df_drop, cbar=0, linewidths=0.2, vmax=limit, vmin=-limit, square=True, cmap='coolwarm')  #  yticklabels=False
    plt.savefig('%s/heatmap_%s_log2_survival_seq.pdf' % (output_dir, feature))
    plt.show()

    """pd_ctrl['ctrl'] = [1 if pd_ctrl['treatment'][i] == 'DMSO'else 0 for i in range(len(pd_ctrl))]
    pd_ctrl = pd_ctrl.drop(['index', 'treatment', 'seq', 'point5uM_24hr_%s' % feature, 'point5uM_48hr_%s' % feature, '5uM_24hr_%s' % feature, '5uM_48hr_%s' % feature], axis=1)
    
    norm = matplotlib.colors.Normalize(-1,1)
    colors = [[norm(-1.0), "darkblue"],
              [norm(-0.6), "#d0d2d3"],
              [norm( 0.6), "#d0d2d3"],
              [norm( 1.0), rainboo_colors[6]]]
    
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", colors)
    
    fig, ax = plt.subplots(figsize=(6, 40))
    fig.subplots_adjust(left=0.4)
    ax1 = sns.heatmap(pd_ctrl, cbar=0, linewidths=0.2, vmax=1, vmin=-1, square=True, cmap=cmap, yticklabels=False)  #  yticklabels=False
    plt.savefig('%s/heatmap_pos-neg_ctrl_log2_survival_seq.pdf' % output_dir)
    plt.show()"""