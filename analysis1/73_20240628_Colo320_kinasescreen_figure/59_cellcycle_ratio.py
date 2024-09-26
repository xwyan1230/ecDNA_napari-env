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
exclude_index = [154, 188]
rainboo_colors = ['#9e4747', '#bc4d4a', '#ce6a6b', '#ecada3', '#f0d7c7', '#bfd3c3', '#669daa', '#3e7287', '#395b77',
                  '#1e2f54']
exclude = ['Alisertib', 'AZD2811', 'Barasertib', 'PF-03814735', 'BDP5290', 'GSK269962A', 'Y-33075', 'Axitinib',
           'Midostaurin', 'MRT67307']
target_df = pd.read_excel('%s/kinase_inhibitor_target.xlsx' % master_folder, na_values=['.'])
ratio_seq = pd.read_csv('%s/log2_ratio_sort_index.txt' % (output_dir), na_values=['.'], sep='\t')

features = ['mean_delta_n_pos_G1', 'mean_delta_n_pos_G1S', 'mean_delta_n_pos_S', 'mean_delta_n_pos_G2M', 'mean_delta_n_pos_G2MG1']
for feature in features:
    df = pd.DataFrame()
    for batch in batches:
        data = pd.read_csv('%s/%s/%s_cellcycle_average_update1.txt' % (data_dir, batch, batch), na_values=['.'], sep='\t')
        df['index'] = data.index
        df['treatment'] = data['treatment']
        df['%s_%s' % (batch, feature)] = data[feature]

    print(len(df))
    df_drop = df.copy()
    df_drop = df_drop[~df_drop['treatment'].isin(exclude)].copy().reset_index(drop=True)
    print(len(df_drop))
    df_drop['seq'] = ratio_seq['seq1'].tolist()
    df_drop = df_drop[~df_drop['index'].isin(exclude_index)].reset_index(drop=True)
    print(len(df_drop))
    df_drop = df_drop.sort_values(by='seq').reset_index(drop=True)
    pd_ctrl = df_drop.copy()
    df_drop.index = df_drop['treatment']
    df_drop = df_drop.drop(['index', 'treatment', 'seq'], axis=1)

    limit = 0.5