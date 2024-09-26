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
exclude = ['Alisertib', 'AZD2811', 'Barasertib', 'PF-03814735', 'BDP5290', 'GSK269962A', 'Y-33075', 'Axitinib',
           'Midostaurin', 'MRT67307']
rainboo_colors = ['#9e4747', '#bc4d4a', '#ce6a6b', '#ecada3', '#f0d7c7', '#bfd3c3', '#669daa', '#3e7287', '#395b77',
                  '#1e2f54']
target_df = pd.read_excel('%s/kinase_inhibitor_target.xlsx' % master_folder, na_values=['.'])
survival_seq = pd.read_csv('%s/average_survival_hoechst_scale_log2.txt' % (output_dir), na_values=['.'], sep='\t')


def gene_table(data):
    gene_lst = []
    for i in range(len(data)):
        gene_lst = gene_lst + list(str(data['target'][i]).split(', '))
    gene_lst = list(set(gene_lst))
    n_gene = []
    for i in range(len(gene_lst)):
        n_gene.append(len(data[data['target'].str.contains(gene_lst[i])]))
    df = pd.DataFrame()
    df['gene'] = gene_lst
    df['n'] = n_gene
    return df


features = ['log2_ratio']
mean_vals = [[-0.007709193112894106, -0.014924431761693634, -0.00923016335318805, -0.011552983080901559]]
std_vals = [[0.14974719048439, 0.20726761398568413, 0.16388695015427213, 0.1823688326368574]]
for k in range(len(features)):
    df = pd.DataFrame()
    feature = features[k]
    for b in range(4):
        batch = batches[b]
        data = pd.read_csv('%s/%s/%s_average_update1.txt' % (data_dir, batch, batch), na_values=['.'], sep='\t')
        df['index'] = data.index
        df['treatment'] = data['treatment']
        df['%s_%s' % (batch, feature)] = data[feature]
        df['%s_%s_hit' % (batch, feature)] = \
            [1 if (df['%s_%s' % (batch, feature)][i] > (mean_vals[k][b] + 3 * std_vals[k][b])) |
                  (df['%s_%s' % (batch, feature)][i] < (mean_vals[k][b] - 3 * std_vals[k][b])) else 0
             for i in range(len(df))]
    df['%s_sum' % feature] = df['point5uM_24hr_%s' % feature] + df['point5uM_48hr_%s' % feature] + df[
        '5uM_24hr_%s' % feature] + df['5uM_48hr_%s' % feature]
    df['%s_hit_sum' % feature] = df['point5uM_24hr_%s_hit' % feature] + df['point5uM_48hr_%s_hit' % feature] + df[
        '5uM_24hr_%s_hit' % feature] + df['5uM_48hr_%s_hit' % feature]
    df['target'] = ['DMSO' if df['treatment'][i] == 'DMSO' else
                    target_df[target_df['compound name'] == df['treatment'][i]]['TargetGene'].tolist()[0] for i in
                    range(len(df))]

    df['seq'] = survival_seq['seq']
    print(len(df))
    df_drop = df.copy()
    df_drop = df_drop[~df_drop['treatment'].isin(exclude)].copy().reset_index(drop=True)
    print(len(df_drop))
    df_drop = df_drop.sort_values(by='seq').reset_index(drop=True)
    # pd_ctrl = df_drop.copy()

    df_flt = df_drop.copy()
    # df_flt = df_flt[df_flt['%s_hit_sum' % feature] >= 2].copy()
    # df_flt = df_flt.sort_values(by=['%s_hit_sum' % feature, '%s_sum' % feature], ascending=[False, True])
    df_temp1 = df_flt[df_flt['%s_sum' % feature]>=0].copy().reset_index(drop=True)
    df_temp2 = df_flt[df_flt['%s_sum' % feature]<0].copy().reset_index(drop=True)
    df_temp1 = df_temp1.sort_values(by=['%s_hit_sum' % feature, '%s_sum' % feature], ascending=[False, False])
    df_temp2 = df_temp2.sort_values(by=['%s_hit_sum' % feature, '%s_sum' % feature], ascending=[True, False])
    df_flt = pd.concat([df_temp1, df_temp2], axis=0).reset_index(drop=True)
    df_flt_save = df_flt.copy()
    df_flt_save['seq1'] = df_flt_save.index
    df_flt_save = df_flt_save.sort_values(by='index')
    df_flt_save.to_csv('%s/log2_ratio_sort_index.txt' % (output_dir), index=False, sep='\t')

    print(df_flt[:33])
    df_gene_total = gene_table(df_flt)
    df_gene_DM = gene_table(df_flt[:33])
    gene_percentage = []
    for i in range(len(df_gene_DM)):
        gene_percentage.append(
            int(df_gene_DM['n'][i]) / df_gene_total[df_gene_total['gene'] == df_gene_DM['gene'][i]]['n'].tolist()[0])
    df_gene_DM['percentage'] = gene_percentage
    print(df_flt[252:])
    df_gene_HSR = gene_table(df_flt[252:].reset_index(drop=True))
    gene_percentage = []
    for i in range(len(df_gene_HSR)):
        gene_percentage.append(
            int(df_gene_HSR['n'][i]) / df_gene_total[df_gene_total['gene'] == df_gene_HSR['gene'][i]]['n'].tolist()[0])
    df_gene_HSR['percentage'] = gene_percentage
    df_temp = df_gene_HSR[~((df_gene_HSR['n'] == 1) & (df_gene_HSR['percentage'] == 1))]
    df_temp = df_temp[df_temp['percentage'] > 0.3]
    print(df_temp)

    df_flt.index = df_flt['treatment']
    df_flt = df_flt.drop(['index', 'treatment', 'seq', 'target', '%s_hit_sum' % feature, '%s_sum' % feature], axis=1)
    df_flt1 = df_flt.copy()
    df_flt1 = df_flt1.drop(['point5uM_24hr_%s_hit' % feature, 'point5uM_48hr_%s_hit' % feature, '5uM_24hr_%s_hit' % feature, '5uM_48hr_%s_hit' % feature], axis=1)
    df_flt2 = df_flt.copy()
    df_flt2 = df_flt2.drop(['point5uM_24hr_%s' % feature, 'point5uM_48hr_%s' % feature, '5uM_24hr_%s' % feature, '5uM_48hr_%s' % feature], axis=1)

    # heat map
    fig, ax = plt.subplots(figsize=(6, 15))
    fig.subplots_adjust(left=0.4)
    ax1 = sns.heatmap(df_flt1[252:], cbar=0, linewidths=0.2, vmax=3, vmin=-3, square=True,
                      cmap='coolwarm')  # yticklabels=False
    plt.savefig('%s/heatmap_%s_HSR.pdf' % (output_dir, feature))
    plt.show()

    norm = matplotlib.colors.Normalize(-1, 1)
    colors = [[norm(-1.0), "darkblue"],
              [norm(-0.6), "#d0d2d3"],
              [norm(0.6), "#d0d2d3"],
              [norm(1.0), rainboo_colors[0]]]

    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", colors)

    fig, ax = plt.subplots(figsize=(6, 15))
    fig.subplots_adjust(left=0.4)
    ax1 = sns.heatmap(df_flt2[252:], cbar=0, linewidths=0.2, vmax=1, vmin=-1, square=True,
                      cmap=cmap)  # yticklabels=False
    plt.savefig('%s/heatmap_%s_hit_HSR.pdf' % (output_dir, feature))
    plt.show()

    fig, ax = plt.subplots(figsize=(6, 15))
    fig.subplots_adjust(left=0.4)
    ax1 = sns.heatmap(df_flt1[:33], cbar=0, linewidths=0.2, vmax=3, vmin=-3, square=True,
                      cmap='coolwarm')  # yticklabels=False
    plt.savefig('%s/heatmap_%s_DM.pdf' % (output_dir, feature))
    plt.show()

    norm = matplotlib.colors.Normalize(-1, 1)
    colors = [[norm(-1.0), "darkblue"],
              [norm(-0.6), "#d0d2d3"],
              [norm(0.6), "#d0d2d3"],
              [norm(1.0), rainboo_colors[0]]]

    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", colors)

    fig, ax = plt.subplots(figsize=(6, 15))
    fig.subplots_adjust(left=0.4)
    ax1 = sns.heatmap(df_flt2[:33], cbar=0, linewidths=0.2, vmax=1, vmin=-1, square=True,
                      cmap=cmap)  # yticklabels=False
    plt.savefig('%s/heatmap_%s_hit_DM.pdf' % (output_dir, feature))
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