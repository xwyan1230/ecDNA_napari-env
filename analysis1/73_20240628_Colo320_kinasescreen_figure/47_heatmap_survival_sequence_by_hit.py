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
rainboo_colors = ['#9e4747', '#bc4d4a', '#ce6a6b', '#ecada3', '#f0d7c7', '#bfd3c3', '#669daa', '#3e7287', '#395b77', '#1e2f54']
target_df = pd.read_excel('%s/kinase_inhibitor_target.xlsx' % master_folder, na_values=['.'])

df = pd.DataFrame()
features = ['log2_mean_hoechst']
mean_vals = [[-6.125465398936461e-05, -4.115266751604684e-05, -3.1213083382731086e-05, -4.0784426923211045e-05]]
std_vals = [[0.00834228398227882, 0.00819877007504514, 0.0066035016188801645, 0.008849233347603581]]
for k in range(len(features)):
    feature = features[k]
    for b in range(4):
        batch = batches[b]
        data = pd.read_csv('%s/%s/%s_average_update1.txt' % (data_dir, batch, batch), na_values=['.'], sep='\t')
        df['index'] = data.index
        df['treatment'] = data['treatment']
        df['%s_%s' % (batch, feature)] = data[feature]
        df['%s_%s_hit' % (batch, feature)] = \
            [1 if (df['%s_%s' % (batch, feature)][i]>(mean_vals[k][b]+3*std_vals[k][b])) |
                  (df['%s_%s' % (batch, feature)][i]<(mean_vals[k][b]-3*std_vals[k][b])) else 0
             for i in range(len(df))]
    df['%s_sum' % feature] = df['point5uM_24hr_%s' % feature] + df['point5uM_48hr_%s' % feature] + df['5uM_24hr_%s' % feature] + df['5uM_48hr_%s' % feature]
    df['%s_hit_sum' % feature] = df['point5uM_24hr_%s_hit' % feature] + df['point5uM_48hr_%s_hit' % feature] + df['5uM_24hr_%s_hit' % feature] + df['5uM_48hr_%s_hit' % feature]
    df['target'] = ['DMSO' if df['treatment'][i] == 'DMSO' else target_df[target_df['compound name'] == df['treatment'][i]]['TargetGene'].tolist()[0] for i in range(len(df))]

for feature in features:
    df_heatmap = pd.DataFrame()
    for batch in batches:
        df_heatmap['treatment'] = df['treatment']
        df_heatmap['%s_%s' % (batch, feature)] = df['%s_%s' % (batch, feature)]
    df_heatmap['%s_sum' % feature] = df['%s_sum' % feature]
    df_heatmap['%s_hit_sum' % feature] = df['%s_hit_sum' % feature]
    df_heatmap['target'] = ['DMSO' if df['treatment'][i] == 'DMSO' else
                            target_df[target_df['compound name'] == df_heatmap['treatment'][i]]['TargetGene'].tolist()[
                                0] for i in range(len(df_heatmap))]
    df_heatmap = df_heatmap[
        (df_heatmap['%s_hit_sum' % feature] >= 1) & (df_heatmap['%s_sum' % feature] < 0)].reset_index(drop=True)
    df_heatmap = df_heatmap.sort_values(by=['%s_hit_sum' % feature, '%s_sum' % feature], ascending=[False, True])

    df_heatmap.to_csv('%s/label_hoechst_hit_by_hit.txt' % (output_dir), index=False, sep='\t')
    df_heatmap.index = df_heatmap['treatment']
    df_heatmap = df_heatmap.drop(['treatment', 'target', '%s_sum' % feature, '%s_hit_sum' % feature], axis=1)

    # heat map
    fig, ax = plt.subplots(figsize=(6, 15))
    fig.subplots_adjust(left=0.4)
    ax1 = sns.heatmap(df_heatmap, cbar=0, linewidths=0.2, vmax=0.2, vmin=-0.2, square=True, cmap='coolwarm')  #  yticklabels=False
    plt.savefig('%s/heatmap_hoechst_by_hit.pdf' % output_dir)
    plt.show()

    """fig, ax = plt.subplots(figsize=(6, 10))
    fig.subplots_adjust(left=0.4)
    ax1 = sns.heatmap(df_num_cluster_sort[0:53], linewidths=0.2, vmax=0.2, vmin=-0.2, square=True, cmap='coolwarm')  #  yticklabels=False
    plt.savefig('%s/heatmap_rescale_survival_hoechst_enlarge_log2.pdf' % output_dir)
    plt.show()"""

norm = matplotlib.colors.Normalize(-1,1)
colors = [[norm(-1.0), "darkblue"],
          [norm(-0.6), "#d0d2d3"],
          [norm( 0.6), "#d0d2d3"],
          [norm( 1.0), rainboo_colors[0]]]

cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", colors)

for feature in features:
    df_heatmap = pd.DataFrame()
    for batch in batches:
        df_heatmap['treatment'] = df['treatment']
        df_heatmap['%s_%s_hit' % (batch, feature)] = df['%s_%s_hit' % (batch, feature)]

    df_heatmap['%s_sum' % feature] = df['%s_sum' % feature]
    df_heatmap['%s_hit_sum' % feature] = df['%s_hit_sum' % feature]
    df_heatmap = df_heatmap[
        (df_heatmap['%s_hit_sum' % feature] >= 1) & (df_heatmap['%s_sum' % feature] < 0)].reset_index(drop=True)
    print(len(df_heatmap))
    df_heatmap = df_heatmap.sort_values(by=['%s_hit_sum' % feature, '%s_sum' % feature], ascending=[False, True])
    df_heatmap.index = df_heatmap['treatment']
    df_heatmap = df_heatmap.drop(['treatment', '%s_sum' % feature, '%s_hit_sum' % feature], axis=1)

    # heat map
    fig, ax = plt.subplots(figsize=(6, 15))
    fig.subplots_adjust(left=0.4)
    ax1 = sns.heatmap(df_heatmap, cbar=0, linewidths=0.2, vmax=1, vmin=-1, square=True, cmap=cmap)  #  yticklabels=False
    plt.savefig('%s/heatmap_hoechst_hit_by_hit.pdf' % output_dir)
    plt.show()

# pie chart
plt.subplots(figsize=(9, 7))
y = np.array([10, 10, 35, 20, 165])
plt.pie(y, colors=rainboo_colors[:5])
plt.savefig('%s/pie_hoechst_hit.pdf' % output_dir)
plt.show()


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


for feature in features:
    df_heatmap = pd.DataFrame()
    df_heatmap['treatment'] = df['treatment']
    df_heatmap['target'] = df['target']
    df_heatmap['%s_sum' % feature] = df['%s_sum' % feature]
    df_heatmap['%s_hit_sum' % feature] = df['%s_hit_sum' % feature]
    df_gene_total = gene_table(df_heatmap)
    print(df_gene_total)
    df_heatmap = df_heatmap[
        (df_heatmap['%s_hit_sum' % feature] >= 2) & (df_heatmap['%s_sum' % feature] < 0)].reset_index(drop=True)
    df_gene_hit = gene_table(df_heatmap)
    gene_percentage = []
    for i in range(len(df_gene_hit)):
        gene_percentage.append(int(df_gene_hit['n'][i])/df_gene_total[df_gene_total['gene'] == df_gene_hit['gene'][i]]['n'].tolist()[0])
    df_gene_hit['percentage'] = gene_percentage
    print(len(df_gene_hit))
    df_test = df_gene_hit[df_gene_hit['percentage']>=0.5]
    print(len(df_test))
    print(df_test)
    """df_test = df_gene_hit[df_gene_hit['n'] == 1]
    print(len(df_test))
    print(df_test)
    df_test = df_gene_hit[(df_gene_hit['n'] == 1) & (df_gene_hit['percentage'] == 1)]
    print(len(df_test))
    print(df_test)
    df_gene_hit = df_gene_hit[((df_gene_hit['n'] > 1) & (df_gene_hit['percentage']>0.5))]
    df_gene_hit = df_gene_hit.sort_values(by=['percentage', 'n'], ascending=[False, False]).reset_index(drop=True)
    df_gene_hit.index = df_gene_hit['gene']
    print(len(df_gene_hit))
    print(df_gene_hit)"""
    """df_heatmap = df_heatmap.sort_values(by=['%s_hit_sum' % feature, '%s_sum' % feature], ascending=[False, True]).reset_index(drop=True)
    
    # gene_lst = ['PLK1', 'MTOR', 'PTK2',  'PIK3CA', 'CHK1', 'FLT3', 'ALK', 'CDK2', 'CDK4', 'CDK7', 'CDK9', 'CDK12']
    gene_lst = []
    for i in range(len(df_heatmap)):
        gene_lst = gene_lst + list(str(df_heatmap['target'][i]).split(', '))
    gene_lst = list(set(gene_lst))
    print(gene_lst)
    # single_chemical_lst = ['CRT 0066101', 'OSU-T315', 'MRT199665', 'NCB-0846', 'Adavosertib', 'PF-3758309', 'CCT241533', 'URMC-099', 'YM-201636', 'BSJ-04-122', 'KU-60019', 'XL413']
    pd_genes = pd.DataFrame()
    pd_genes['treatment'] = df_heatmap['treatment']
    pd_genes['target'] = df_heatmap['target']

    for i in range(len(gene_lst)):
        gene = gene_lst[i]
        n_gene_lst = []
        for k in range(len(pd_genes)):
            if gene in str(pd_genes['target'][k]):
                n_gene_lst.append(1)
            else:
                n_gene_lst.append(0)
        pd_genes[gene] = n_gene_lst

    pd_genes.index = pd_genes['treatment']
    pd_genes = pd_genes.drop(['treatment', 'target'], axis=1)

    fig, ax = plt.subplots(figsize=(20, 10))
    fig.subplots_adjust(left=0.4)
    ax1 = sns.heatmap(pd_genes, cbar=0, linewidths=0.2, vmax=1, vmin=-1, square=True, cmap='coolwarm')  #  yticklabels=False
    plt.savefig('%s/heatmap_hoechst_hit_gene.pdf' % output_dir)
    plt.show()
"""