import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA as sklearnPCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
import scipy.cluster.hierarchy as shc
from pandas.plotting import parallel_coordinates
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
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


hit_lists = []
features = ['log2_mean_n_neg', 'log2_mean_n_pos']
mean_vals = [[-0.02337924154156307, -0.008967950850693634, -0.017017078281940476, -0.010755872968990585],
             [-0.04356061175326381, -0.029352313214499957, -0.04423469973669149, -0.029881275098349855]]
std_vals = [[0.17685021665897127, 0.12496598282510274, 0.12633541871010792, 0.14610200717457078],
            [0.2533806378651046, 0.18007797517772448, 0.2376410432445826, 0.22690028387274297]]
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
    df_flt = df_flt[df_flt['%s_hit_sum' % feature] >= 2].copy()
    print(len(df_flt))
    df_flt = df_flt.sort_values(by=['%s_hit_sum' % feature, '%s_sum' % feature], ascending=[False, True])
    df_flt.index = df_flt['treatment']
    hit_lists.append(df_flt['treatment'].tolist())
    df_flt = df_flt.drop(['index', 'treatment', 'seq', 'target', '%s_hit_sum' % feature, '%s_sum' % feature], axis=1)
    df_flt1 = df_flt.copy()
    df_flt1 = df_flt1.drop(['point5uM_24hr_%s_hit' % feature, 'point5uM_48hr_%s_hit' % feature, '5uM_24hr_%s_hit' % feature, '5uM_48hr_%s_hit' % feature], axis=1)
    df_flt2 = df_flt.copy()
    df_flt2 = df_flt2.drop(['point5uM_24hr_%s' % feature, 'point5uM_48hr_%s' % feature, '5uM_24hr_%s' % feature, '5uM_48hr_%s' % feature], axis=1)
    # pie chart
    """plt.subplots(figsize=(9, 7))
    y = np.array([10, 3, 46, 19, 152])
    plt.pie(y, colors=[rainboo_colors[9], rainboo_colors[8], rainboo_colors[7], rainboo_colors[6], rainboo_colors[5]])
    plt.savefig('%s/pie_hoechst_neg_hit.pdf' % output_dir)
    plt.show()

    plt.subplots(figsize=(9, 7))
    y = np.array([10, 5, 28, 28, 159])
    plt.pie(y, colors=rainboo_colors[:5])
    plt.savefig('%s/pie_hoechst_pos_hit.pdf' % output_dir)
    plt.show()"""

    # heat map
    fig, ax = plt.subplots(figsize=(6, 20))
    fig.subplots_adjust(left=0.4)
    ax1 = sns.heatmap(df_flt1, cbar=0, linewidths=0.2, vmax=6, vmin=-6, square=True,
                      cmap='coolwarm')  # yticklabels=False
    # plt.savefig('%s/heatmap_%s_log2_survival_seq_hit_log2.pdf' % (output_dir, feature))
    plt.show()

    norm = matplotlib.colors.Normalize(-1, 1)
    colors = [[norm(-1.0), "darkblue"],
              [norm(-0.6), "#d0d2d3"],
              [norm(0.6), "#d0d2d3"],
              [norm(1.0), rainboo_colors[0]]]

    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", colors)

    fig, ax = plt.subplots(figsize=(6, 20))
    fig.subplots_adjust(left=0.4)
    ax1 = sns.heatmap(df_flt2, cbar=0, linewidths=0.2, vmax=1, vmin=-1, square=True,
                      cmap=cmap)  # yticklabels=False
    # plt.savefig('%s/heatmap_%s_log2_survival_seq_hit.pdf' % (output_dir, feature))
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

# setA = hits[(hits['batch'] == conA[0]) & (hits['category'] == conA[1])]['hits'].tolist()[0]
# setB = hits[(hits['batch'] == conB[0]) & (hits['category'] == conB[1])]['hits'].tolist()[0]

setA = hit_lists[0]
setB = hit_lists[1]

print(set(setA) - set(setB))
print(set(setB) - set(setA))
print(set(setA) & set(setB))

plt.subplots(figsize=(9, 9))
venn2([set(setA), set(setB)], set_labels=['neg', 'pos'])
plt.savefig('%s/venn2_neg_vs_pos.pdf' % (output_dir))
plt.show()

gene_lst = list(set(setA) & set(setB))
pd_gene = pd.DataFrame()
pd_gene['treatment'] = gene_lst
pd_gene['target'] = ['DMSO' if pd_gene['treatment'][i] == 'DMSO'else target_df[target_df['compound name'] == pd_gene['treatment'][i]]['TargetGene'].tolist()[0] for i in range(len(pd_gene))]
print(pd_gene)

df_gene_total = gene_table(df_drop)
df_gene_DM = gene_table(pd_gene)
gene_percentage = []
for i in range(len(df_gene_DM)):
    gene_percentage.append(
        int(df_gene_DM['n'][i]) / df_gene_total[df_gene_total['gene'] == df_gene_DM['gene'][i]]['n'].tolist()[0])
df_gene_DM['percentage'] = gene_percentage
df_temp = df_gene_DM[~((df_gene_DM['n'] == 1)&(df_gene_DM['percentage']==1))]
df_temp = df_temp[df_temp['percentage']>0.3]
print(df_temp)