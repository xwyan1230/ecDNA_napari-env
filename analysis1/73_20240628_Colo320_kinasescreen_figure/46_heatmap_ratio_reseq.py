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
# exclude_index = [78, 85, 134, 154, 188]
exclude_index = []
rainboo_colors = ['#9e4747', '#bc4d4a', '#ce6a6b', '#ecada3', '#f0d7c7', '#bfd3c3', '#669daa', '#3e7287', '#395b77', '#1e2f54']
target_df = pd.read_excel('%s/kinase_inhibitor_target.xlsx' % master_folder, na_values=['.'])
survival_seq = pd.read_csv('%s/average_survival_hoechst_scale_log2.txt' % (output_dir), na_values=['.'], sep='\t')

df = pd.DataFrame()
features = ['mean_n', 'log2_ratio', 'log2_mean_hoechst', 'log2_mean_n', 'log2_mean_n_neg', 'log2_mean_n_pos']
for feature in features:
    for batch in batches:
        data = pd.read_csv('%s/%s/%s_average_update1.txt' % (data_dir, batch, batch), na_values=['.'], sep='\t')
        df['index'] = data.index
        df['treatment'] = data['treatment']
        df['%s_%s' % (batch, feature)] = data[feature]

df_drop = df.copy()
df_drop = df_drop.drop(exclude_index).reset_index(drop=True)
print(len(df_drop))
"""for batch in batches:
    df_drop = df_drop[df_drop['%s_mean_n' % batch] > 1000].reset_index(drop=True)"""
print(len(df_drop))

df_label = pd.DataFrame()
df_label['ori_index'] = df_drop['index']
df_label['treatment'] = df_drop['treatment']
df_label['index'] = df_label.index
print(df_label)

df_num = pd.DataFrame()
df_num['log2_ratio'] = df_drop['point5uM_24hr_log2_ratio'].tolist() + df_drop['point5uM_48hr_log2_ratio'].tolist() + df_drop['5uM_24hr_log2_ratio'].tolist() + df_drop['5uM_48hr_log2_ratio'].tolist()
scaler = StandardScaler()
df_num_scale = pd.DataFrame(scaler.fit_transform(df_num))
df_num_scale.columns = df_num.columns

feature = 'log2_ratio'
n_len = len(df_drop)
df_num_cluster = pd.DataFrame()
for i in range(4):
    batch = batches[i]
    df_num_cluster['%s_%s' % (batch, feature)] = df_num_scale.iloc[i*n_len:(i+1)*n_len][feature].reset_index(drop=True)
print(df_num_cluster)

plt.figure(figsize=(20, 3))
clusters = shc.linkage(df_num_cluster, method='ward', metric="euclidean")
R = shc.dendrogram(Z=clusters)
plt.savefig('%s/ahc_ratio_noflt.pdf' % output_dir)
plt.show()
nodes = R['ivl']

df_label_sort = pd.DataFrame(columns=df_label.columns)
for i in range(len(df_label)):
    df_label_sort.loc[len(df_label_sort.index)] = df_label.iloc[int(dat.list_invert(nodes)[i])]
df_label_sort['seq'] = df_label_sort.index
"""seq = []
for i in range(len(df_label_sort)):
    if df_label_sort['seq1'][i] >= (n_len-17):
        seq.append(i-(n_len-17))
    elif df_label_sort['seq1'][i] <= 35:
        seq.append(52-i)
    else:
        seq.append(i+17)
df_label_sort['seq'] = seq"""
df_label_sort = df_label_sort.sort_values(by='index').reset_index(drop=True)
df_label_sort['target'] = ['DMSO' if df_label_sort['treatment'][i] == 'DMSO'else target_df[target_df['compound name'] == df_label_sort['treatment'][i]]['TargetGene'].tolist()[0] for i in range(len(df_label_sort))]
print(df_label_sort)
df_label_sort.to_csv('%s/label_ratio_noflt.txt' % (output_dir), index=False, sep='\t')

df_label_sort1 = df_label_sort.copy()
df_label_sort1 = df_label_sort1.sort_values(by='seq')
df_label_sort1.to_csv('%s/label_ratio_noflt_sort.txt' % (output_dir), index=False, sep='\t')

"""df_cluster = df_cluster_log2.copy()
df_cluster['index'] = df_label_sort['index']
df_cluster['treatment'] = df_label_sort['treatment']
df_cluster['seq'] = df_label_sort['seq']
df_cluster['target'] = df_label_sort['target']
df_cluster.to_csv('%s/average_survival_hoechst_scale_log2.txt' % (output_dir), index=False, sep='\t')

df_cluster_sort = df_cluster.copy()
df_cluster_sort = df_cluster_sort.sort_values(by='seq').reset_index(drop=True)
df_cluster_sort.to_csv('%s/average_survival_hoechst_scale_sort_log2.txt' % (output_dir), index=False, sep='\t')"""

features = ['log2_ratio', 'log2_mean_hoechst', 'log2_mean_n', 'log2_mean_n_neg', 'log2_mean_n_pos']
mean_vals = [[-0.007709193112894106, -0.014924431761693634, -0.00923016335318805, -0.011552983080901559],
             [-6.125465398936461e-05, -4.115266751604684e-05, -3.1213083382731086e-05, -4.0784426923211045e-05],
             [-0.029512062154734148, -0.01025642796528373, -0.024416822422841142, -0.013235743526684733],
             [-0.02337924154156307, -0.008967950850693634, -0.017017078281940476, -0.010755872968990585],
             [-0.04356061175326381, -0.029352313214499957, -0.04423469973669149, -0.029881275098349855]]
std_vals = [[0.14974719048439, 0.20726761398568413, 0.16388695015427213, 0.1823688326368574],
            [0.00834228398227882, 0.00819877007504514, 0.0066035016188801645, 0.008849233347603581],
            [0.20196940317110101, 0.11260042275606479, 0.16174595859458207, 0.161922769066963],
            [0.17685021665897127, 0.12496598282510274, 0.12633541871010792, 0.14610200717457078],
            [0.2533806378651046, 0.18007797517772448, 0.2376410432445826, 0.22690028387274297]]

for k in range(len(features)):
    feature = features[k]
    df_num_cluster_sort = pd.DataFrame()
    for b in range(4):
        batch = batches[b]
        df_num_cluster_sort['%s_%s' % (batch, feature)] = df_drop['%s_%s' % (batch, feature)]
        """df_num_cluster_sort['%s_%s_hit' % (batch, feature)] = \
            [1 if (df_num_cluster_sort['%s_%s' % (batch, feature)][i]>(mean_vals[k][b]+3*std_vals[k][b])) |
                  (df_num_cluster_sort['%s_%s' % (batch, feature)][i]<(mean_vals[k][b]-3*std_vals[k][b])) else 0
             for i in range(len(df_num_cluster_sort))]
        df_num_cluster_sort = df_num_cluster_sort.drop('%s_%s' % (batch, feature), axis=1)"""
    df_num_cluster_sort['treatment'] = df_label_sort['treatment'].tolist()
    df_num_cluster_sort['seq'] = df_label_sort['seq'].tolist()
    df_num_cluster_sort = df_num_cluster_sort.sort_values(by='seq')
    df_num_cluster_sort.index = df_num_cluster_sort['treatment']
    df_num_cluster_sort = df_num_cluster_sort.drop(['seq', 'treatment'], axis=1)
    print(df_num_cluster_sort)

    # limit = 1
    if feature == 'log2_mean_hoechst':
        limit = 0.2
    elif feature == 'log2_ratio':
        limit = 2
    else:
        limit = 3
    # heat map
    fig, ax = plt.subplots(figsize=(6, 40))
    fig.subplots_adjust(left=0.4)
    ax1 = sns.heatmap(df_num_cluster_sort, cbar=0, linewidths=0.2, vmax=limit, vmin=-limit, square=True, cmap='coolwarm')  #  yticklabels=False
    plt.savefig('%s/heatmap_ratio_%s_noflt.pdf' % (output_dir, feature))
    plt.show()

    """fig, ax = plt.subplots(figsize=(6, 10))
    fig.subplots_adjust(left=0.4)
    ax1 = sns.heatmap(df_num_cluster_sort[0:53], cbar=0, linewidths=0.2, vmax=limit, vmin=-limit, square=True, cmap='coolwarm')  #  yticklabels=False
    plt.savefig('%s/heatmap_ratio_enlarge_%s.pdf' % (output_dir, feature))
    plt.show()"""