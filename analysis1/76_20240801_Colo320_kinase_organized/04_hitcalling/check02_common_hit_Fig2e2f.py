# Fig2e, 2f

import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import pandas as pd
import utilities as uti
import seaborn as sns
import matplotlib

# INPUT PARAMETERS
data_dir = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240801_analysis/summary/"
data_dir1 = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240801_analysis/notes/"
output_dir = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240801_analysis/summary/"

# PARAMETERS
rainboo_colors = ['#9e4747', '#bc4d4a', '#ce6a6b', '#ecada3', '#f0d7c7', '#bfd3c3', '#669daa', '#3e7287', '#395b77', '#1e2f54']
data_neg = pd.read_csv('%s/04_others/01_survival_n_neg.txt' % data_dir, na_values=['.'], sep='\t')
data_pos = pd.read_csv('%s/04_others/01_survival_n_pos.txt' % data_dir, na_values=['.'], sep='\t')
target_df = pd.read_csv('%s/targets.txt' % data_dir1, na_values=['.'], sep='\t')
n_target_limit = 5

setA = data_neg[(data_neg['hit_sum'] <= -2) & (data_neg['treatment'] != 'DMSO')]['treatment'].tolist()
setB = data_pos[(data_pos['hit_sum'] <= -2) & (data_pos['treatment'] != 'DMSO')]['treatment'].tolist()

plt.subplots(figsize=(9, 9))
venn2([set(setA), set(setB)], set_labels=['neg', 'pos'])
# plt.savefig('%s/venn2_neg_vs_pos.pdf' % (output_dir))
plt.show()

"""print(set(setA) - set(setB))
print(set(setB) - set(setA))
print(set(setA) & set(setB))"""
common_hits = list(set(setA) & set(setB))
data, data_n, data_percentage = uti.get_gene_stats(common_hits, target_df, n_target_limit)
print(data)

norm = matplotlib.colors.Normalize(-1, 1)
colors = [[norm(-1.0), "#d0d2d3"],
          [norm(0), "#d0d2d3"],
          [norm(1.0), rainboo_colors[0]]]

cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", colors)

limit = 4
fig, ax = plt.subplots(figsize=(6, 10))
fig.subplots_adjust(left=0.4)
ax1 = sns.heatmap(data_n, cbar=0, linewidths=0.2, vmax=limit, vmin=-limit, square=True,
                  cmap=cmap)  # yticklabels=False
# plt.savefig('%s/heatmap_%s_log2_survival_seq_hit.pdf' % (output_dir, feature))
plt.show()

fig, ax = plt.subplots(figsize=(6, 10))
fig.subplots_adjust(left=0.4)
ax1 = sns.heatmap(data_percentage, cbar=0, linewidths=0.2, vmax=1, vmin=-1, square=True,
                  cmap=cmap)  # yticklabels=False
# plt.savefig('%s/heatmap_%s_log2_survival_seq_hit.pdf' % (output_dir, feature))
plt.show()