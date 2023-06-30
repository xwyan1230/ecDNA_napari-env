import skimage.io as skio
import pandas as pd
import matplotlib.pyplot as plt
from shared.sinaplot import sinaplot
import shared.dataframe as dat
import seaborn as sns
import numpy as np
from skimage.measure import label, regionprops

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230609_analysis_DM_KOscreen_48hr/"
data_dir = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

n_dilation = 4

# samples
sample = 'C2'
hue_order = ['GFP', 'mCherry']

df = pd.read_csv('%s/%s_n4_simplified.txt' % (data_dir, sample), na_values=['.'], sep='\t')
line_colors = [(0.30, 0.30, 0.30), (0.85, 0.35, 0.25)]  # black, red
df_sample = df[df['group'].isin(hue_order)].copy().reset_index(drop=True)
df = df_sample
feature = ['area_ecDNA_individual', 'area_ratio_ecDNA_individual', 'mean_int_ecDNA_individual']
for f in feature:
    df[f] = [dat.str_to_float(df[f][i]) for i in range(len(df))]

data = pd.DataFrame(columns=['group', 'area_ecDNA_individual', 'area_ratio_ecDNA_individual', 'mean_int_ecDNA_individual'])
for i in range(len(df)):
    for j in range(len(df['area_ecDNA_individual'][i])):
        data.loc[len(data.index)] = [df['group'][i], df['area_ecDNA_individual'][i][j], df['area_ratio_ecDNA_individual'][i][j], df['mean_int_ecDNA_individual'][i][j]]

for f in feature:
    fig, ax = plt.subplots(figsize=(2*len(hue_order)+1, 9))
    fig.subplots_adjust(left=0.2)
    sinaplot(data=data, x='group', y=f, order=hue_order, violin=False, scale='area', point_size=2)
    # plt.ylim([0, 0.04])
    plt.savefig('%s/%s_%s.pdf' % (output_dir, sample, f))
    plt.show()

data_GFP = data[data['group'] == 'GFP'].copy().reset_index(drop=True)
data_mCherry = data[data['group'] == 'mCherry'].copy().reset_index(drop=True)

"""f = 'area_ecDNA_individual'
plt.subplots(figsize=(6, 4))
weights1 = np.ones_like(data_GFP[f]) / len(data_GFP)
weights2 = np.ones_like(data_mCherry[f]) / len(data_mCherry)
plt.hist([data_GFP[f], data_mCherry[f]], weights=[weights1, weights2], color=line_colors,
         edgecolor=(0.2, 0.2, 0.2), label=hue_order, bins=20, range=[20, 100])
plt.xlabel(f)
plt.ylabel('Probability')
plt.legend()
plt.savefig('%s/%s_%s_hist.pdf' % (output_dir, sample, f))
plt.close()"""

plt.subplots(figsize=(9, 6))
sns.scatterplot(data=data_GFP, y='mean_int_ecDNA_individual', x='area_ecDNA_individual', color=line_colors[0], s=5, alpha=1)
sns.scatterplot(data=data_mCherry, y='mean_int_ecDNA_individual', x='area_ecDNA_individual', color=line_colors[1], s=5, alpha=1)
plt.legend()
plt.savefig('%s/%s_int_vs_area_scatter.pdf' % (output_dir, sample))
plt.show()

print(np.mean(data_GFP['mean_int_ecDNA_individual']))
print(np.mean(data_mCherry['mean_int_ecDNA_individual']))