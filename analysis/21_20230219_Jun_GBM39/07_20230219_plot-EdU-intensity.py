import skimage.io as skio
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from skimage.measure import label, regionprops
import shared.dataframe as dat

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230219_analysis_Jun_EGFR_RPAs33p_Edu/"
data_dir = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

n_dilation = 0

# samples
sample1 = 'gbm39ec con'
sample2 = 'gbm39hsr con'
figure_name = 'gbm39ec_vs_hsr'

hue_order = [sample1, sample2]
line_colors = [(0.30, 0.30, 0.30), (0.85, 0.35, 0.25)]  # black, red

# load data
df1 = pd.read_csv(("%s%s_n%s.txt" % (data_dir, sample1, n_dilation)), na_values=['.'], sep='\t')
df2 = pd.read_csv(("%s%s_n%s.txt" % (data_dir, sample2, n_dilation)), na_values=['.'], sep='\t')
df1['sample'] = [sample1] * len(df1)
df2['sample'] = [sample2] * len(df2)
print(len(df1))
df1 = df1[df1['circ_nuclear'] > 0.7].copy().reset_index(drop=True)
print(len(df1))
print(len(df2))
df2 = df2[df2['circ_nuclear'] > 0.7].copy().reset_index(drop=True)
print(len(df2))

df = pd.concat([df1, df2], axis=0).reset_index(drop=True)

cutoff = 30

"""sns.set_palette(sns.color_palette(line_colors))
plt.subplots(figsize=(9, 6))
sns.histplot(data=df, x='mean_int_EdU', bins=40, hue='sample')
plt.axvline(x=cutoff, color='r', linestyle='--')
plt.savefig('%s%s_hist_EdU_mean.tiff' % (output_dir, figure_name))
plt.show()"""

print(len(df[df['mean_int_EdU'] < cutoff]))
print(len(df[df['mean_int_EdU'] > cutoff]))

print(len(df1[df1['mean_int_EdU'] < cutoff]))
print(len(df1[df1['mean_int_EdU'] > cutoff]))

print(len(df2[df2['mean_int_EdU'] < cutoff]))
print(len(df2[df2['mean_int_EdU'] > cutoff]))