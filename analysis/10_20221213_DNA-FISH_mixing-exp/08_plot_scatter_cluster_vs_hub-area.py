import skimage.io as skio
import pandas as pd
import matplotlib.pyplot as plt
import shared.dataframe as dat
import seaborn as sns
from skimage.measure import label, regionprops

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221213_analysis_DNA-FISH_mixing_exp/20221017_mixing-test_mCherry-series_after-heating/"
data_dir = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

df = pd.read_csv(("%sDM_mix_DM-H2B-mCherry_ecDNA.txt" % data_dir), na_values=['.'], sep='\t')

hist_colors = [(0.95, 0.50, 0.50), (0.90, 0.90, 0.90), (0.50, 0.90, 0.90)]
line_colors = [(0.85, 0.35, 0.25), (0.30, 0.30, 0.30), (0.30, 0.70, 0.70)]

sample = []
for i in range(len(df)):
    if df['mCherry_mean'].tolist()[i] < 5000:
        sample.append('neg')
    elif df['mCherry_mean'].tolist()[i] > 10000:
        sample.append('pos')
    else:
        sample.append('NA')
df['sample'] = sample

df_sort = df[df['sample'].isin(['neg', 'pos'])].copy().reset_index(drop=True)
df_pos = df[df['sample'].isin(['pos'])].copy().reset_index(drop=True)
df_neg = df[df['sample'].isin(['neg'])].copy().reset_index(drop=True)

sns.set_palette(sns.color_palette(line_colors))
plt.subplots(figsize=(9,6))
sns.scatterplot(data=df_sort, x='total_area_ecDNA', y='dis_to_hub_area', hue='sample')
plt.savefig('%s/dis_to_hub_area_vs_total_area.pdf' % output_dir)
plt.show()
