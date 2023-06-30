import skimage.io as skio
import pandas as pd
import matplotlib.pyplot as plt
import shared.dataframe as dat
import seaborn as sns
from shared.sinaplot import sinaplot
import numpy as np
from skimage.measure import label, regionprops

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230216_analysis_BRD4-series/"
data_dir = "%sfigures/DNAFISH/" % master_folder
output_dir = "%sfigures/DNAFISH/" % master_folder

# samples
sample1 = 'DM-Ctrl_mix_mCh-Ctrl'
figure_name = sample1
sample1_pos = 'DM H2B-mCherry'
sample1_neg = 'DM'

hue_order = [sample1_neg, sample1_pos]
line_colors = [(0.30, 0.30, 0.30), (0.85, 0.35, 0.25)]  # black, red
feature = 'bg'

# load data
df1 = pd.read_csv(("%s%s_background.txt" % (data_dir, sample1)), na_values=['.'], sep='\t')
# df1['bg'] = [float(i) for i in df1['bg'].tolist()]
df = df1

sns.set_palette(sns.color_palette(line_colors))
fig, ax = plt.subplots(figsize=(2*len(hue_order)+1, 9))
fig.subplots_adjust(left=0.2)
# sns.violinplot(data=df, x='sample', y=feature, order=hue_order)
sinaplot(data=df, x='sample', y=feature, order=hue_order, violin=False, scale='area')
plt.savefig('%s/%s_%s.pdf' % (output_dir, figure_name, feature))
plt.show()

for i in hue_order:
    print(i)
    temp = df[df['sample'] == i].copy().reset_index(drop=True)
    print(np.mean(temp[feature]))