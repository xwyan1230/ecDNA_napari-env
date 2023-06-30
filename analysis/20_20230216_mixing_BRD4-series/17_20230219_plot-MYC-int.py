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
data_dir = "%sfigures/IF/" % master_folder
output_dir = "%sfigures/IF/" % master_folder

# samples
sample1 = 'DM-CTCF-KO'
control = 'DM-Ctrl'
antibody = 'newMYC'
figure_name = sample1

hue_order = [control, sample1]
line_colors = [(0.30, 0.30, 0.30), (0.85, 0.35, 0.25)]  # black, red
feature = 'mean_int_MYC'

# load data
df1 = pd.read_csv(("%s%s_%s_MYC.txt" % (data_dir, sample1, antibody)), na_values=['.'], sep='\t')
dfc = pd.read_csv(("%s%s_%s_MYC.txt" % (data_dir, control, antibody)), na_values=['.'], sep='\t')

df1['sample'] = [sample1] * len(df1)
dfc['sample'] = [control] * len(dfc)

df = pd.concat([df1, dfc], axis=0)

sns.set_palette(sns.color_palette(line_colors))
fig, ax = plt.subplots(figsize=(2*len(hue_order)+1, 9))
fig.subplots_adjust(left=0.2)
# sns.violinplot(data=df, x='sample', y=feature, order=hue_order)
sinaplot(data=df, x='sample', y=feature, order=hue_order, violin=False, scale='area')
plt.savefig('%s/%s_%s_%s.pdf' % (output_dir, figure_name, antibody, feature))
plt.show()

for i in hue_order:
    print(i)
    temp = df[df['sample'] == i].copy().reset_index(drop=True)
    print(np.mean(temp[feature]))