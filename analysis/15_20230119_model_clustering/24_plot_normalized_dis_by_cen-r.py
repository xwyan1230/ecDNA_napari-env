import skimage.io as skio
import pandas as pd
import matplotlib.pyplot as plt
import shared.dataframe as dat
import seaborn as sns
import matplotlib as mpl
import numpy as np
import math
from matplotlib import cm
from skimage.measure import label, regionprops

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230119_model_clustering/"
data_dir = "%stxt/dataset3/" % master_folder
output_dir = "%sfigures1/" % master_folder

ac = 5000

seps = [0, 10, 20, 30, 40, 50, 60, 70]
sep_feature = 'cen_r'
cmap = mpl.cm.get_cmap('Spectral')
x = np.arange(0, 1, 1/len(seps))
line_color = [cmap(i) for i in x]

data = pd.DataFrame()
for j in seps:
    df = pd.read_csv(("%s%s/%s_5_cluster.txt" % (data_dir, j, ac)), na_values=['.'], sep='\t')
    feature = ['radial_curve_DNAFISH']
    for f in feature:
        df[f] = [dat.str_to_float(df[f][i]) for i in range(len(df))]
    data = pd.concat([data, df], axis=0)

data['r'] = np.sqrt(data['area_nuclear']/math.pi)
data['total_area_ecDNA_sqrt'] = np.sqrt(data['total_area_ecDNA']/math.pi)
data['total_area_ecDNA_sqrt_normalized'] = data['total_area_ecDNA_sqrt']/data['r']
data['dis_to_hub_area_v2_normalized'] = data['dis_to_hub_area_v2']/data['r']

plt.subplots(figsize=(12, 9))
sns.violinplot(data=data, x='cen_r', y='dis_to_hub_area_v2_normalized')
plt.savefig('%sac%s_dis_to_hub_area_v2_normalized.pdf' % (output_dir, ac))
plt.show()