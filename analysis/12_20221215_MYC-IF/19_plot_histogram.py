import skimage.io as skio
import pandas as pd
import matplotlib.pyplot as plt
import shared.dataframe as dat
import seaborn as sns
import matplotlib as mpl
import numpy as np
from matplotlib import cm
from skimage.measure import label, regionprops

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221215_analysis_MYC-IF/20221117_immunoFISH_acid/"
data_dir = "%sfigures/txt/03_radial_calibrated/" % master_folder
output_dir = "%sfigures/" % master_folder

sample = 'MYC-new'
n = 8
df = pd.read_csv(("%s%s_radial_n%s.txt" % (data_dir, sample, n)), na_values=['.'], sep='\t')

hist_colors = [(0.90, 0.90, 0.90), (0.95, 0.50, 0.50),  (0.50, 0.90, 0.90)]
line_colors = [(0.30, 0.30, 0.30), (0.85, 0.35, 0.25),  (0.30, 0.70, 0.70)]

feature = ['radial_curve_nuclear', 'radial_curve_DNAFISH', 'radial_curve_normalized']
for f in feature:
    df[f] = [dat.str_to_float(df[f][i]) for i in range(len(df))]

df['total_area_ecDNA_sqrt'] = np.sqrt(df['total_area_ecDNA'])

feature = 'MYC_mean'
plt.subplots(figsize=(12, 9))
sns.histplot(data=df, x=feature)
plt.savefig('%s%s_%s.pdf' % (output_dir, feature, sample))
plt.show()

