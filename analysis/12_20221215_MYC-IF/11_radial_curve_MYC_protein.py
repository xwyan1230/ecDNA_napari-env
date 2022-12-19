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
data_dir = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

sample = 'MYC-new'
df = pd.read_csv(("%s%s_MYC.txt" % (data_dir, sample)), na_values=['.'], sep='\t')

hist_colors = [(0.90, 0.90, 0.90), (0.95, 0.50, 0.50),  (0.50, 0.90, 0.90)]
line_colors = [(0.30, 0.30, 0.30), (0.85, 0.35, 0.25),  (0.30, 0.70, 0.70)]

feature = ['radial_curve_nuclear', 'radial_curve_MYC', 'radial_curve_normalized']
for f in feature:
    df[f] = [dat.str_to_float(df[f][i]) for i in range(len(df))]

# radial curve
print("Plotting radial curve...")
x = np.arange(0.05, 1, 0.1)
x_label = 'relative r'

data = df
number_nuclear = len(data)

mean_curve1, ci_lower1, ci_higher1 = dat.mean_list(data['radial_curve_MYC'].tolist())
mean_curve2, ci_lower2, ci_higher2 = dat.mean_list(data['radial_curve_nuclear'].tolist())
mean_curve3, ci_lower3, ci_higher3 = dat.mean_list(data['radial_curve_normalized'].tolist())

plt.subplots(figsize=(12, 9))
for i in range(len(data)):
    plt.plot(x, data['radial_curve_MYC'][i], alpha=0.05, color=[line_colors[1][j] + 0.05 for j in range(len(line_colors[1]))])
for i in range(len(data)):
    plt.plot(x, data['radial_curve_nuclear'][i], alpha=0.05, color=[line_colors[0][j]+0.05 for j in range(len(line_colors[0]))])
plt.plot(x, mean_curve1, color=line_colors[1], label='%s, n=%s' % ('MYC IF', number_nuclear))
plt.plot(x, mean_curve2, color=line_colors[0], label='%s, n=%s' % ('DNA (hoechst stain)', number_nuclear))
plt.plot(x, ci_lower1, color=line_colors[1], linestyle='--', linewidth=0.5)
plt.plot(x, ci_higher1, color=line_colors[1], linestyle='--', linewidth=0.5)
plt.plot(x, ci_lower2, color=line_colors[0], linestyle='--', linewidth=0.5)
plt.plot(x, ci_higher2, color=line_colors[0], linestyle='--', linewidth=0.5)
plt.xlabel(x_label)
plt.ylim([0.5, 1.5])
plt.ylabel('radial_curve')
plt.legend()
plt.savefig('%sradial_curve_%s.pdf' % (output_dir, sample))
plt.show()

plt.subplots(figsize=(12, 9))
for i in range(len(data)):
    plt.plot(x, data['radial_curve_normalized'][i], alpha=0.05, color=[line_colors[1][j] + 0.05 for j in range(len(line_colors[1]))])
plt.plot(x, mean_curve3, color=line_colors[1], label='%s, n=%s' % ('MYC IF (normalized)', number_nuclear))
plt.plot(x, ci_lower3, color=line_colors[1], linestyle='--', linewidth=0.5)
plt.plot(x, ci_higher3, color=line_colors[1], linestyle='--', linewidth=0.5)
plt.axhline(y=1, color='black', linestyle='--')
plt.xlabel(x_label)
plt.ylim([0.5, 1.5])
plt.ylabel('radial_curve')
plt.legend()
plt.savefig('%sradial_curve_normalized_%s.pdf' % (output_dir, sample))
plt.show()