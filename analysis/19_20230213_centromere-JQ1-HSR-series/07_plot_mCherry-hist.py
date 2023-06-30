import skimage.io as skio
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from skimage.measure import label, regionprops

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230206_analysis_centromere-JQ1-DM-HSR/"
data_dir = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

sample = 'DM-JQ1_mix_mCh-Ctrl'
n_dilation = 4

df = pd.read_csv(("%s%s_n%s.txt" % (data_dir, sample, n_dilation)), na_values=['.'], sep='\t')

low = 6000
high = 10000

plt.subplots(figsize=(9, 6))
sns.histplot(data=df, x='mCherry_mean', bins=40)
plt.axvline(x=low, color='r', linestyle='--')
plt.axvline(x=high, color='r', linestyle='--')
plt.savefig('%s%s_hist_mCherry_mean.tiff' % (output_dir, sample))
plt.show()

print(len(df[df['mCherry_mean'] < low]))
print(len(df[df['mCherry_mean'] > high]))
print(len(df[(df['mCherry_mean'] >= low) & (df['mCherry_mean'] <= high)]))