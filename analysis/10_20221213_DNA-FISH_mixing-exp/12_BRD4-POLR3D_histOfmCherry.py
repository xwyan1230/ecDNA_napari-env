import skimage.io as skio
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from skimage.measure import label, regionprops

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221213_analysis_DNA-FISH_mixing_exp/20221121_H2B-series_POLR3D-BRD4-BRD1-DAPK2/"
data_dir = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

sample = 'H2B+BRD1'
df = pd.read_csv(("%s%s_ecDNA.txt" % (data_dir, sample)), na_values=['.'], sep='\t')

plt.subplots(figsize=(9, 6))
sns.histplot(data=df, x='mCherry_mean', bins=40)
plt.savefig('%shist_mCherry_mean_%s.tiff' % (output_dir, sample))
plt.show()

print(len(df))
print(len(df[df['mCherry_mean'] < 20000]))
print(len(df[df['mCherry_mean'] > 30000]))