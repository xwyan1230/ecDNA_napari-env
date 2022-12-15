import skimage.io as skio
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from skimage.measure import label, regionprops

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20221213_analysis_DNA-FISH_mixing_exp/20221017_mixing-test_mCherry-series_after-heating/"
data_dir = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

df_Ctrl = pd.read_csv(("%sDM.txt" % data_dir), na_values=['.'], sep='\t')
df_mCherry = pd.read_csv(("%sDM-H2B-mCherry.txt" % data_dir), na_values=['.'], sep='\t')
df_mix = pd.read_csv(("%sDM_mix_DM-H2B-mCherry_ecDNA.txt" % data_dir), na_values=['.'], sep='\t')

df_Ctrl['sample'] = ['Ctrl'] * len(df_Ctrl)
df_mCherry['sample'] = ['H2B-mCherry'] * len(df_mCherry)
df_mix['sample'] = ['Ctrl_mix_H2B-mCherry'] * len(df_mix)

df = pd.concat([df_Ctrl, df_mCherry, df_mix], axis=0).reset_index(drop=True)

plt.subplots(figsize=(9, 6))
sns.histplot(data=df, x='mCherry_mean', hue='sample', bins=40)
# plt.savefig('%shist_mCherry_mean.tiff' % output_dir)
plt.show()

print(len(df_Ctrl))
print(len(df_Ctrl[df_Ctrl['mCherry_mean'] >10000]))
print(len(df_mCherry))
print(len(df_mCherry[df_mCherry['mCherry_mean'] >10000]))
print(len(df_mix))
print(len(df_mix[df_mix['mCherry_mean'] >10000]))

print(len(df_Ctrl))
print(len(df_Ctrl[df_Ctrl['mCherry_mean'] <5000]))
print(len(df_mCherry))
print(len(df_mCherry[df_mCherry['mCherry_mean'] <5000]))
print(len(df_mix))
print(len(df_mix[df_mix['mCherry_mean'] <5000]))