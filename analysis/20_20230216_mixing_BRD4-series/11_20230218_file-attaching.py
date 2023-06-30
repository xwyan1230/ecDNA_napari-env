import skimage.io as skio
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from skimage.measure import label, regionprops

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230216_analysis_BRD4-series/"
data_dir = "%sfigures/DNAFISH/" % master_folder
output_dir = "%sfigures/DNAFISH/" % master_folder

sample = 'DM-Ctrl_mix_mCh-Ctrl'
n_dilation = 4
total_batch = 2

data = pd.DataFrame()
for i in range(total_batch):
    df = pd.read_csv(("%s%s_%s_n%s_1.txt" % (data_dir, sample, i+1, n_dilation)), na_values=['.'], sep='\t')
    data = pd.concat([data, df], axis=0)

data.to_csv('%s%s_n%s_1.txt' % (output_dir, sample, n_dilation), index=False, sep='\t')
print("DONE!")
