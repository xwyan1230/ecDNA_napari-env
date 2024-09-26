import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import shared.dataframe as dat
from scipy.optimize import curve_fit
import numpy as np
import seaborn as sns
import os
from venny4py.venny4py import *

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240422_analysis_Colo320_kinaseScreen/"
data_dir = "%sprocessed/" % master_folder
output_dir = "%sfigures/" % master_folder

conA = ['5uM_48hr', 'survival']
conB = ['5uM_48hr', 'survival_n']

hits = pd.read_csv('%s/hit.txt' % master_folder, na_values=['.'], sep='\t')
hits['hits'] = [dat.str_to_list(hits['hits'][i]) for i in range(len(hits))]

setA = hits[(hits['batch'] == conA[0]) & (hits['category'] == conA[1])]['hits'].tolist()[0]
setB = hits[(hits['batch'] == conB[0]) & (hits['category'] == conB[1])]['hits'].tolist()[0]

print(set(setA) - set(setB))
print(set(setB) - set(setA))

plt.subplots(figsize=(9, 9))
venn2([set(setA), set(setB)], set_labels=['%s_%s' % (conA[0], conA[1]), '%s_%s' % (conB[0], conB[1])])
plt.savefig('%s/venn2_%s_%s_vs_%s_%s.pdf' % (output_dir, conA[0], conA[1], conB[0], conB[1]))
plt.show()
