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

conA = ['point5uM_24hr', 'survival_all']
conB = ['point5uM_48hr', 'survival_all']
conC = ['5uM_24hr', 'survival_all']
conD = ['5uM_48hr', 'survival_all']

hits = pd.read_csv('%s/hit.txt' % master_folder, na_values=['.'], sep='\t')
hits['hits'] = [dat.str_to_list(hits['hits'][i]) for i in range(len(hits))]

setA = hits[(hits['batch'] == conA[0]) & (hits['category'] == conA[1])]['hits'].tolist()[0]
setB = hits[(hits['batch'] == conB[0]) & (hits['category'] == conB[1])]['hits'].tolist()[0]
setC = hits[(hits['batch'] == conC[0]) & (hits['category'] == conC[1])]['hits'].tolist()[0]
setD = hits[(hits['batch'] == conD[0]) & (hits['category'] == conD[1])]['hits'].tolist()[0]

print(setA)
print(setB)

# dict of sets
sets = {
    conA[0]: set(setA),
    conB[0]: set(setB),
    conC[0]: set(setC),
    conD[0]: set(setD)}

# plt.subplots(figsize=(9, 9))
venny4py(sets=sets)
plt.savefig('%s/venn4_%s_%s_vs_%s_%s_vs_%s_%s_vs_%s_%s.pdf' % (output_dir, conA[0], conA[1], conB[0], conB[1], conC[0], conC[1], conD[0], conD[1]))
plt.show()
