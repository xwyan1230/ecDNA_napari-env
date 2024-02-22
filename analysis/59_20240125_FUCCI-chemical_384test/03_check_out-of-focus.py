import skimage.io as skio
import napari
import shared.display as dis
import matplotlib.pyplot as plt
from skimage.measure import regionprops, label
from skimage.morphology import binary_dilation, disk
import numpy as np
import shared.segmentation as seg
import shared.objects as obj
import seaborn as sns
import pandas as pd
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240125_analysis_FUCCI-chemical_384test/"
data_dir = "%sfigures/" % master_folder
output_dir = "%sfigures/" % master_folder

folder = '20x_3x3_2500'

for i in range(10):
    i=1
    print('i%s' % i)
    if i == 9:
        data = pd.read_csv('%s/%s/XY%s.txt' % (data_dir, folder, i+1), na_values=['.'], sep='\t')
    else:
        data = pd.read_csv('%s/%s/XY0%s.txt' % (data_dir, folder, i + 1), na_values=['.'], sep='\t')

    data['log10_emiRFP670'] = np.log10(data['emiRFP670'])

    hc = 17000
    for j in range(9):
        j=0
        temp = data[(data['hoechst'] > hc) & (data['fov'] == j)].copy().reset_index(drop=True)
        print(len(temp))
        print(np.mean(data[data['fov'] == j]['hoechst']))
        print(len(data[(data['fov']==j)&(data['hoechst']<hc)])/len(data[data['fov']==j]))

        lc = 2.95
        uc = 2.95
        print(len(temp[temp['log10_emiRFP670'] < lc]))
        print(len(temp[temp['log10_emiRFP670'] > uc]))
        print(len(temp[temp['log10_emiRFP670'] < lc])/len(temp))
        print(len(temp[temp['log10_emiRFP670'] > uc])/len(temp))
        fig, ax = plt.subplots(figsize=(9, 6))
        # fig.subplots_adjust(right=0.8)
        sns.histplot(data=data[data['fov']==j], x='hoechst', bins=50) # annot=True
        plt.axvline(hc, 0, 1000, c='r')
        # plt.legend(loc=(1.04, 0))
        # plt.savefig('%s/%s/hoechst_hist/XY%s.pdf' % (output_dir, folder, i+1))
        plt.show()


