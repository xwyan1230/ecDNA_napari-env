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
    print('i%s' % i)
    if i == 9:
        data = pd.read_csv('%s/%s/XY%s.txt' % (data_dir, folder, i+1), na_values=['.'], sep='\t')
    else:
        data = pd.read_csv('%s/%s/XY0%s.txt' % (data_dir, folder, i + 1), na_values=['.'], sep='\t')

    data['log10_emiRFP670'] = np.log10(data['emiRFP670'])

    hc = 17000
    print(np.mean(data[data['hoechst']>hc]['hoechst']))
    fig, ax = plt.subplots(figsize=(9, 6))
    # fig.subplots_adjust(right=0.8)
    sns.histplot(data=data, x='hoechst') # annot=True
    plt.axvline(hc, 0, 1000, c='r')
    # plt.legend(loc=(1.04, 0))
    # plt.savefig('%s/%s/hoechst_hist/XY%s.pdf' % (output_dir, folder, i+1))
    plt.show()

    lc = 2.95
    uc = 2.95
    print(len(data[data['hoechst'] > hc]))
    print(len(data[(data['hoechst'] > hc) & (data['log10_emiRFP670'] < lc)]))
    print(len(data[(data['hoechst'] > hc) & (data['log10_emiRFP670'] > uc)]))
    print(len(data[(data['hoechst'] > hc) & (data['log10_emiRFP670'] < lc)])/len(data[data['hoechst'] > hc]))
    print(len(data[(data['hoechst'] > hc) & (data['log10_emiRFP670'] > uc)]) / len(data[data['hoechst'] > hc]))

    fig, ax = plt.subplots(figsize=(9, 6))
    sns.histplot(data=data[data['hoechst'] > hc], x='log10_emiRFP670', bins=20)  # annot=True
    plt.xlim([2, 6])
    plt.axvline(lc, 0, 1000, c='r')
    plt.axvline(uc, 0, 1000, c='r')
    # plt.axvline(16000, 0, 1000, c='r')
    # plt.legend(loc=(1.04, 0))
    # plt.savefig('%s/%s/hoechst_hist/XY%s.pdf' % (output_dir, folder, i+1))
    plt.show()

