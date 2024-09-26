import nd2
import napari
import pandas as pd
import numpy as np
import shared.image as ima
import tifffile as tif
import matplotlib.pyplot as plt
import skimage.io as skio
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20240715_analysis_PCC_cellcycle/"
data_dir = '%sprocessed/' % master_folder
output_dir = '%sfigures/' % master_folder

G1 = pd.read_csv('%s/02_G1.txt' % data_dir, na_values=['.'], sep='\t')
M = pd.read_csv('%s/01_M.txt' % data_dir, na_values=['.'], sep='\t')
S = pd.read_csv('%s/03_S.txt' % data_dir, na_values=['.'], sep='\t')

G1['total'] = G1['DM_n']+G1['HSR_n']
M['total'] = M['DM_n']+M['HSR_n']
S['total'] = S['DM_n']+S['HSR_n']

plt.subplots(figsize=(9, 7))
plt.hist(M['total'], weights=np.ones(len(M)) / len(M), range=[0, 150], bins=40, color='#dd933c', edgecolor='w', alpha=0.6)
plt.hist(S['total'], weights=np.ones(len(S)) / len(S), range=[0, 150], bins=40, color='#669daa', edgecolor='w', alpha=0.6)
plt.hist(G1['total'], weights=np.ones(len(G1)) / len(G1), range=[0, 150], bins=40, color='#bc4d4a', edgecolor='w', alpha=0.6)
plt.xlim([0, 150])
plt.ylim([0, 0.2])
if not os.path.exists('%s/' % (output_dir)):
    os.makedirs('%s/' % (output_dir))
plt.savefig('%s/cellcycle.pdf' % (output_dir))
plt.show()

plt.subplots(figsize=(9, 7))
plt.hist(M['total'], weights=np.ones(len(M)) / len(M), range=[0, 150], bins=40, color='#dd933c', edgecolor='w', alpha=1)
plt.xlim([0, 150])
plt.ylim([0, 0.2])
if not os.path.exists('%s/' % (output_dir)):
    os.makedirs('%s/' % (output_dir))
plt.savefig('%s/cellcycle_M.pdf' % (output_dir))
plt.show()

plt.subplots(figsize=(9, 7))
plt.hist(S['total'], weights=np.ones(len(S)) / len(S), range=[0, 150], bins=40, color='#669daa', edgecolor='w', alpha=1)
plt.xlim([0, 150])
plt.ylim([0, 0.2])
if not os.path.exists('%s/' % (output_dir)):
    os.makedirs('%s/' % (output_dir))
plt.savefig('%s/cellcycle_S.pdf' % (output_dir))
plt.show()

plt.subplots(figsize=(9, 7))
plt.hist(G1['total'], weights=np.ones(len(G1)) / len(G1), range=[0, 150], bins=40, color='#bc4d4a', edgecolor='w', alpha=1)
plt.xlim([0, 150])
plt.ylim([0, 0.2])
if not os.path.exists('%s/' % (output_dir)):
    os.makedirs('%s/' % (output_dir))
plt.savefig('%s/cellcycle_G1.pdf' % (output_dir))
plt.show()