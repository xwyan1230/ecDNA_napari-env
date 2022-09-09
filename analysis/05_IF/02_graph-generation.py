import pandas as pd
import numpy as np
import shared.dataframe as dat
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA as sklearnPCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
from pandas.plotting import parallel_coordinates
import shared.display as dis
import random
import shared.math as mat
import seaborn as sns
import os

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220922_EdU_test/"
sample_lst = ['CNTRL', '10uM_1hr', '2uM_ON', '10uM_ON']
save_folder = master_folder

data = pd.DataFrame()

for sample in sample_lst:
