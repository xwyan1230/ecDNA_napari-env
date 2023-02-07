import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats

# INPUT PARAMETERS
# file info
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20230119_model_clustering/"
data_dir = "%s" % master_folder
output_dir = "%s" % master_folder

# x_rvs = pd.Series(poisson.rvs(1.2, size=1000, random_state=2))
# data = x_rvs.value_counts().sort_index().to_dict()

data = pd.read_csv('%sJun16.txt' % data_dir, na_values=['.'], sep='\t')
DM_copy = data['DM_copy'].tolist()
print(DM_copy)

"""dist = scipy.stats.binom
bounds = [(0, 300), (0.5, 0.5)]
res = scipy.stats.fit(dist, DM_copy, bounds)
print(res.params)
print(np.mean(DM_copy))

res.plot()
plt.show()"""

mu = np.mean(DM_copy)
std = np.std(DM_copy)
print(mu)
print(std)

plt.subplots(figsize=(12, 9))
sns.histplot(data=data, x='DM_copy', bins=15)
plt.savefig('%sDM_copy.pdf' % output_dir)
plt.show()



