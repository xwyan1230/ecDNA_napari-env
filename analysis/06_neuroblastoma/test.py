import matplotlib
import numpy as np

cmap = matplotlib.cm.get_cmap('Spectral')

rgba = [cmap(i) for i in np.arange(0, 1, 0.1)]
print(len(rgba))