import pandas as pd
import numpy as np

data = pd.DataFrame({'a': [1,2,3], 'b': [1,2,3]})
data_target = data.iloc[[1]]
data_rest = data.copy().drop(1).reset_index(drop=True)
data_current = pd.concat([data_target, data_rest], axis=0).reset_index(drop=True)

test=np.array([[1,2,3],[1,2,3],[1,2,3]])
print(test.max())