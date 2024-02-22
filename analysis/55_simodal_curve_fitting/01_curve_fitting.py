import pandas as pd
import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt
import seaborn as sns

master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20231107_simodal_curve_fitting/"
samples = ['GBM39ec', 'GBM39hsr']


def ll4(x,b,c,d,e):
    '''This function is basically a copy of the LL.4 function from the R drc package with
     - b: hill slope
     - c: min response
     - d: max response
     - e: EC50'''
    return(c+(d-c)/(1+np.exp(b*(np.log(x)-np.log(e)))))


data = pd.DataFrame()
for sample in samples:
    temp = pd.read_csv('%s/%s_Wee1_inhibitor.txt' % (master_folder, sample), na_values=['.'], sep='\t')
    data_temp = pd.DataFrame({'compound': ['%s' % sample] * 48,
                         'dose': [0.015625, 0.03125, 0.0625, 0.125, 0.25, 0.5, 1, 2, 4, 8, 16, 32] * 4,
                         'logDose': list(np.arange(-6, 6, 1)) * 4})
    data_temp['response'] = temp['response']
    data = pd.concat([data, data_temp], axis=0)
print(len(data))

compoundData = data.groupby(['compound'])
fitData = []
for name,group in compoundData:
    fitCoefs, covMatrix = opt.curve_fit(ll4, group.dose, group.response)
    resids = group.response-group.dose.apply(lambda x: ll4(x,*fitCoefs))
    curFit = dict(zip(['b','c','d','e'],fitCoefs))
    curFit['compound']=name
    curFit['residuals']=sum(resids**2)
    fitData.append(curFit)
fitCompound = [item['compound'] for item in fitData]
fitTable = pd.DataFrame(fitData).set_index('compound')
print(fitTable)

refDose = np.linspace(-7, 7, 256)

plt.subplots(figsize=(9, 6))
sns.lmplot('logDose', 'response', data=data, hue='compound', hue_order=samples, fit_reg=False, height=6, aspect=1.25)
for fit in fitData:
    plt.plot([np.log2(i) for i in refDose],[ll4(i,*[fit[i] for i in ['b','c','d','e']]) for i in refDose])
    plt.savefig('%s/%s_and_%s_Wee1_inhibitor.pdf' % (master_folder, samples[0], samples[1]))
plt.show()