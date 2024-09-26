import pandas as pd
import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt
import seaborn as sns

master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20231107_simodal_curve_fitting/20240912/"
samples = ['Colo320DM', 'Colo320HSR']
chemical = 'Ponatinib'
mode = 'viability'


def ll4(x,b,c,d,e):
    '''This function is basically a copy of the LL.4 function from the R drc package with
     - b: hill slope
     - c: min response
     - d: max response
     - e: EC50'''
    return(c+(d-c)/(1+np.exp(b*(np.log(x)-np.log(e)))))


data = pd.DataFrame()
for sample in samples:
    temp = pd.read_csv('%s/%s_%s_%s.txt' % (master_folder, sample, chemical, mode), na_values=['.'], sep='\t')
    data_temp = pd.DataFrame({'compound': ['%s' % sample] * 40,
                         'dose': [0.03125, 0.0625, 0.125, 0.25, 0.5, 1, 2, 4, 8, 16] * 4,
                         'logDose': list(np.arange(-5, 5, 1)) * 4})
    data_temp['response'] = temp['response']
    data = pd.concat([data, data_temp], axis=0).reset_index(drop=True)
print(len(data))
print(data)

refDose = np.linspace(-7, 7, 256)

sns.lmplot(data=data, x='logDose', y='response', hue='compound', hue_order=samples, fit_reg=False, height=6, aspect=1.25)
plt.savefig('%s/%s_and_%s_%s_%s_wofit.pdf' % (master_folder, samples[0], samples[1], chemical, mode))
plt.show()

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

data['response_ori'] = data['response']
if mode == 'viability':
    data['response'] = [data['response_ori'][i]/float(fitTable['d'][0]) if data['compound'][i] == samples[0] else data['response_ori'][i]/float(fitTable['d'][1]) for i in range(len(data))]
else:
    data['response'] = [data['response_ori'][i] / float(fitTable['c'][0]) if data['compound'][i] == samples[0] else
                        data['response_ori'][i] / float(fitTable['c'][1]) for i in range(len(data))]

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

"""sns.lmplot(data=data, x='logDose', y='response', hue='compound', hue_order=samples, fit_reg=False, height=6, aspect=1.25)
for fit in fitData:
    plt.plot([np.log2(i) for i in refDose],[ll4(i,*[fit[i] for i in ['b','c','d','e']]) for i in refDose])
    plt.savefig('%s/%s_and_%s_%s_%s.pdf' % (master_folder, samples[0], samples[1], chemical, mode))
plt.show()"""

data['response_ori'] = data['response']
if mode == 'viability':
    data['response'] = [data['response_ori'][i]/float(fitTable['d'][0]) if data['compound'][i] == samples[0] else data['response_ori'][i]/float(fitTable['d'][1]) for i in range(len(data))]
else:
    data['response'] = [data['response_ori'][i] / float(fitTable['c'][0]) if data['compound'][i] == samples[0] else
                        data['response_ori'][i] / float(fitTable['c'][1]) for i in range(len(data))]

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

sns.lmplot(data=data, x='logDose', y='response', hue='compound', hue_order=samples, fit_reg=False, height=6, aspect=1.25)
for fit in fitData:
    plt.plot([np.log2(i) for i in refDose],[ll4(i,*[fit[i] for i in ['b','c','d','e']]) for i in refDose])
    plt.savefig('%s/%s_and_%s_%s_%s.pdf' % (master_folder, samples[0], samples[1], chemical, mode))
plt.show()