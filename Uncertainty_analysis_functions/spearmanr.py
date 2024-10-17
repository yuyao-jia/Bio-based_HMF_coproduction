# -*- coding: utf-8 -*-
"""
#Compute Spearman's rank correlation coefficient and p-value
"""
import pandas as pd
import scipy
from scipy import stats
from model_uncertainty_analysis import *
from datetime import date

today = date.today()
date_formatted = today.strftime("%m%d%Y")

total_parameters = len(model.parameters)
def get_spearmanr(simulations_results = 'results.xlsx',rho_file='rho_file.xlsx',p_file='p_file.xlsx',n_parameters = total_parameters):
    data=pd.read_excel(simulations_results,skiprows=[0,2],index_col=[0])
    # parameters = ['Lipid content in dry oilcane','Feedstock capacity','Oilcane unit price','Furfural unit price','HMF unit price',
    #           'Biodiesel unit price','Glycerine unit price','Hexane unit price',
    #             'Nanofiltration membrane lifetime','Nanofiltration membrane unit cost','Microbial lipid yield from glucose and xylose',
    #           'Loss of vegetative lipid after pretreatment & saccharification',
    #             'Retention of furfural in 1st nanofitlration','Retention of HMF in 1st nanofitlration','Retention of HMF in 2nd nanofitlration',]
    # metric = ['HMF MPSP','MFPP',
    #       'furfural yield','HMF yield','biodiesel yield','crude glycerol yield',
    #       'furfural purity','hmf purity','Overall vegetative lipid recovery%',
    #       'Total capital investment','Annual operating cost','Annual material cost','Annual Heating demand',
    #       'Annual cooling demand','Annual electricity usage','NPV','IRR','Net GWP','Economic HMF GWP','Energy HMF GWP']
    # data.columns = ['Element',
    #             'Lipid content in dry oilcane','Feedstock capacity','Oilcane unit price','Furfural unit price','HMF unit price',
    #             'Biodiesel unit price','Glycerine unit price','Hexane unit price',
    #             'Nanofiltration membrane lifetime','Nanofiltration membrane unit cost','Microbial lipid yield from glucose and xylose',
    #             'Loss of vegetative lipid after pretreatment & saccharification',
    #             'Retention of furfural in 1st nanofitlration','Retention of HMF in 1st nanofitlration','Retention of HMF in 2nd nanofitlration',
    #             'HMF MPSP','MFPP',
    #             'furfural yield','HMF yield','biodiesel yield','crude glycerol yield',
    #             'furfural purity','hmf purity','Overall vegetative lipid recovery%',
    #             'Total capital investment','Annual operating cost','Annual material cost','Annual Heating demand',
    #             'Annual cooling demand','Annual electricity usage','NPV','IRR','Net GWP','Economic HMF GWP','Energy HMF GWP']
    col_index = data.columns
    parameters = col_index[0:n_parameters] # starting index is included, end index is not included
    metric = col_index[n_parameters:]
    rho = []
    p_value = []
    values_empty = {}
    rho_values = pd.DataFrame(values_empty)
    p_values = pd.DataFrame(values_empty)
    for i in metric:
        met = data[i].tolist()
        for j in parameters:
            par = data[j].tolist()
            res = stats.spearmanr(par,met)
            rho.append(res.correlation)
            p_value.append(res.pvalue)
        rho_values[i] = rho
        p_values[i] = p_value
        rho = []
        p_value = []
    rho_values.index = parameters
    p_values.index = parameters
    rho_values.to_excel(rho_file)
    p_values.to_excel(p_file)

    #quantile values
    q_values = pd.DataFrame(values_empty)
    q = [0.05,0.25,0.5,0.75,0.95]
    res_data = data[metric]
    for i in range(5):
        qi = res_data.quantile(q[i]).values
        q_values.insert(i,column=f"{q[i]}",value=qi)
    q_values.index = metric
    return rho_values, p_values

# rho, p = get_spearmanr(simulations_results = 'results.xlsx')

def get_samples_from_excel(sample_file = 'samples.xlsx'):
    data = pd.read_excel(sample_file,usecols='B:U')
    samples_all = data.values.tolist()

# path = '/Users/yuyaojia/Documents/Fermentation/Yuyao_TEA_LCA/Bioindustrial-Park/biorefineries/HMF/parameter_tests/test1_500'
path = '/tmp/pycharm_project_525/parameter_tests/test1_500'
data=pd.read_excel(f'{path}/06172024_N500_results.xlsx',skiprows=[0,2],index_col=[0])
cols = data.columns
parameters = cols[0:total_parameters]  # starting index is included, end index is not included
metric = cols[total_parameters:]

failed = pd.DataFrame(columns=cols)
for i in data.index:
    if np.isnan(data.loc[i]).any():
        failed.loc[i] = data.loc[i]
# failed.to_excel('failed_samples.xlsx')

failed_samples = failed[parameters].values
model.load_samples(failed_samples)
model.evaluate(notify=1)

results_failed = model.table.copy()
results_failed.index = failed.index

for i in results_failed.index:
    if np.isnan(results_failed.loc[i]).any() == False:
        failed = failed.drop(i)
        data.loc[i] = results_failed.loc[i].to_list()
    else:
        print(f"Failed to evaluate {i}")

#test if there is any NaN remaining in the data
# np.isnan(data.values).any()