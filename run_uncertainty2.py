import biosteam as bst
import numpy as np
import chaospy
from chaospy import distributions as shape
import matplotlib.pyplot as plt
from datetime import date
from model_uncertainty_analysis import *
import matplotlib.pyplot as plt
from scipy import stats
import pandas as pd
import argparse
from tqdm import tqdm
# parser = argparse.ArgumentParser()
# parser.add_argument('job_number', type=int)
# JOB_ID = parser.parse_args().job_number
# N_samples = parser.parse_args().job_number # 2000 samples are enough
N_samples = 2000

today = date.today()
date_formatted = today.strftime("%m%d%Y")
# date_formatted = '06132024'
# results_path = "/home/rodrigo/olive_results"
results_path='N2000_current_06242024/08212024_N2000'

# data = pd.read_excel('05282024_baseline_parameters.xlsx',usecols="B:U")
# samples = data.values
# model.load_samples(np.array(([samples[1]])))
# model.evaluate()

# #Generate samples & save them to an excel file
# N_samples = parser.parse_args().job_number # 1000 samples are enough?
# model.get_distribution_summary()['Uniform'].to_excel(f'{results_path}/{date_formatted}_N{N_samples}_parameter_distribution.xlsx')
# baseline_metrics.to_excel(f'{results_path}/{date_formatted}_baseline_metrics_verification.xlsx')


# rule = 'L'  # For Latin-Hypercube sampling
# np.random.seed(1234)  # For consistent results
# samples = model.sample(N_samples, rule)
# pd.DataFrame(samples).to_excel(f'{results_path}/{date_formatted}_N{N_samples}_samples.xlsx')

# data = pd.read_excel(f'06132024_N500_results_failed_1.xlsx',usecols='B:X',skiprows=[0,2])
# samples = data.values [59:]
# model.load_samples(samples)

# RUNS_PER_JOB = 1000
# samples_this_job = samples[RUNS_PER_JOB * JOB_ID : RUNS_PER_JOB * (JOB_ID + 1)]

# try:
#     model.evaluate(notify=50)
#     model.table.to_excel(f"{results_path}/{date_formatted}_N{N_samples}_results.xlsx",)
# except:
#     model.table.to_excel(f"{results_path}/{date_formatted}_N{N_samples}_results_failed.xlsx")

# model.evaluate(notify=50)
# df_rho, df_p = model.spearman_r()
# results = model.table
# results.to_excel(f"{results_path}/{date_formatted}_N{N_samples}_results.xlsx",)
# df_rho.to_excel(f"{results_path}/{date_formatted}_N{N_samples}_rho.xlsx",)
# data = model.table.copy()

# data = pd.read_excel(f'{results_path}/08212024_N2000_results.xlsx',skiprows=[0,2],index_col=[0])
data = pd.read_excel(f'{results_path}/08222024_N2000_results.xlsx',index_col=[0])
total_parameters = len(model.parameters)
cols = data.columns
parameters = cols[0:total_parameters]  # starting index is included, end index is not included
metric = cols[total_parameters:]
# results = data.copy(deep=True)
# def calculate_quantiles():
#     q_values = pd.DataFrame(index=metric)
#     q = [0.05,0.25,0.5,0.75,0.95]
#     res_data = data[metric]
#     for i in range(5):
#         qi = res_data.quantile(q[i]).values
#         q_values.insert(i,column=f"{q[i]}",value=qi)
#     return q_values
#
# q_values = calculate_quantiles()
# q_values.to_excel(f"{results_path}/{date_formatted}_N{N_samples}_quantiles.xlsx",)
#
# def run_failed_samples():
#     failed = pd.DataFrame(columns=cols)
#     for i in data.index:
#         if np.isnan(data.loc[i]).any():
#             failed.loc[i] = data.loc[i]
#     failed_samples = failed[parameters].values
#     model.load_samples(failed_samples)
#     model.evaluate(notify=1)
#
#     results_failed = model.table.copy()
#     results_failed.index = failed.index
#     for i in results_failed.index:
#         if np.isnan(results_failed.loc[i]).any() == False:
#             data.loc[i] = results_failed.loc[i].to_list()
#         else:
#             print(f"Failed to evaluate {i}")
#
# def calculate_spearmanrho():
#     rho = []
#     p_value = []
#     rho_values = pd.DataFrame(index=parameters, columns=metric)
#     p_values = pd.DataFrame(index=parameters, columns=metric)
#     for i in metric:
#         met = data[i].tolist()
#         for j in parameters:
#             par = data[j].tolist()
#             res = stats.spearmanr(par, met,nan_policy='omit')
#             rho.append(res.correlation)
#             p_value.append(res.pvalue)
#         rho_values[i] = rho
#         p_values[i] = p_value
#         rho = []
#         p_value = []
#     return rho_values, p_values
#
# calculated_rho_values, calculated_p_values = calculate_spearmanrho()
# calculated_rho_values.to_excel(f"{results_path}/{date_formatted}_N{N_samples}_rho.xlsx",)
#
# if np.isnan(data.values).any():
#     print('Failed to evaluate all samples & retrying failed samples')
#     count = 0
#     max_guesses_allowed = 10
#     while np.isnan(data.values).any() and count < max_guesses_allowed:
#         run_failed_samples()
#         count += 1
#     data.to_excel(f"{results_path}/{date_formatted}_N{N_samples}_results_retried.xlsx",)
    # calculated_rho_values, calculated_p_values = calculate_spearmanrho()
    # calculated_rho_values.to_excel(f"{results_path}/{date_formatted}_N{N_samples}_rho.xlsx",)

# try:
#     model.evaluate(notify=2)
#     model.table.to_excel(f"{date_formatted}_N{N_samples}_results_{JOB_ID}.xlsx",)
# except:
#     model.table.to_excel(f"{date_formatted}_N{N_samples}_results_failed_{JOB_ID}.xlsx")

# for i in tqdm(samples[RUNS_PER_JOB * JOB_ID : RUNS_PER_JOB * (JOB_ID + 1)]):
#     s = np.array(([i]))
#     # s = np.array(([samples[JOB_NUMBER]))
#     model.load_samples(s)
#     try:
#         model.evaluate(notify = 2)
#         for j in model.table.values:
#             results.loc[len(results)] = j
#         # results.to_excel(f"01182024_results_{JOB_ID}.xlsx")
#         # print('Simulated',len(results))
#     except:
#         failed_samples.loc[len(failed_samples)] = i
#         # failed_samples.to_excel(f"01182024_failed_{JOB_ID}.xlsx")
#         # print('Failed',len(failed_samples))
#
# results.to_excel(f"01192024_failed_results_{JOB_ID}.xlsx")
# failed_samples.to_excel(f"01192024_failedagain{JOB_ID}.xlsx")

# #Load samples from an excel file
# data = pd.read_excel('05282024_samples.xlsx',usecols="B:U")
# samples = data.values
# #if want to run individual parameters to try:
# # model.load_samples(np.array([all[2]]))
# # model.evaluate()
#
# #Create dataframe to store results and failed samples
# cols = model.parameters + model.metrics
# results = pd.DataFrame(columns=cols)
# failed_samples = pd.DataFrame(columns=model.parameters)





# #Use parser to get the job number
# parser = argparse.ArgumentParser()
# parser.add_argument('job_number', type=int)
# JOB_ID = parser.parse_args().job_number
#
# RUNS_PER_JOB = 50

# samples_this_job = samples[RUNS_PER_JOB * JOB_ID : RUNS_PER_JOB * (JOB_ID + 1)]
# model.load_samples(samples_this_job)
#
# try:
#     model.evaluate()
#     model.table.to_excel(f"05282024_results_{JOB_ID}.xlsx")
# except:
#     model.table.to_excel(f"05282024_results_{JOB_ID}_with_error.xlsx")

# for i in tqdm(samples[RUNS_PER_JOB * JOB_ID : RUNS_PER_JOB * (JOB_ID + 1)]):
#     s = np.array(([i]))
#     model.load_samples(s)
#     try:
#         model.evaluate()
#         for j in model.table.values:
#             results.loc[len(results)] = j
#         print('Simulated',len(results))
#     except:
#         failed_samples.loc[len(failed_samples)] = i
#         print('Failed',len(failed_samples))
#
# results.to_excel(f"05282024_results_{JOB_ID}.xlsx")
# failed_samples.to_excel(f"05282024_failed_{JOB_ID}.xlsx")

def calcualte_quantiles_percentage(df):
    qm = [0.05, 0.5, 0.95]
    q_values = [df.quantile(i) for i in qm]
    return f"{q_values[1]*100:.2f}% [{q_values[0]*100:.2f} – {q_values[2]*100:.2f}%]"

mlp_iec = data['Mlp IEC [MM$]']*(1+0.04+0.045+0.09)
tdc = data['Total direct cost [MM$]']
mlpiec_percent = calcualte_quantiles_percentage(mlp_iec/tdc)

fp_mc = data['Feedstock processing MC [$/hr]']*get_annual_factor()/1000000
amc = data['Annual material cost [MM$/yr]']
fpmc_percent = calcualte_quantiles_percentage(fp_mc/amc)

brf_heat = data['Brf heating demand [GJ/hr]']
brf_cool = data['Brf cooling demand [GJ/hr]']
total_heat = data['Total heating demand [GJ/hr]']
total_cool = data['Total cooling demand [GJ/hr]']
brfh_percent = calcualte_quantiles_percentage(brf_heat/total_heat)
brfc_percent = calcualte_quantiles_percentage(brf_cool/total_cool)

fac_elec = data['Facilities electricity consumption [MW]']
mlp_elec = data['Mlp electricity consumption [MW]']
total_elec = data['Total electricity usage [MW]']
fac_percent = calcualte_quantiles_percentage(fac_elec/total_elec)
mlp_percent = calcualte_quantiles_percentage(mlp_elec/total_elec)
fac_cool = data['Facilities cooling demand [GJ/hr]']
mlp_cool = data['Mlp cooling demand [GJ/hr]']
fac_cool_percent = calcualte_quantiles_percentage(fac_cool/total_cool)
mlp_cool_percent = calcualte_quantiles_percentage(mlp_cool/total_cool)

brfh_gwp = data['GWP from bioproducts recovery heating demand [million kg CO2-eq/yr]']
resth_gwp = data['GWP from remaining heating demand [million kg CO2-eq/yr]']
total_gwp = data['Annual system GWP [million kg CO2-eq/yr]']
brfh_gwp_percent = calcualte_quantiles_percentage(brfh_gwp/total_gwp)
resth_gwp_percent = calcualte_quantiles_percentage(resth_gwp/total_gwp)
elec_gwp = data['GWP from cooling demand [million kg CO2-eq/yr]'] + data['GWP from non-cooling demand [million kg CO2-eq/yr]']
elec_gwp_percent = calcualte_quantiles_percentage(elec_gwp/total_gwp)