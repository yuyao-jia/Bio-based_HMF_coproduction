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
parser = argparse.ArgumentParser()
parser.add_argument('job_number', type=int)
# JOB_ID = parser.parse_args().job_number
N_samples = parser.parse_args().job_number # 2000 samples are enough

today = date.today()
date_formatted = today.strftime("%m%d%Y")
# date_formatted = '06132024'
# results_path = "/home/rodrigo/olive_results"
results_path='parameter_tests/test3/test3_repeat'

productivity_base = F.unit.R401.productivity
# lb_pro = 0.17 #from hydrolysate study
lb_pro = productivity_base * 0.75
ub_pro = productivity_base * 1.25
@model.parameter(
    name='Microbial lipid productivity',
    element=F.R401,
    units='g/L/hr',
    kind='coupled',
    distribution=shape.Uniform(lb_pro, ub_pro),
    baseline=productivity_base,
)
def set_lipid_production_productivity(prod):
    F.unit.R401.productivity = prod

# data = pd.read_excel('05282024_baseline_parameters.xlsx',usecols="B:U")
# samples = data.values
# model.load_samples(np.array(([samples[1]])))
# model.evaluate()

# #Generate samples & save them to an excel file
# N_samples = parser.parse_args().job_number # 1000 samples are enough?
model.get_distribution_summary()['Uniform'].to_excel(f'{results_path}/{date_formatted}_N{N_samples}_parameter_distribution.xlsx')
rule = 'L'  # For Latin-Hypercube sampling
np.random.seed(1234)  # For consistent results
samples = model.sample(N_samples, rule)
pd.DataFrame(samples).to_excel(f'{results_path}/{date_formatted}_N{N_samples}_samples.xlsx')

# baseline_metrics.to_excel(f"06122024_baseline_metrics_buyelec_noburnbagasse_N500_{JOB_ID}.xlsx")
# data = pd.read_excel(f'06132024_N500_results_failed_1.xlsx',usecols='B:X',skiprows=[0,2])
# samples = data.values [59:]
# model.load_samples(samples)

# RUNS_PER_JOB = 1000
#
# samples_this_job = samples[RUNS_PER_JOB * JOB_ID : RUNS_PER_JOB * (JOB_ID + 1)]

model.load_samples(samples)
model.evaluate(notify=50)
df_rho, df_p = model.spearman_r()
results = model.table
results.to_excel(f"{results_path}/{date_formatted}_N{N_samples}_results.xlsx",)
df_rho.to_excel(f"{results_path}/{date_formatted}_N{N_samples}_rho.xlsx",)
if np.isnan(results.values).any():
    print('Failed to evaluate all samples & retrying failed samples')
    # try:
    #     data = pd.read_excel(f"{results_path}/{date_formatted}_N{N_samples}_results.xlsx", skiprows=[0, 2], index_col=[0])
    #     total_parameters = len(model.parameters)
    #     cols = data.columns
    #     parameters = cols[0:total_parameters]  # starting index is included, end index is not included
    #     metric = cols[total_parameters:]
    #
    # # path = '/Users/yuyaojia/Documents/Fermentation/Yuyao_TEA_LCA/Bioindustrial-Park/biorefineries/HMF/parameter_tests/test1_500'
    # # path = '/tmp/pycharm_project_525/parameter_tests/test1_500'
    #
    #     failed = pd.DataFrame(columns=cols)
    #     for i in data.index:
    #         if np.isnan(data.loc[i]).any():
    #             failed.loc[i] = data.loc[i]
    # # failed.to_excel('failed_samples.xlsx')
    #
    #     failed_samples = failed[parameters].values
    #     model.load_samples(failed_samples)
    #     model.evaluate(notify=1)
    #
    #     results_failed = model.table.copy()
    #     results_failed.index = failed.index
    #
    #     for i in results_failed.index:
    #         if np.isnan(results_failed.loc[i]).any() == False:
    #             failed = failed.drop(i)
    #             data.loc[i] = results_failed.loc[i].to_list()
    #         else:
    #             print(f"Failed to evaluate {i}")
    #     data.to_excel(f"{results_path}/{date_formatted}_N{N_samples}_results_withfailed.xlsx",)
    #
    #     rho = []
    #     p_value = []
    #     rho_values = pd.DataFrame(index=parameters, columns=metric)
    #     p_values = pd.DataFrame(index=parameters, columns=metric)
    #     for i in metric:
    #         met = data[i].tolist()
    #         for j in parameters:
    #             par = data[j].tolist()
    #             res = stats.spearmanr(par, met)
    #             rho.append(res.correlation)
    #             p_value.append(res.pvalue)
    #         rho_values[i] = rho
    #         p_values[i] = p_value
    #         rho = []
    #         p_value = []
    #     rho_values.to_excel(f"{results_path}/{date_formatted}_N{N_samples}_rho.xlsx",)
    #
    # except:
    #     print('Failed to evaluate failed samples')


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
