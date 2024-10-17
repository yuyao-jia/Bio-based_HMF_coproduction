import biosteam as bst
import numpy as np
import chaospy
from chaospy import distributions as shape
import matplotlib.pyplot as plt
import numpy as np
from Unit_functions.Evaluation_functions.cellulosic_ethanol_tea_yj import create_cellulosic_ethanol_tea as create_tea
from biosteam.evaluation import Model, Metric
from biorefineries.tea.cellulosic_ethanol_tea import CellulosicEthanolTEA, create_cellulosic_ethanol_tea
from biosteam import preferences
from biosteam import report
from biosteam.plots.utils import CABBI_green_colormap
from thermosteam.utils import GG_colors,GG_light_colors
import contourplots
from contourplots import stacked_bar_plot, box_and_whiskers_plot, animated_contourplot
from model_uncertainty_analysis import *
import matplotlib.pyplot as plt
import pandas as pd
import scipy
from scipy import stats
import argparse
from tqdm import tqdm

# steps = 200
# xdatab = np.linspace(0.05,0.15,steps)#oilcane lipid content
# w_data1b = []
# w_data2b = []
# w_data3b = []
# w_data4b = []
# w_data5b = []
# w_data6b = []
# w_data7b = []
# w_data8b = []
#
# def MPSP_at_x1(x):
#     # y = lipid content in dry oilcane
#     feedstock_copy = bst.Stream('feedstock_copy')
#     feedstock_copy.copy_like(feedstock)
#     a = set_oil_fraction(lipid_fraction=x, stream=feedstock_copy,
#                      FFA_fraction=0.1,
#                      z_mass_carbs_baseline=0.149,
#                      PL_fraction=0.1)
#     try:
#         baseline_sys.simulate()
#     except:
#         F.D301.P = 3900
#         baseline_sys.simulate()
#     hmf_mpsp = oilcane_tea.solve_price(purified_hmf)
#     mfpp = oilcane_tea.solve_price(feedstock)
#     biodiesel_mpsp = oilcane_tea.solve_price(biodiesel)
#     gwpenergy = get_hmf_GWP_energy()
#     gwpeconomic = get_hmf_GWP_economic()
#     gwpdis = get_hmf_GWP_all_displaced()
#     gwpbioenergy = get_biodiesel_GWP_energy()
#     gwpbioeconomic = get_biodiesel_GWP_economic()
#     return [hmf_mpsp,mfpp,biodiesel_mpsp,
#             gwpenergy,gwpeconomic,gwpdis,
#             gwpbioenergy, gwpbioeconomic]
# xb=[]
# for i in xdatab:
#         xb.append(i)
# x_oilcane= pd.DataFrame(columns=['Oilcane lipid content [%]',])
# x_oilcane['Oilcane lipid content [%]'] = xb
# x_oilcane.to_excel('x_oilcane_lipid_content.xlsx')
#
# for i in xdatab:
#     a = MPSP_at_x1(i)
#     w_data1b.append(a[0])
#     w_data2b.append(a[1])
#     w_data3b.append(a[2])
#     w_data4b.append(a[3])
#     w_data5b.append(a[4])
#     w_data6b.append(a[5])
#     w_data7b.append(a[6])
#     w_data8b.append(a[7])
#     print(i)
# xydata = x_oilcane
# xydata.insert(1,'HMF MPSP [$/kg]',w_data1b)
# xydata.insert(2,'MFPP [$/kg]',w_data2b)
# xydata.insert(3,'Biodiesel MPSP [$/kg]',w_data3b)
# xydata.insert(4,'HMF GWP (Energy) [kg CO2-eq/kg]',w_data4b)
# xydata.insert(5,'HMF GWP (Economic) [kg CO2-eq/kg]',w_data5b)
# xydata.insert(6,'HMF GWP (Displacement) [kg CO2-eq/MJ]',w_data6b)
# xydata.insert(7,'Biodiesel GWP (Energy) [kg CO2-eq/kg]',w_data7b)
# xydata.insert(8,'Biodiesel GWP (Economic) [kg CO2-eq/kg]',w_data8b)
#
# a = xydata['MFPP [$/kg]']
# xydata['MFPP [$/kg]'] = a*1000
# xydata.columns=['Oilcane lipid content [%]', 'HMF MPSP [$/kg]','MFPP [$/kg]',
#                 'Biodiesel MPSP [$/MT]','HMF GWP (Energy) [kg CO2-eq/kg HMF]','HMF GWP (Economic) [kg CO2-eq/kg HMF]',
#                 'HMF GWP (Displacement) [kg CO2-eq/HMF]', 'Biodiesel GWP (Energy) [kg CO2-eq/kg biodiesel]',
#                 'Biodiesel GWP (Economic) [kg CO2-eq/kg biodiesel]']
# xydata.to_excel('xy_oilcane.xlsx')

# steps = 400
# # cap_b = 3.33e+05
# lb_cap = 1e+9/get_annual_factor() #1 million MT/year
# ub_cap = 5e+9/get_annual_factor() #5 million MT/year
# xdatab = np.linspace(lb_cap,ub_cap,steps)
# w_data1b = []
# w_data2b = []
# w_data3b = []
# w_data4b = []
# w_data5b = []
# w_data6b = []
# w_data7b = []
# w_data8b = []
#
# def MPSP_at_x2(x):
#     # x = feedstock capa
#     feedstock.F_mass = x
#     try:
#         baseline_sys.simulate()
#     except:
#         F.D301.P = 3900
#         baseline_sys.simulate()
#     hmf_mpsp = oilcane_tea.solve_price(purified_hmf)
#     mfpp = oilcane_tea.solve_price(feedstock)*1000
#     biodiesel_mpsp = oilcane_tea.solve_price(biodiesel)
#     gwpenergy = get_hmf_GWP_energy()
#     gwpeconomic = get_hmf_GWP_economic()
#     gwpdis = get_hmf_GWP_all_displaced()
#     gwpbioenergy = get_biodiesel_GWP_energy()
#     gwpbioeconomic = get_biodiesel_GWP_economic()
#     return [hmf_mpsp,mfpp,biodiesel_mpsp,
#             gwpenergy,gwpeconomic,gwpdis,
#             gwpbioenergy, gwpbioeconomic]
# # xb=xdatab
# # x_capacity= pd.DataFrame(columns=['Feedstock flow (kg/hr)','Feedstock capacity (million MT/yr)',])
# # x_capacity['Feedstock flow (kg/hr)'] = xb
# # x_capacity['Feedstock capacity (million MT/yr)'] = xb*get_annual_factor()/1000
# # x_capacity.to_excel('x_capacity.xlsx')
# data = pd.read_excel('x_capacity.xlsx',usecols='B')
# xdatab = data.values
#
# parser = argparse.ArgumentParser()
# parser.add_argument('job_number', type=int)
# JOB_ID = parser.parse_args().job_number
#
# RUNS_PER_JOB =20
#
# for i in tqdm(xdatab[RUNS_PER_JOB * JOB_ID : RUNS_PER_JOB * (JOB_ID + 1)]):
#     a = MPSP_at_x2(i)
#     w_data1b.append(a[0])
#     w_data2b.append(a[1])
#     w_data3b.append(a[2])
#     w_data4b.append(a[3])
#     w_data5b.append(a[4])
#     w_data6b.append(a[5])
#     w_data7b.append(a[6])
#     w_data8b.append(a[7])
#     # print(i)
# xydata = pd.DataFrame(xdatab[RUNS_PER_JOB * JOB_ID : RUNS_PER_JOB * (JOB_ID + 1)])
# xydata.columns = ['Feedstock flow (kg/hr)']
# # xydata['Feedstock flow (kg/hr)'] = xdatab[RUNS_PER_JOB * JOB_ID : RUNS_PER_JOB * (JOB_ID + 1)]
# # xydata['Feedstock capacity (million MT/yr)'] = xdatab[RUNS_PER_JOB * JOB_ID : RUNS_PER_JOB * (JOB_ID + 1)]*get_annual_factor()/1000
# xydata.insert(1,'Feedstock capacity (million MT/yr)',value= xdatab[RUNS_PER_JOB * JOB_ID : RUNS_PER_JOB * (JOB_ID + 1)]*get_annual_factor()/1000/1000000)
# xydata.insert(2,'HMF MPSP [$/kg]',w_data1b)
# xydata.insert(3,'MFPP [$/MT]',w_data2b)
# xydata.insert(4,'Biodiesel MPSP [$/kg]',w_data3b)
# xydata.insert(5,'HMF GWP (Energy) [kg CO2-eq/kg HMF]',w_data4b)
# xydata.insert(6,'HMF GWP (Economic) [kg CO2-eq/kg HMF]',w_data5b)
# xydata.insert(7,'HMF GWP (Displacement) [kg CO2-eq/kg HMF]',w_data6b)
# xydata.insert(8,'Biodiesel GWP (Energy) [kg CO2-eq/kg biodiesel]',w_data7b)
# xydata.insert(9,'Biodiesel GWP (Economic) [kg CO2-eq/kg biodiesel]',w_data8b)
#
# # xydata.columns=['Feedstock capacity', 'HMF MPSP [$/kg]','MFPP [$/kg]',
# #                 'Biodiesel MPSP [$/MT]','HMF GWP (Energy) [kg CO2-eq/kg HMF]','HMF GWP (Economic) [kg CO2-eq/kg HMF]',
# #                 'HMF GWP (Displacement) [kg CO2-eq/HMF]', 'Biodiesel GWP (Energy) [kg CO2-eq/kg biodiesel]',
# #                 'Biodiesel GWP (Economic) [kg CO2-eq/kg biodiesel]']
# xydata.to_excel(f"01202024_xy_capacity_{JOB_ID}.xlsx")

# steps = 400
# hmf_base = purified_hmf.price
# xdatab = np.linspace(0.5,hmf_base,steps)
# w_data1b = []
# w_data2b = []
# # w_data3b = []
# # w_data4b = []
# # w_data5b = []
# # w_data6b = []
# # w_data7b = []
# # w_data8b = []
#
# def MPSP_at_x3(x):
#     purified_hmf.price = x
#     # try:
#     #     baseline_sys.simulate()
#     # except:
#     #     F.D301.P = 3900
#     #     baseline_sys.simulate()
#     # hmf_mpsp = oilcane_tea.solve_price(purified_hmf)
#     mfpp = oilcane_tea.solve_price(feedstock)
#     biodiesel_mpsp = oilcane_tea.solve_price(biodiesel)
#     # gwpenergy = get_hmf_GWP_energy()
#     # gwpeconomic = get_hmf_GWP_economic()
#     # gwpdis = get_hmf_GWP_all_displaced()
#     # gwpbioenergy = get_biodiesel_GWP_energy()
#     # gwpbioeconomic = get_biodiesel_GWP_economic()
#     return [mfpp,biodiesel_mpsp,]
#             # gwpenergy,gwpeconomic,gwpdis,
#             # gwpbioenergy, gwpbioeconomic]
# xb=[]
# for i in xdatab:
#         xb.append(i)
# x_hmf_price= pd.DataFrame(columns=['HMF price [$/kg]',])
# x_hmf_price['HMF price [$/kg]'] = xb
# x_hmf_price.to_excel('x_hmf_price.xlsx')
#
# for i in xdatab:
#     a = MPSP_at_x3(i)
#     w_data1b.append(a[0])
#     w_data2b.append(a[1])
#     # w_data3b.append(a[2])
#     # w_data4b.append(a[3])
#     # w_data5b.append(a[4])
#     # w_data6b.append(a[5])
#     # w_data7b.append(a[6])
#     # w_data8b.append(a[7])
#     print(i)
# xydata = x_hmf_price
# xydata.insert(1,'MFPP [$/MT]',w_data1b)
# xydata.insert(3,'Biodiesel MPSP [$/kg]',w_data2b)
#
# a = xydata['MFPP [$/MT]']
# xydata['MFPP [$/MT]'] = a*1000
# xydata.columns=['HMF price [$/kg]', 'MFPP [$/MT]',
#                 'Biodiesel MPSP [$/kg]']
# xydata.to_excel('xy_hmf_price.xlsx')


# steps = 400
# b_base = biodiesel.price
# xdatab = np.linspace(b_base*0.5,b_base*1.5,steps)
# w_data1b = []
# w_data2b = []
# # w_data3b = []
# # w_data4b = []
# # w_data5b = []
# # w_data6b = []
# # w_data7b = []
# # w_data8b = []
#
# def MPSP_at_x4(x):
#     biodiesel.price = x
#     # try:
#     #     baseline_sys.simulate()
#     # except:
#     #     F.D301.P = 3900
#     #     baseline_sys.simulate()
#     mfpp = oilcane_tea.solve_price(feedstock)
#     hmf_mpsp = oilcane_tea.solve_price(purified_hmf)
#     # gwpenergy = get_hmf_GWP_energy()
#     # gwpeconomic = get_hmf_GWP_economic()
#     # gwpdis = get_hmf_GWP_all_displaced()
#     # gwpbioenergy = get_biodiesel_GWP_energy()
#     # gwpbioeconomic = get_biodiesel_GWP_economic()
#     return [mfpp,hmf_mpsp,]
#             # gwpenergy,gwpeconomic,gwpdis,
#             # gwpbioenergy, gwpbioeconomic]
# xb=[]
# for i in xdatab:
#         xb.append(i)
# x_biodiesel_price= pd.DataFrame(columns=['Biodiesel price [$/kg]',])
# x_biodiesel_price['Biodiesel price [$/kg]'] = xb
# x_biodiesel_price.to_excel('x_biodiesel_price.xlsx')
#
# for i in xdatab:
#     a = MPSP_at_x4(i)
#     w_data1b.append(a[0])
#     w_data2b.append(a[1])
#     # w_data3b.append(a[2])
#     # w_data4b.append(a[3])
#     # w_data5b.append(a[4])
#     # w_data6b.append(a[5])
#     # w_data7b.append(a[6])
#     # w_data8b.append(a[7])
#     print(i)
# xydata = x_biodiesel_price
# xydata.insert(1,'MFPP [$/MT]',w_data1b)
# xydata.insert(2,'HMF MPSP [$/kg]',w_data2b)
#
# a = xydata['MFPP [$/MT]']
# xydata['MFPP [$/MT]'] = a*1000
# xydata.columns=['Biodisel price [$/kg]', 'MFPP [$/MT]',
#                 'HMF MPSP [$/kg]']
# xydata.to_excel('xy_biodiesel_price.xlsx')

steps =100
xdatab = np.linspace(0.8,0.99,steps)
w_data1b = []
w_data2b = []
w_data3b = []
w_data4b = []
w_data5b = []
w_data6b = []
w_data7b = []
w_data8b = []

def MPSP_at_x5(x):
    F.BT701.boiler_efficiency = x
    try:
        baseline_sys.simulate()
    except:
        F.D301.P = 3900
        baseline_sys.simulate()
    hmf_mpsp = oilcane_tea.solve_price(purified_hmf)
    mfpp = oilcane_tea.solve_price(feedstock)*1000
    biodiesel_mpsp = oilcane_tea.solve_price(biodiesel)
    gwpenergy = get_hmf_GWP_energy()
    gwpeconomic = get_hmf_GWP_economic()
    gwpdis = get_hmf_GWP_all_displaced()
    gwpbioenergy = get_biodiesel_GWP_energy()
    gwpbioeconomic = get_biodiesel_GWP_economic()
    return [hmf_mpsp,mfpp,biodiesel_mpsp,
            gwpenergy,gwpeconomic,gwpdis,
            gwpbioenergy, gwpbioeconomic]

# x_boilereff= pd.DataFrame(columns=['Boiler efficiency (%)',])
# x_boilereff['Boiler efficiency (%)'] = xdatab
# x_boilereff.to_excel('x_boiler_efficiency.xlsx')
data = pd.read_excel('x_boiler_efficiency.xlsx',usecols='B')
xdatab = data.values

parser = argparse.ArgumentParser()
parser.add_argument('job_number', type=int)
JOB_ID = parser.parse_args().job_number

RUNS_PER_JOB =3

for i in tqdm(xdatab[RUNS_PER_JOB * JOB_ID : RUNS_PER_JOB * (JOB_ID + 1)]):
    a = MPSP_at_x5(i)
    w_data1b.append(a[0])
    w_data2b.append(a[1])
    w_data3b.append(a[2])
    w_data4b.append(a[3])
    w_data5b.append(a[4])
    w_data6b.append(a[5])
    w_data7b.append(a[6])
    w_data8b.append(a[7])
    # print(i)
xydata =  pd.DataFrame(xdatab[RUNS_PER_JOB * JOB_ID : RUNS_PER_JOB * (JOB_ID + 1)])
xydata.columns = ['Boilier efficiency [%]']
# xydata['Boilier efficiency (%)'] = xdatab[RUNS_PER_JOB * JOB_ID : RUNS_PER_JOB * (JOB_ID + 1)]
xydata.insert(1,'HMF MPSP [$/kg]',w_data1b)
xydata.insert(2,'MFPP [$/MT]',w_data2b)
xydata.insert(3,'Biodiesel MPSP [$/kg]',w_data3b)
xydata.insert(4,'HMF GWP (Energy) [kg CO2-eq/kg HMF]',w_data4b)
xydata.insert(5,'HMF GWP (Economic) [kg CO2-eq/kg HMF]',w_data5b)
xydata.insert(6,'HMF GWP (Displacement) [kg CO2-eq/kg HMF]',w_data6b)
xydata.insert(7,'Biodiesel GWP (Energy) [kg CO2-eq/kg biodiesel]',w_data7b)
xydata.insert(8,'Biodiesel GWP (Economic) [kg CO2-eq/kg biodiesel]',w_data8b)

# xydata.columns=['Feedstock capacity', 'HMF MPSP [$/kg]','MFPP [$/kg]',
#                 'Biodiesel MPSP [$/MT]','HMF GWP (Energy) [kg CO2-eq/kg HMF]','HMF GWP (Economic) [kg CO2-eq/kg HMF]',
#                 'HMF GWP (Displacement) [kg CO2-eq/HMF]', 'Biodiesel GWP (Energy) [kg CO2-eq/kg biodiesel]',
#                 'Biodiesel GWP (Economic) [kg CO2-eq/kg biodiesel]']
xydata.to_excel(f"01212024_xy_boilereff_{JOB_ID}.xlsx")
