# -*- coding: utf-8 -*-
"""
Created on Wed Sep  6 14:19:25 2023

@author: lavan
"""

import biosteam as bst
import numpy as np
# import chaospy
# from chaospy import distributions as shape
import matplotlib.pyplot as plt
import numpy as np
# from Unit_functions.Evaluation_functions.cellulosic_ethanol_tea_yj import create_cellulosic_ethanol_tea as create_tea
# from biosteam.evaluation import Model, Metric
# from biorefineries.tea.cellulosic_ethanol_tea import CellulosicEthanolTEA, create_cellulosic_ethanol_tea
# from biosteam import preferences
# from biosteam import report
# from biosteam.plots.utils import CABBI_green_colormap
# from thermosteam.utils import GG_colors,GG_light_colors
# import contourplots
# from contourplots import stacked_bar_plot, box_and_whiskers_plot, animated_contourplot
from model_uncertainty_analysis import *
import matplotlib.pyplot as plt
import pandas as pd
import argparse
from tqdm import tqdm
from datetime import date

parser = argparse.ArgumentParser()
parser.add_argument('step_number', type=int)
steps = parser.parse_args().step_number

today = date.today()
date_formatted = today.strftime("%m%d%Y")
results_path='contourplots/vegetative_lipid/test2'

# feedstock lipid content vs. vegetative lipid recovery after pretreatment & saccharification
# steps = 30
x_data = np.linspace(0.05,0.15,steps)#lipid content in oilcane
y_data = np.linspace(0.5,0.95,steps)#vegetative lipid recovery after pretreatment & saccharification

x_plot=[]
y_plot=[]
for j in y_data:
    for i in x_data:
        x_plot.append(i)
        y_plot.append(j)

def get_main_indicators():
    #main economic
    mfpp = get_MFPP()
    hmf_mpsp = get_MPSP_hmf()
    biodiesel_mpsp = get_MPSP_biodiesel()

    #main environmental
    gwpenergy = get_hmf_GWP_energy()
    gwpenergy_grave = get_hmf_GWP_energy_grave()

    #other environmental (cradle to gate)
    gwpeconomic = get_hmf_GWP_economic()
    gwpdis = get_hmf_GWP_displacement()
    gwpbioenergy = get_biodiesel_GWP_energy()
    gwpbioeconomic = get_biodiesel_GWP_economic()

    #other environmental (cradle to grave)
    gwpeco_grave = get_hmf_GWP_economic_grave()
    gwpdis_grave = get_hmf_GWP_displacement_grave()
    gwpbioenergy_grave = get_biodiesel_GWP_energy_grave()
    gwpbioeconomic_grave = get_biodiesel_GWP_economic_grave()

    return ([mfpp,hmf_mpsp, biodiesel_mpsp,
            gwpenergy,gwpenergy_grave,
            gwpeconomic,gwpdis,gwpbioenergy,gwpbioeconomic,
            gwpeco_grave,gwpdis_grave,gwpbioenergy_grave,gwpbioeconomic_grave])
def get_production_and_utilities():
    furfural_yield = get_total_yield_furfural()
    hmf_yield = get_total_yield_hmf()
    biodiesel_yield = get_total_yield_biodiesel()
    crude_glycerol_yield = get_total_yield_crude_glycerol()
    electricity_yield = get_net_electricity_yield()

    heating = get_system_heating_demand()
    cooling = get_system_cooling_demand()
    electricity = get_total_electricity_demand()
    purchased_elec = get_purchased_electricity()
    ng_usage = get_natural_gas_demand()
    return ([furfural_yield,hmf_yield,biodiesel_yield,crude_glycerol_yield,electricity_yield,
            heating,cooling,electricity,purchased_elec,ng_usage])
def get_technical_metrics():
    furfural_purity = get_purity_furfural()
    hmf_purity = get_purity_hmf()

    veglipid_rec = get_vegetative_lipid_recovery()
    veglipid_extractioneff = get_vegetative_lipid_extraction_efficiency()
    veglipid_extraction_verifier = verify_vegetative_lipid_extraction_efficiency()

    miclipid_rec = get_microbial_lipid_recovery()
    miclipid_extractioneff = get_microbial_lipid_extraction_efficiency()
    miclipid_extraction_verifier = verify_microbial_lipid_extraction_efficiency()
    micrlipid_titer = get_actual_titer()
    micrlipid_productivity = get_actual_productivity()
    return ([furfural_purity,hmf_purity,
             veglipid_rec,veglipid_extractioneff,veglipid_extraction_verifier,
             miclipid_rec,miclipid_extractioneff,miclipid_extraction_verifier,
             micrlipid_titer,micrlipid_productivity])

def get_economic_metrics():
    tci = get_overall_TCI()
    tdc = get_overall_TDC()
    aoc = get_operating_cost()
    amc = get_material_cost()
    npv = get_NPV()
    irr = solve_IRR()
    return ([tci,tdc,aoc,amc,npv,irr])

def get_gwp_total_with_breakdown():
    annual_sysgwp = get_net_GWP()
    feedstock_gwp = get_feedstock_GWP()
    othermaterial_gwp = get_other_material_GWP()
    brf_heating_gwp = get_brf_heating_GWP()
    remaining_heating_gwp = get_remaining_heating_GWP()
    cooling_gwp = get_cooling_GWP()
    noncooling_gwp = get_non_cooling_GWP()
    gwp_verifier = verify_GWP()

    gwp_gate = [annual_sysgwp,feedstock_gwp,othermaterial_gwp,
                brf_heating_gwp,remaining_heating_gwp,cooling_gwp,
                noncooling_gwp,gwp_verifier]

    annual_sysgwp_grave = get_net_GWP_grave()
    feedstock_gwp_grave = get_feedstock_GWP_grave()
    othermaterial_gwp_grave = get_other_material_GWP_grave()
    brf_heating_gwp_grave = get_brf_heating_GWP_grave()
    remaining_heating_gwp_grave = get_remaining_heating_GWP_grave()
    cooling_gwp_grave = get_cooling_GWP_grave()
    noncooling_gwp_grave = get_non_cooling_GWP_grave()
    gwp_verifier_grave = verify_GWP_grave()
    gwp_grave = [annual_sysgwp_grave,feedstock_gwp_grave,othermaterial_gwp_grave,
                brf_heating_gwp_grave,remaining_heating_gwp_grave,cooling_gwp_grave,
                 noncooling_gwp_grave, gwp_verifier_grave]
    return gwp_gate,gwp_grave

def get_section_metrics():
    fp_iec = get_fp_cost()
    fp_mc = get_fp_material_cost()
    fp_heating = get_fp_heating_demand()
    fp_cooling = get_fp_cooling_demand()
    fp_electricity = get_fp_electricity_consumption()
    fp = [fp_iec,fp_mc,fp_heating,fp_cooling,fp_electricity]

    pt_iec = get_pt_cost()
    pt_mc = get_pt_material_cost()
    pt_heating = get_pt_heating_demand()
    pt_cooling = get_pt_cooling_demand()
    pt_electricity = get_pt_electricity_consumption()
    pt = [pt_iec,pt_mc,pt_heating,pt_cooling,pt_electricity]

    brf_iec = get_brf_cost()
    brf_mc = get_brf_material_cost()
    brf_heating = get_brf_heating_demand()
    brf_cooling = get_brf_cooling_demand()
    brf_electricity = get_brf_electricity_consumption()
    brf = [brf_iec,brf_mc,brf_heating,brf_cooling,brf_electricity]

    mlp_iec = get_mlp_cost()
    mlp_mc = get_mlp_material_cost()
    mlp_heating = get_mlp_heating_demand()
    mlp_cooling = get_mlp_cooling_demand()
    mlp_electricity = get_mlp_electricity_consumption()
    mlp = [mlp_iec,mlp_mc,mlp_heating,mlp_cooling,mlp_electricity]

    bp_iec = get_bp_cost()
    bp_mc = get_bp_material_cost()
    bp_heating = get_bp_heating_demand()
    bp_cooling = get_bp_cooling_demand()
    bp_electricity = get_bp_electricity_consumption()
    bp = [bp_iec,bp_mc,bp_heating,bp_cooling,bp_electricity]

    wwt_iec = get_wwt_cost()
    wwt_mc = get_wwt_material_cost()
    wwt_heating = get_wwt_heating_demand()
    wwt_cooling = get_wwt_cooling_demand()
    wwt_electricity = get_wwt_electricity_consumption()
    wwt = [wwt_iec,wwt_mc,wwt_heating,wwt_cooling,wwt_electricity]

    chp_iec = get_chp_cost()
    chp_mc = get_chp_material_cost()
    chp_heating = get_chp_heating_demand()
    chp_cooling = get_chp_cooling_demand()
    chp_electricity = get_chp_electricity_consumption()
    chp = [chp_iec,chp_mc,chp_heating,chp_cooling,chp_electricity]

    facilities_iec = get_facilities_cost()
    facilities_mc = get_facilities_material_cost()
    facilities_heating = get_facilities_heating_demand()
    facilities_cooling = get_facilities_cooling_demand()
    facilities_electricity = get_facilities_electricity_consumption()
    facilities = [facilities_iec,facilities_mc,facilities_heating,facilities_cooling,facilities_electricity]
    return fp,pt,brf,mlp,bp,wwt,chp,facilities

def metrics_at_x_and_y_lipid(x,y):
    #x = lipid content in dry oilcane
    feedstock_copy = bst.Stream('feedstock_copy')
    feedstock_copy.copy_like(feedstock)
    a = set_oil_fraction(lipid_fraction=x, stream=feedstock_copy,
                         FFA_fraction=0.1,
                         z_mass_carbs_baseline=0.149,
                         PL_fraction=0.1)
    feedstock.copy_like(a)
    #y = vegetative lipid recovery after pretreatment & saccharification
    F.U402.isplit["Lipid"] = y
    baseline_sys.simulate()
    #get main indicators
    main_indicators = get_main_indicators()
    #get production and utilities
    production_utilities = get_production_and_utilities()
    #get technical metrics
    technical_metrics = get_technical_metrics()
    #get economic metrics
    economic_metrics = get_economic_metrics()
    #get gwp total with breakdown
    gwp_total_gate,gwp_total_grave = get_gwp_total_with_breakdown()
    #get section metrics
    fp,pt,brf,mlp,bp,wwt,chp,facilities = get_section_metrics()
    section_metrics = fp+pt+brf+mlp+bp+wwt+chp+facilities
    return (main_indicators + production_utilities +
            technical_metrics + economic_metrics +
            gwp_total_gate + gwp_total_grave +
            section_metrics)

metric_names = [i[1] for i in baseline_metrics.index]
results = pd.DataFrame(columns = ['Feedstock lipid content [%]','Vegetative lipid recovery after pretreatment & saccharification [%]'] + metric_names)

for i in tqdm(np.arange(len(x_plot))):
    xi = x_plot[i]
    yi = y_plot[i]
    ri = metrics_at_x_and_y_lipid(xi, yi)
    di = [xi, yi] + ri
    results.loc[i] = di

results.to_excel(f"{results_path}/{date_formatted}_veglipid_{steps}_results.xlsx")
baseline_metrics.to_excel(f"{results_path}/{date_formatted}_baseline_metrics_for_verification.xlsx")
##generate color bar
# from biosteam.utils import colors
# import matplotlib.colors as mcolors
# cstep=6
# def mkstr(steps):
#     cgen = CABBI_green_colormap(N_levels=steps)
#     ccs = [[int(_c * 255) for _c in cgen(c)] for c in range(cgen.N)]
#     cmap_hex = [mcolors.rgb2hex(cgen(_c)) for _c in range(cgen.N)]
#     out_str = f"AddCustom('colormap', 'CABBI-step{steps}', ((-1,0,0,0)"
#     for _color in ccs:
#         out_str += f", ({', '.join([str(c) for c in _color])})"
#     out_str += "))"
#     print(f"Add this to the file:\n{out_str}")
#     return out_str
# CABBI_colors = (colors.CABBI_yellow.RGBn,
#                 colors.CABBI_green.RGBn,
#                 colors.CABBI_teal_green.shade(75).RGBn)

# To combine data from different files
# a = pd.read_excel('01182024_results_0.xlsx')
# results = pd.DataFrame(columns=a.columns)
# for i in range(50):
#     ri = pd.read_excel(f"01182024_results_{i}.xlsx")
#     r = pd.concat([results,ri])
#     results = r

#start
# ssh ubuntu
# tums a

#end
# ctrl+b d to detach anywhere
# ctrl+d to exit

#convert economic HMF GWP to float
#a = results[cols[51]] #economic HMF GWP
# b = []
# for i in range(982):
#     if type(a[i]) is not float:
#         bi = float(a[i][1:len(a[i])-1])
#     else:
#         bi = a[i]
#     b.append(bi)
# results[cols[51]]= b

#download files from ubuntu to local in terminal located in local: last folder is where fileds are downloaded to
# scp 'ubuntu:/tmp/pycharm_project_525/01222024*.xlsx' Downloads/

#calculate %changes in metrics per 1% change in x or y
# cmpsp1 = []
# cmfpp1 = []
# for i in range(20):
#     di = data.values[i*20:(i+1)*20]
#     for j in range(19):
#         xi1 = di[j][0] #change to 1 if for y values, same for below
#         xi2 = di[j+1][0]
#         mpspi1 = di[j][2]
#         mpspi2 = di[j+1][2]
#         mfppi1 = di[j][3]
#         mfppi2 = di[j+1][3]
#         delta_mpsp = (mpspi2-mpspi1)/(xi2-xi1)*0.01
#         delta_mfpp = (mfppi2-mfppi1)/(xi2-xi1)*0.01
#         cmpsp1.append(delta_mpsp)
#         cmfpp1.append(delta_mfpp)