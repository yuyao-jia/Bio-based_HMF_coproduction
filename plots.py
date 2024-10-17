# -*- coding: utf-8 -*-
"""
Created on Wed Sep  6 14:19:25 2023

@author: lavan
"""

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
from system_setup.lca_characterization_factors_yj import set_GWPCF, GWP
import scipy
from scipy import stats
# from Uncertainty_analysis_functions.spearmanr import get_spearmanr
#TODO: try sankey plot for NF rejection factors
#TODO: ask about heating & cooling duty
#stacked bar plot
# %% Plot GWP breakdown stacked plot
#TODO: plot system (without product) GWP breakdown stacked plot?

# set_GWPCF(biodiesel, 'biodiesel displacement')
# set_GWPCF(crude_glycerol, 'glycerine displacement', dilution=0.8)
# set_GWPCF(purified_furfural, 'furfural displacement')
# set_GWPCF(purified_hmf,'HMF displacement')
# coproducts1=[F.purified_hmf,F.purified_furfural, F.crude_glycerol] #assume HMF is main product

# cost_breakdown = bst.plots.plot_unit_groups(area_groups, fraction=True)
# cost_breakdown[0].set_figwidth(12)
# cost_breakdown[0].savefig(fname='cost_breakdown_by_area',bbox_inches='tight',dpi=100)

#sum of GWP based on hmf
breakdown=[get_feedstock_GWP(), get_other_materials_impact(),
                             get_hmf_purification_heating_demand_GWP(),get_heating_demand_GWP_non_hmf_purification(),
                             get_electricity_demand_GWP(),
                             get_total_non_BT_direct_emissions_GWP()]
a = sum(breakdown)
# % of GWP of each category based on hmf
b=[]
for i in breakdown:
    b.append(i*100/a)
df_gwp =pd.DataFrame({"GWP breakdown":b})
df_gwp.index = ['Feedstock',
                'Other materials',
                'Heating demand (HMF& Furfural Purification)',
                'Heating demand (non HMF& Furfural Purification)',
                'Electricity demand',
                'Total direct emissions (non BT)',]

# stacked_bar_plot(df_gwp,y_ticks = [-20,0,50,100],
#                  y_label = 'Annual system GWP100',y_units='kg CO2-eq/yr',
#                  fig_width=4,
#                  colors =contourplots.utils.defaults_dict['colors']['Guest_Group_TEA_Breakdown'],
#                  filename = 'Net_GWP_breakdown_v1',)
#
# c = sum([get_feedstock_GWP(),
#          get_other_materials_impact(),
#          get_ng_GWP(),#natrual gas
#          get_total_direct_BT_emissions_GWP(),
#          get_total_non_BT_direct_emissions_GWP(),
#          ])
# #GWP breakdown based on hmf
# d = [
#      get_feedstock_GWP()*100/c,
#      get_other_materials_impact()*100/c,
#      get_ng_GWP()*100/c,
#      get_total_direct_BT_emissions_GWP()*100/c,
#      get_total_non_BT_direct_emissions_GWP()*100/c,
#      ]
# #
# df_gwp2 =pd.DataFrame({"GWP breakdown":d})
# df_gwp2.index = ['Feedstock',
#                 'Other feed inputs',
#                 'Natural gas',
#                 'Total direct BT emissions GWP',
#                 'Total non BT direct emissions GWP',]

# stacked_bar_plot(df_gwp2,y_ticks = [-25,0,25,50,75,100],
#                  y_label='Annual system GWP100',y_units='kg CO2-eq/yr',
#                  colors=contourplots.utils.defaults_dict['colors']['Guest_Group_TEA_Breakdown'],
#                  filename = 'Net_GWP_breakdown_v2',)
#
#HMW MPSP & MFPP & HMF GWP sensitivity plots
# data=pd.read_excel('results_12032023.xlsx',header=1,skiprows=[2])
# col_index = data.columns
# parameters = col_index[1:23]  # starting index is included, end index is not included
# metric = col_index[23:]

# # get_spearmanr(simulations_results = 'results.xlsx')
# rho = pd.read_excel('calculated_rho_12032023.xlsx')
# parameter_names = parameters

# MPSP_HMF_fig = bst.plots.plot_spearman_1d(rhos=rho[metric[0]].tolist(),
#                            index=parameter_names,
#                            name='minimum product selling price [$/kg] of HMF',
#                            color = GG_colors.blue.RGBn)
# MPSP_HMF_fig[0].set_figwidth(10)
# MPSP_HMF_fig[0].set_figheight(10)
# MPSP_HMF_fig[0].savefig(fname='MPSP_HMF_rho',bbox_inches='tight',dpi=100)
#
# MFPP_fig = bst.plots.plot_spearman_1d(rhos=rho[metric[1]].tolist(),
#                            index=parameter_names,
#                            name='maximum feedstock purchase price [$/kg] of oilcane',color = GG_colors.blue.RGBn)
# MFPP_fig[0].set_figwidth(10)
# MFPP_fig[0].set_figheight(10)
# # MFPP_fig[0].font = 'Arial'
# MFPP_fig[0].savefig(fname='MFPP_rho',bbox_inches='tight',dpi=100)
#
# GWP_economic_fig = bst.plots.plot_spearman_1d(rhos=rho[metric[21]].tolist(),
#                            index=parameter_names,
#                            name='HMF global warming potential using economic allocation [kg CO2-eq/kg HMF]',
#                             color = GG_colors.blue.RGBn)
# GWP_economic_fig[0].set_figwidth(14)
# GWP_economic_fig[0].set_figheight(10)
# GWP_economic_fig[0].savefig(fname='GWP_HMF_econnomic_rho',bbox_inches='tight',dpi=100)
#
# GWP_energy_fig = bst.plots.plot_spearman_1d(rhos=rho[metric[22]].tolist(),
#                            index=parameter_names,
#                            name='HMF global warming potential using energy allocation [kg CO2-eq/kg HMF]',
#                             color = GG_colors.blue.RGBn)
# GWP_energy_fig[0].set_figwidth(14)
# GWP_energy_fig[0].set_figheight(10)
# GWP_energy_fig[0].savefig(fname='GWP_HMF_energy_rho',bbox_inches='tight',dpi=100)
#
# GWP_displacement_fig = bst.plots.plot_spearman_1d(rhos=rho[metric[23]].tolist(),
#                            index=parameter_names,
#                            name='HMF global warming potential using displacement allocation [kg CO2-eq/kg HMF]',
#                             color = GG_colors.blue.RGBn)
# GWP_displacement_fig[0].set_figwidth(14)
# GWP_displacement_fig[0].set_figheight(10)
# GWP_displacement_fig[0].savefig(fname='GWP_HMF_displacement_rho',bbox_inches='tight',dpi=100)

# #Perform analysis for contourplots
# #range of x and y data should be extreme values
# x_data = biodiesel_price = np.linspace(price['Biodiesel']*0.75,price['Biodiesel']*1.25,25)
# y_data = furfural_price = np.linspace(price['Purified_furfural']*0.75,price['Purified_furfural']*1.25,25)
# steps = 10
# x_data = np.linspace(0.1,0.5,steps)
# y_data = np.linspace(0.1,0.6,steps)
# w_data1 = []
# w_data2 = []
# z_data=[1,]
# def MPSP_at_x_and_y_1(x,y):
#     F.U301.rejection_factors['HMF'] = x
#     F.U302.rejection_factors['HMF'] = y
#     baseline_sys.simulate()
#     return [oilcane_tea.solve_price(purified_hmf),oilcane_tea.solve_price(feedstock)]
#
# for i in x_data:
#     w_data1.append([])
#     w_data2.append([])
#     for j in y_data:
#         w_data1[-1].append(MPSP_at_x_and_y_1(i,j)[0])
#         w_data2[-1].append(MPSP_at_x_and_y_1(i, j)[1])
#
# contourplots.animated_contourplot(w_data_vs_x_y_at_multiple_z=[w_data1], # shape = z * x * y # values of the metric you want to plot on the color axis; e.g., MPSP
#                                 x_data=x_data, # x axis values
#                                 y_data=y_data, # y axis values
#                                 z_data=z_data, # z axis values
#                                 x_label= "Retention of HMF in 1st nanofitlration",
#                                 y_label= "Retention of HMF in 2nd nanofitlration",
#                                 z_label= "ignore", # title of the z axis
#                                 w_label= "HMF MPSP", # title of the color axis
#                                 x_ticks= np.linspace(0.1,0.5,5),#should be same as range of x_data
#                                 y_ticks = np.linspace(0.1,0.6,6),
#                                 z_ticks=[0,1,2],
#                                 w_levels=np.array([i for i in range(10,45,5)]), #number of different shades of color
#                                   # levels for unlabeled, filled contour areas (labeled and ticked only on color bar)
#                                 w_ticks=np.array([i for i in range(10,45,5)]),#number on the plots on the line that represents MPSP
#                                 x_units="%",
#                                 y_units="%",
#                                 z_units=" ",
#                                 w_units="$/kg",
#                                 fmt_clabel=lambda cvalue: "{:.2f}".format(cvalue), # format of contour labels
#                                 cmap=CABBI_green_colormap(), # can use 'viridis' or other default matplotlib colormaps
#                                 cbar_ticks= np.array([i for i in range(10,45,5)]),#legends color, should be same num of w_ticks
#                                 clabel_fontsize=8,
#                                 z_marker_color='g', # default matplotlib color names
#                                 axis_title_fonts={'size': {'x': 7, 'y':7,
#                                                            'z':7, 'w':7}},
#                                 fps=3, # animation frames (z values traversed) per second
#                                 n_loops='inf', # the number of times the animated contourplot should loop animation over z; infinite by default
#                                 animated_contourplot_filename='a_MPSP_contourplot_retention', # file name to save animated contourplot as (no extensions)
#                                 keep_frames=True, # leaves frame PNG files undeleted after running; False by default
#                                 )


# steps = 20
# x_data2 = np.linspace(0.01,0.04,steps)
# y_data2 = np.linspace(1.0,2.5,steps)
# w_data3 = []
# w_data4 = []
# z_data2=[1,]
# def MPSP_at_x_and_y_2(x,y):
#     feedstock.price = x
#     biodiesel.price = y
#     # baseline_sys.simulate()
#     return [oilcane_tea.solve_price(purified_hmf),oilcane_tea.solve_price(feedstock)]
#     # return oilcane_tea.solve_price(purified_hmf)
#
# for j in y_data2:
#     w_data3.append([])
#     w_data4.append([])
#     for i in x_data2:
#         w_data3[-1].append(MPSP_at_x_and_y_2(i,j)[0])
#         w_data4[-1].append(MPSP_at_x_and_y_2(i,j)[1])

# # %% Plot results
# contourplots.animated_contourplot(w_data_vs_x_y_at_multiple_z=[w_data2], # shape = z * x * y # values of the metric you want to plot on the color axis; e.g., MPSP
#                                 x_data=x_data2, # x axis values
#                                 y_data=y_data2, # y axis values
#                                 z_data=z_data2, # z axis values
#                                 x_label= "Feedstock price",
#                                 y_label= "Biodiesel price",
#                                 z_label= "ignore", # title of the z axis
#                                 w_label= "HMF MPSP", # title of the color axis
#                                 x_ticks= np.linspace(0.01,0.04,5),#should be same as range of x_data
#                                 y_ticks = np.linspace(1,2.5,5),
#                                 z_ticks=[0,1,2],
#                                 w_levels=np.array([i for i in range(0,40,5)]), #number of different shades of color
#                                   # levels for unlabeled, filled contour areas (labeled and ticked only on color bar)
#                                 w_ticks=np.array([i for i in range(0,40,5)]),#number on the plots on the line that represents MPSP
#                                 x_units="$",
#                                 y_units="$",
#                                 z_units=" ",
#                                 w_units="$/kg",
#                                 fmt_clabel=lambda cvalue: "{:.2f}".format(cvalue), # format of contour labels
#                                 cmap=CABBI_green_colormap(), # can use 'viridis' or other default matplotlib colormaps
#                                 cbar_ticks= np.linspace(0,40,num = 5),#legends color, should be same num of w_ticks
#                                 clabel_fontsize=8,
#                                 z_marker_color='g', # default matplotlib color names
#                                 axis_title_fonts={'size': {'x': 7, 'y':7,
#                                                            'z':7, 'w':7}},
#                                 fps=3, # animation frames (z values traversed) per second
#                                 n_loops='inf', # the number of times the animated contourplot should loop animation over z; infinite by default
#                                 animated_contourplot_filename='a_MPSP_contourplot', # file name to save animated contourplot as (no extensions)
#                                 keep_frames=True, # leaves frame PNG files undeleted after running; False by default
#                                 )
# contourplots.animated_contourplot(w_data_vs_x_y_at_multiple_z=[price_data], # shape = z * x * y # values of the metric you want to plot on the color axis; e.g., MPSP
#                                 x_data=np.sort(dol_per_kg_soybean), # x axis values
#                                 y_data=premiums, # y axis values
#                                 z_data=z_data, # z axis values
#                                 x_label= "soybean oilseed price",
#                                 y_label= "Premium on soybean seeds",
#                                 z_label= "ignore", # title of the z axis
#                                 w_label="MPSP", # title of the color axis
#                                 x_ticks= dol_per_kg_soybean,
#                                 y_ticks = premiums,
#                                 z_ticks=[0,1,2],
#                                 w_levels=np.array([i for i in range(2,11)]), # levels for unlabeled, filled contour areas (labeled and ticked only on color bar)
#                                 w_ticks=np.array([i for i in range(3,11,1)]),
#                                 x_units="$/kg",
#                                 y_units="$/kg",
#                                 z_units=" ",
#                                 w_units="$/kg",
#                                 fmt_clabel=lambda cvalue: "{:.2f}".format(cvalue), # format of contour labels
#                                 cmap=CABBI_green_colormap(), # can use 'viridis' or other default matplotlib colormaps
#                                 cbar_ticks= [3,4,5,6,7,8,9],
#                                 clabel_fontsize = 8,
#                                 # np.linspace(3,12,num = 12),
#                                 z_marker_color='g', # default matplotlib color names
#                                 axis_title_fonts={'size': {'x': 7, 'y':7,
#                                                             'z':7, 'w':7}},
#                                 fps=3, # animation frames (z values traversed) per second
#                                 n_loops='inf', # the number of times the animated contourplot should loop animation over z; infinite by default
#                                 animated_contourplot_filename='MPSP_contourplot_', # file name to save animated contourplot as (no extensions)
#                                 keep_frames=True, # leaves frame PNG files undeleted after running; False by default
#                                 )
#star baseline value on contourplots
# #contourplots: membrane lifetime vs membrane purchase cost vs MPSP?
#
#
# box and whisker plots
# #list values of MPSP and GWP from 1000 Monte Carlo simulations
# MPSP_HMF= data[metric[0]].tolist()
# MFPP = data[metric[1]].tolist()
# Net_GWP = (data[metric[21]]).tolist()
# GWP_hmf_economic = data[metric[22]].tolist()
# GWP_hmf_energy = data[metric[23]].tolist()
# GWP_hmf_displacement = data[metric[26]].tolist()
# GWP_hmf_displacement_eco = data[metric[27]].tolist()
# GWP_hmf_displacement_energy = data[metric[28]].tolist()
# bs_metrics = pd.read_excel('baseline_metrics_12032023.xlsx',usecols='C').T
# bs_metrics.columns=metric
#
# box_and_whiskers_plot(MPSP_HMF,baseline_values=[bs_metrics[metric[0]]],y_ticks=[10,20,30,40],y_label = 'Minimum product selling price of HMF',y_units = '$.Kg –1',filename='MPSP_HMF_box_and_whiskers_plot')
# box_and_whiskers_plot(MFPP,baseline_values=[bs_metrics[metric[1]]],y_ticks=[0,0.05,0.10],y_label = 'Maximum feedstock purchasing price',y_units = '$.Kg –1',filename='MFPP_box_and_whiskers_plot')# box_and_whiskers_plot(GWP_hmf_economic,baseline_values=[bs_metrics[metric[21]],y_ticks=[80,100,120,140,160,160],y_label = 'HMF Global warming potential Economic',y_units = 'kg CO2-eq/kg–1',filename='GWP_HMF_economic_box_and_whiskers_plot')
# box_and_whiskers_plot(Net_GWP,baseline_values=[bs_metrics[metric[21]]],y_ticks=[0.5*1000000000,1*1000000000,1.5*1000000000],y_label = 'Annual System Global warming potential',y_units = 'kg CO2-eq/yr',filename='Net_GWP_box_and_whiskers_plot')
# box_and_whiskers_plot(GWP_hmf_economic,baseline_values=[bs_metrics[metric[22]]],y_ticks=[0,25,50,75,100],y_label = 'HMF Global warming potential Economic',y_units = 'kg CO2-eq/kg–1',filename='GWP_HMF_economic_box_and_whiskers_plot')
# box_and_whiskers_plot(GWP_hmf_energy,baseline_values=[bs_metrics[metric[23]]],y_ticks=[0,5,10],y_label = 'HMF Global warming potential Energy',y_units = 'kg CO2-eq/kg–1',filename='GWP_HMF_energy_box_and_whiskers_plot')
# box_and_whiskers_plot(GWP_hmf_displacement,baseline_values=[bs_metrics[metric[26]]],y_ticks=[100,150,200],y_label = 'HMF Global warming potential Displacement',y_units = 'kg CO2-eq/kg–1',filename='GWP_HMF_displacement_box_and_whiskers_plot')
# box_and_whiskers_plot(GWP_hmf_displacement_eco,baseline_values=[bs_metrics[metric[27]]],y_ticks=[0,25,50,75,100],y_label = 'HMF Global warming potential Displacement Economic',y_units = 'kg CO2-eq/kg–1',filename='GWP_HMF_displacement_eco_box_and_whiskers_plot')
# box_and_whiskers_plot(GWP_hmf_displacement_energy,baseline_values=[bs_metrics[metric[28]]],y_ticks=[0,5,10],y_label = 'HMF Global warming potential Displacement Energy',y_units = 'kg CO2-eq/kg–1',filename='GWP_HMF_displacement_energy_box_and_whiskers_plot')

#MPSP at different internal rates of return (IRRs)

#quantile data
# q = pd.DataFrame(columns=['5th','25th','median','75th','95th'])
# col=q.columns
# qdata = [0.05,0.25,0.5,0.75,0.95]
# for i in range(len(qdata)):
#     a = data.quantile(qdata[i])
#     q[col[i]] = a.values