import biosteam as bst
import numpy as np
import chaospy
from math import sqrt, log, exp
from thermosteam import Stream
import pandas as pd
import thermosteam
from chaospy import distributions as shape
from system_setup.chemical import chem
import system_setup.chemical
from biorefinery.biorefinery_system import F, sys
from biorefinery.create_cellulosic_biodiesel_hmf import (
    create_oilcane_to_biodiesel_individual_oil_separation_microbial_vegetative)
import numpy as np
from biosteam.evaluation import Model, Metric
from system_setup.prices_yj import prices_per_Kg as price
from system_setup.prices_yj import set_price
from system_setup import process_settings_yj
from system_setup.lca_characterization_factors_yj import GWP_characterization_factors as GWP_factors
from system_setup.lca_characterization_factors_yj import set_GWPCF
from biosteam import preferences
from biosteam import report
from Uncertainty_analysis_functions.oilcane_composition_funcs import (
    set_lipid_fraction as set_oil_fraction,
    set_line_composition_parameters,
)
from Unit_functions.Evaluation_functions.cellulosic_ethanol_tea_yj import create_cellulosic_ethanol_tea as create_tea
from model_uncertainty_analysis import model,baseline_sys,baseline_metrics
from model_uncertainty_analysis import (get_MPSP_hmf,get_MPSP_biodiesel,get_MFPP, get_total_yield_furfural,get_total_yield_hmf,
    get_total_yield_biodiesel,get_total_yield_crude_glycerol,get_net_electricity_yield,get_purity_furfural,get_purity_hmf,
    get_vegetative_lipid_recovery,get_microbial_lipid_recovery,get_overall_TCI,get_overall_TDC,
    get_operating_cost,get_material_cost,get_system_heating_demand,get_system_cooling_demand, get_total_electricity_demand,
    get_natural_gas_demand,get_mlp_cost,get_mlp_electricity_consumption,get_brf_cooling_duty,get_brf_heating_duty,get_chp_cost,
    get_NPV,solve_IRR,get_net_GWP,get_hmf_purification_heating_demand_GWP,get_heating_demand_GWP_non_hmf_purification,
    get_electricity_demand_GWP,get_hmf_GWP_economic,get_hmf_GWP_energy,get_hmf_gwp_mass,get_hmf_GWP_all_displaced,
    get_biodiesel_GWP_economic,get_biodiesel_GWP_energy,get_biodiesel_gwp_mass)
metric_funcs=[get_MPSP_hmf,get_MPSP_biodiesel,get_MFPP, get_total_yield_furfural,get_total_yield_hmf,
    get_total_yield_biodiesel,get_total_yield_crude_glycerol,get_net_electricity_yield,get_purity_furfural,get_purity_hmf,
    get_vegetative_lipid_recovery,get_microbial_lipid_recovery,get_overall_TCI,get_overall_TDC,
    get_operating_cost,get_material_cost,get_system_heating_demand,get_system_cooling_demand, get_total_electricity_demand,
    get_natural_gas_demand,get_mlp_cost,get_mlp_electricity_consumption,get_brf_cooling_duty,get_brf_heating_duty,get_chp_cost,
    get_NPV,solve_IRR,get_net_GWP,get_hmf_purification_heating_demand_GWP,get_heating_demand_GWP_non_hmf_purification,
    get_electricity_demand_GWP,get_hmf_GWP_economic,get_hmf_GWP_energy,get_hmf_gwp_mass,get_hmf_GWP_all_displaced,
    get_biodiesel_GWP_economic,get_biodiesel_GWP_energy,get_biodiesel_gwp_mass]
feedstock = F.stream.oilcane
purified_furfural = F.stream.purified_furfural
purified_hmf = F.stream.purified_hmf
crude_glycerol = F.stream.crude_glycerol
def set_feedstock_lipid_content(lipid_content):
    feedstock_copy = bst.Stream('feedstock_copy')
    feedstock_copy.copy_like(feedstock)
    a = set_oil_fraction(lipid_fraction=lipid_content, stream=feedstock_copy,
                         FFA_fraction=0.1,
                         z_mass_carbs_baseline=0.149,
                         PL_fraction=0.1)
    feedstock.copy_like(a)
def adjust_veg_lipid_recovery(lipid_recovery):
    F.U402.isplit['Lipid'] = lipid_recovery
def set_rf_nf1_hmf(hmf_rf1):
    F.U301.rejection_factors['HMF'] = hmf_rf1
def set_rf_nf2(hmf_rf2):
    F.U302.rejection_factors['HMF'] = hmf_rf2
def set_rf_nf1_furfural(furfural_rf):
    F.U301.rejection_factors['Furfural'] = furfural_rf
def set_rf_nf1_furfural(furfural_rf2):
    F.U302.rejection_factors['Furfural'] = furfural_rf2
def set_NFmembrane_lifetime(lifetime):
    F.U301.equipment_lifetime = lifetime
    F.U302.equipment_lifetime = lifetime
def set_NFmembrane_cost(membrane_cost):
    F.U301.cost_items['Nanofiltration membrane'].cost = membrane_cost
    F.U302.cost_items['Nanofiltration membrane'].cost = membrane_cost
def set_lipid_production_conversion_glucose(conversion_glu):
    F.unit.R401.cofermentation[0].product_yield('TAG', basis='wt', product_yield=0.33 * conversion_glu)
    F.unit.R401.cofermentation[3].X = 0.999 - F.unit.R401.cofermentation[0].X
def set_lipid_production_conversion_xylose(conversion_xyl):
    F.unit.R401.cofermentation[1].product_yield('TAG', basis='wt', product_yield=0.34 * conversion_xyl)
    F.unit.R401.cofermentation[4].X = 0.999 - F.unit.R401.cofermentation[1].X
def set_lipid_production_conversion_arabinose(conversion_ara):
    F.unit.R401.cofermentation[2].product_yield('TAG', basis='wt', product_yield=0.34 * conversion_ara)
    F.unit.R401.cofermentation[5].X = 0.999 - F.unit.R401.cofermentation[2].X
def set_lipid_production_titer(t):
    F.unit.R401.titer = t
def set_feedstock_input_flow(Cap):
    feedstock.F_mass = Cap
def set_feedstock_price(feedstock_price):
    feedstock.price = feedstock_price
def set_furfural_price(furfural_price):
    purified_furfural.price = furfural_price
def set_hmf_price(hmf_price):
    purified_hmf.price = hmf_price
def set_crude_glycerol_price(crude_glycerol_price):
    crude_glycerol.price = crude_glycerol_price
def set_glycerine_price(glycerine_price):
    F.stream.pure_glycerine.price = glycerine_price
def set_natural_gas_price(natural_gas_price):
    F.stream.natural_gas.price = natural_gas_price
def set_electricity_price(electricity_price):
    bst.settings.electricity_price = electricity_price  # 2023 average price for industrial use (U.S. Energy Information Administration (EIA), 2023)
def set_pre_reactor_cost(pre_reactor_cost):
    F.unit.R201.cost_items['Pretreatment reactor system'].cost = pre_reactor_cost

parameter_funcs = [set_feedstock_lipid_content,adjust_veg_lipid_recovery,
                   set_rf_nf1_hmf,set_rf_nf2,set_rf_nf1_furfural,set_rf_nf1_furfural,set_NFmembrane_lifetime,set_NFmembrane_cost,
                   set_lipid_production_conversion_glucose,set_lipid_production_conversion_xylose,set_lipid_production_conversion_arabinose,
                   set_lipid_production_titer,set_feedstock_input_flow,set_feedstock_price,
                   set_furfural_price,set_hmf_price,set_crude_glycerol_price,set_glycerine_price,set_natural_gas_price,
                   set_electricity_price,set_pre_reactor_cost]

data = pd.read_excel('011520204_failed2.xlsx',usecols='B:V')
samples = data.values
cols = []
for i in model.parameters:
    cols.append(i)
for i in model.metrics:
    cols.append(i)
results = pd.DataFrame(columns=model.metrics)
failed_samples = pd.DataFrame(columns=model.parameters)

# for i in samples:
#     r=[]
#     for j in range(len(i)):
#         parameter_funcs[j](i[j])
#     try:
#         baseline_sys.simulate()
#         for k in range(len(metric_funcs)):
#             r.append(metric_funcs[k]())
#             results.loc[len(results)] = r
#         results.to_excel('01162024_results_failed2.xlsx')
#         print('Simulated',len(results))
#     except:
#         failed_samples.loc[len(failed_samples)] = i
#         failed_samples.to_excel('01162024_failed_failed2.xlsx')
#         print('Failed',len(failed_samples))