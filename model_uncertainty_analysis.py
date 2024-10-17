import biosteam as bst
import numpy as np
import pandas as pd
import argparse
from datetime import date
from chaospy import distributions as shape
from system_setup.chemical import chem
from biorefinery.biorefinery_system import F, sys
from biosteam.evaluation import Model, Metric
from system_setup.prices_yj import prices_per_Kg as price
from system_setup.prices_yj import set_price
from system_setup import process_settings_yj
from system_setup.lca_characterization_factors_yj import (
    GWP_characterization_factors as GWP_factors,
)
from system_setup.lca_characterization_factors_yj import set_GWPCF
from biosteam import preferences
from Uncertainty_analysis_functions.oilcane_composition_funcs import (
    set_lipid_fraction as set_oil_fraction,
)
from Unit_functions.Evaluation_functions.cellulosic_ethanol_tea_yj import (
    create_cellulosic_ethanol_tea as create_tea,
    foc_table as create_foc_table,
    capex_table as create_capex_table,
)
from Unit_functions.Evaluation_functions.lca_oilcane import create_lca_oilcane as create_lca

bst.settings.set_thermo(
    chem, cache=True
)  # takes all chemicals and let biosteam know the system_setup we want to use
process_settings_yj.load_process_settings()
bst.nbtutorial()
preferences.update(flow="kg/hr", T="degC", P="Pa", N=100, composition=True)
# sream volume F_vol = m3/hr, mass = kg/hr, mass/volume = kg/m3 = g/L
preferences.save()
GWP = "GWP100"
bst.settings.set_electricity_CF(
    GWP, GWP_factors["electricity"], basis="kWhr", units="kg*CO2e"
)
bst.settings.electricity_price = price[
    "Electricity"
]  # 2023 average price for industrial use (U.S. Energy Information Administration (EIA), 2023)
#Same as Humbird:
F.BT701.boiler_efficiency = 0.8
F.BT701.turbogen_efficiency = 0.85
# F.U602.Y_biomass = 0.02
baseline_sys = sys
baseline_sys.simulate()

UnitGroup = bst.process_tools.UnitGroup
area_names = {
    "100": "Feedstock processing",
    "200": "Pretreatment",
    "300": "Bioproducts recovery and purification",
    "400": "Microbial lipids production",
    "500": "Biodiesel production",
    "600": "Wastewater treatment",
    "700": "Coheat and power",  # Coheat & power system
    "800": "Utilities",  # (eg.,ADP,PWC,CT,CWP)
}
#rename utilities units to 800 group
utilities = [F.FWT701, F.CWP701, F.ADP701, F.CT701, F.CIP701, F.PWC701]
for i in utilities:
    bst.rename_unit(unit=i, area=800)
area_groups = UnitGroup.group_by_area(baseline_sys.units)
for i in area_groups:
    i.name = area_names[i.name]

def get_cost_area():
    cost_area = UnitGroup.df_from_groups(area_groups)
    return cost_area
# cost_area = get_cost_area()
############################################
# Only for conventional WWT:
set_GWPCF(F.stream.caustic, 'NaOH', 0.5)
set_price(F.stream.caustic, 'NaOH', 0.5)

# # High rate WWT with AnMBR only
# set_price(F.stream.citric_R602, "Citric_acid")
# set_price(
#     F.stream.bisulfite_R602, "Bisulfite"
# )  # price is for 38% solution, no need to dilute
#GWP#########################################################################################################
GWP = "GWP100"
# set_GWPCF(F.stream.PEG, "PEG")
# set_GWPCF(F.stream.MTBE, "MTBE")
set_GWPCF(F.stream.oilcane, "sugarcane")
set_GWPCF(F.stream.cellulase, "cellulase", dilution=0.05)
set_GWPCF(F.stream.H3PO4, "H3PO4")
set_GWPCF(F.stream.lime, "lime", dilution=0.046)  # Diluted with water
set_GWPCF(F.stream.FGD_lime, "lime", dilution=0.451)
# set_GWPCF(F.stream.hexane1, "hexane")
# set_GWPCF(F.stream.hexane2, "hexane")
set_GWPCF(F.stream.pure_glycerine, "pure-glycerol")
set_GWPCF(F.stream.HCl, "HCl")
set_GWPCF(F.stream.NaOH, "NaOH")
set_GWPCF(F.stream.methanol, "methanol")
set_GWPCF(F.stream.catalyst, "methanol catalyst mixture")
set_GWPCF(F.stream.makeup_process_water, "pwc_makeup_water")
set_GWPCF(F.stream.natural_gas, "CH4")
set_GWPCF(F.stream.cooling_tower_chemicals, "Cooling_tower_chemicals")
set_GWPCF(F.stream.boiler_chemicals, "Boiler_chemicals")
# Baseline prices#########################################################################################################
# Raw material prices:
set_price(F.stream.oilcane, "Oilcane")
set_price(F.stream.H3PO4, "H3PO4", dilution=0.588)  # 85% to 50%
set_price(F.stream.lime, "Lime", dilution=0.046)
set_price(F.stream.HCl, "HCl")
set_price(F.stream.polymer, "Polymer")
# set_price(F.stream.hexane1, "Hexane")
# set_price(F.stream.hexane2, "Hexane")
set_price(F.stream.pure_glycerine, "Pure_glycerine")
set_price(F.stream.makeup_process_water, "PWC_makeup_water")
F.stream.biodiesel_wash_water.price = 0
set_price(F.stream.cellulase, "Cellulase", dilution=0.05)
set_price(F.stream.NaOH, "NaOH")
set_price(F.stream.methanol, "Methanol")
set_price(F.stream.catalyst, "methanol catalyst mixture")
set_price(F.stream.cooling_tower_chemicals, "Cooling_tower_chemicals")
set_price(F.stream.boiler_chemicals, "Boiler_chemicals")
set_price(F.stream.FGD_lime, "Lime", dilution=0.451)
# set_price(F.stream.PEG, "PEG")
# set_price(F.stream.MTBE, "MTBE")
set_price(F.stream.NaOH, "NaOH")
set_price(F.stream.N2, "N2")  # only N2 for glycerolysis reactor
set_price(F.stream.natural_gas, "Natural_gas")
# coproduct prices:
set_price(F.stream.crude_glycerol, "Crude_glycerol", dilution=0.80)
set_price(F.stream.ash, "Ash_disposal")
F.stream.brine.price = 0
# product of interests:
# ratio = cellulosic_biomass_based_diesel_ratio()
biomass_diesel_price = price["Biomass based diesel"]
F.stream.biodiesel.price = biomass_diesel_price
# cellulosic_diesel_price = price['Cellulosic based diesel']
# F.stream.biodiesel.price = biomass_diesel_price* ratio[0] + cellulosic_diesel_price * ratio[1]
set_price(F.stream.purified_hmf, "Purified_HMF")
set_price(F.stream.purified_furfural, "Purified_furfural")

purified_hmf = F.stream.purified_hmf
purified_furfural = F.stream.purified_furfural
biodiesel = F.stream.biodiesel
crude_glycerol = F.stream.crude_glycerol
feedstock = F.stream.oilcane
products = [purified_hmf, biodiesel, purified_furfural, crude_glycerol]
product_streams_chemical_IDs = ['HMF', 'Biodiesel', 'Furfural', 'Glycerol']
product_stream_names = ['HMF', 'Biodiesel', 'Furfural', 'Crude glycerol'] #for allocation factor & GWP calculation

# TEA########################################################################################################## utility units
OSBL_unit = [F.FWT801, F.BT701, F.CWP801, F.ADP801, F.CT801, F.CIP801, F.PWC801]
for i in area_groups[5].units:
    if i not in OSBL_unit:
        OSBL_unit.append(i)
tea_baseline = create_tea(baseline_sys, OSBL_units=OSBL_unit)
# tea_baseline.operating_days = 300
tea_baseline.operating_days = 330
tea_baseline.IRR = 0.10
tea_baseline.duration = (2023, 2053)
tea_baseline.income_tax = 0.21
# tea.DPI =  tea.ISBL_installed_equipment_cost+tea.OSBL_installed_equipment_cost
# a = sum(i.installed_cost for i in ISBL units)
# tea.ISBL_installed_equipment_cost = a + a*(warehouse+site_development+additional piping)
# tea.OSBL_intalled_equipment_cost  = sum(i.installed_cost for i in OSBL units) =  cost_area.values[:,0][7] (=installed equipment cost of facilities)
# print(tea.get_cashflow_table())
sys_op_hours = baseline_sys.operating_hours = tea_baseline.operating_days * 24
oilcane_tea = baseline_sys.TEA


# baseline_sys.diagram(file='baseline_diagram',format='pdf')
# baseline_sys.save_report(file='baseline_sys_report.xlsx')

def get_MPSP_hmf():  # minimum product selling price of a product of interest assuming prices of other products are fixed
    MPSP_hmf = oilcane_tea.solve_price(purified_hmf)
    return MPSP_hmf


def get_MPSP_biodiesel():  # minimum product selling price of a product of interest assuming prices of other products are fixed
    MPSP_biodiesel = oilcane_tea.solve_price(biodiesel)
    return MPSP_biodiesel


def get_MFPP():  # maximum feedstock purchase price (oilcane) per MT
    MFPP = oilcane_tea.solve_price(feedstock) * 1000
    return MFPP


def get_annual_factor():
    annual_operating_hours = oilcane_tea.operating_days * 24
    return annual_operating_hours


def get_overall_TCI():
    TCI_sys = oilcane_tea.TCI / 1e6
    return TCI_sys


def get_overall_TDC():
    TDC_sys = oilcane_tea.TDC / 1e6
    return TDC_sys


def get_material_cost():
    material_cost = oilcane_tea.material_cost / 1e6
    return material_cost  # includes the money spent on ash disposal and gypsum sale


def get_operating_cost():
    annual_operating_cost = oilcane_tea.AOC / 1e6
    return annual_operating_cost


def get_NPV():
    NPV_sys = oilcane_tea.NPV
    return NPV_sys


# bst.plots.plot_unit_groups(area_groups, fraction=True)
def solve_IRR():
    IRR_sys = oilcane_tea.solve_IRR()
    return IRR_sys


def get_TDC():  # [MM$]
    TDC_sys = oilcane_tea.TDC / 1e6
    return TDC_sys

def get_purity_furfural():
    purity_furfural = purified_furfural.imass["Furfural"] / purified_furfural.F_mass
    return purity_furfural


def get_purity_hmf():
    purified_hmf = F.stream.purified_hmf
    purity_hmf = purified_hmf.imass["HMF"] / purified_hmf.F_mass
    return purity_hmf


def get_total_yield_hmf():
    yield_hmf = purified_hmf.F_mass * get_annual_factor() / 1e6
    return yield_hmf  # To get total yield in 10^6 Kg


def get_total_yield_furfural():
    yield_furfural = purified_furfural.F_mass * get_annual_factor() / 1e6
    return yield_furfural  # To get total yield in 10^6 Kg


def get_total_yield_biodiesel():
    yield_biodiesel = biodiesel.F_mass * get_annual_factor() / 1e6
    return yield_biodiesel  # To get total yield in 10^6 Kg


def get_total_yield_crude_glycerol():
    yield_glycerol = crude_glycerol.F_mass * get_annual_factor() / 1e6
    return yield_glycerol  # To get total yield in 10^6 Kg


def get_system_heating_demand(): #[GJ/hr]
    # total_heating_demand = -baseline_sys.operating_hours * sum([i.duty for i in F.BT801.heat_utilities if i.flow * i.duty > 0])
    total_heating_demand = sum([sum([i.duty for i in unit.heat_utilities if i.flow > 0 and i.duty > 0])for unit in baseline_sys.units])/1e6
    return total_heating_demand  # KJ/yr

def get_system_cooling_demand():
    # CT901_duty = -baseline_sys.operating_hours * sum([i.duty for i in F.CT901.heat_utilities if i.flow * i.duty < 0])
    # CWP901_duty = -baseline_sys.operating_hours * sum([i.duty for i in F.CWP901.heat_utilities if i.flow * i.duty < 0])
    # total_cooling_duty = -(CT901_duty+CWP901_duty)
    # CT = cooling tower, CWP = chilled water package
    total_cooling_duty = -sum([sum([i.duty for i in unit.heat_utilities if i.flow > 0 and i.duty < 0])for unit in baseline_sys.units])/1e6
    return total_cooling_duty  # KJ/yr


# LCA#########################################################################################################
# print(
#     report.lca_inventory_table(
#         systems=[baseline_sys],
#         key=GWP,
#         items=[purified_furfural,purified_hmf], # For including products without characterization factors
#     ))
#
# print(#TODO:Check lca_displacement_allocation_table
#     report.lca_displacement_allocation_table(
#         systems=[baseline_sys],
#         key=GWP,
#         items=[purified_furfural,purified_hmf],
#     ))

lca_baseline = create_lca(baseline_sys,
                          sys_chemcials=chem,
                          group_units=area_groups,
                          main_feedstock=feedstock,
                          products_for_allocation=products,
                          product_chemical_IDs = product_streams_chemical_IDs,
                          product_names = product_stream_names,
                          EOL_GWP=False,
                          )
# oilcane_lca =  baseline_sys.LCA
oilcane_lca = lca_baseline

#Test
lca_baseline_grave = create_lca(baseline_sys,
                          sys_chemcials=chem,
                          group_units=area_groups,
                          main_feedstock=feedstock,
                          products_for_allocation=products,
                          product_chemical_IDs = product_streams_chemical_IDs,
                          product_names = product_stream_names,
                          EOL_GWP=True,
                          )
oilcane_lca_grave = lca_baseline_grave

# for energy & economic allocation: assume all products have 0 GWP, use inputs to calculate overall GWP
# and then allocate it to each product based on allocation factors
# for displacement allocation: assume one product have 0 GWP, assume GWP of product = chemicals displaced by that product

#Product GWP
def get_hmf_GWP_energy():
    return oilcane_lca.HMF_GWP_energy
def get_hmf_GWP_economic():
    return oilcane_lca.HMF_GWP_economic
def get_hmf_GWP_displacement():
    return oilcane_lca.HMF_GWP_displacement
def get_biodiesel_GWP_energy():
    return oilcane_lca.Biodiesel_GWP_energy
def get_biodiesel_GWP_economic():
    return oilcane_lca.Biodiesel_GWP_economic

#Net GWP and breakdown
def get_net_GWP(): #million kg CO2-eq/yr
    return oilcane_lca.net_annual_GWP/1e6
# def get_hmf_purification_heating_demand_GWP(): #million kg CO2-eq/yr
#     return oilcane_lca.hmf_purification_heating_GWP/1e6
# def get_heating_demand_GWP_non_hmf_purification(): #million kg CO2-eq/yr
#     return oilcane_lca.non_hmf_purification_heating_GWP/1e6
# def get_electricity_demand_GWP(): #million kg CO2-eq/yr
#     return oilcane_lca.total_electricity_demand_GWP/1e6

#Utilitiesn usage
def get_natural_gas_demand():  # natural gas demand in 10^6 kg/yr
    natural_gas_demand = F.stream.natural_gas.F_mass * get_annual_factor()/1e6
    return natural_gas_demand

def get_total_electricity_demand(): #[MW]
    return oilcane_lca.net_electricity_demand/oilcane_lca.annual_factor/1000
def get_purchased_electricity(): #[MW]
    return oilcane_lca.purchase_electricity/oilcane_lca.annual_factor/1000

def get_net_electricity_yield():
    return oilcane_lca.electricity_yield

#Breakdown of GWP
def get_feedstock_GWP():
    gwp_breakdown = oilcane_lca.net_GWP_breakdown_alternative2
    return gwp_breakdown[0]/1e6
def get_other_material_GWP():
    gwp_breakdown = oilcane_lca.net_GWP_breakdown_alternative2
    return gwp_breakdown[1]/1e6
def get_total_direct_nonBT_emissions_GWP():
    gwp_breakdown = oilcane_lca.net_GWP_breakdown_alternative2
    return gwp_breakdown[2]/1e6
def get_brf_heating_GWP():
    gwp_breakdown = oilcane_lca.net_GWP_breakdown_alternative2
    return gwp_breakdown[3]/1e6
def get_remaining_heating_GWP():
    gwp_breakdown = oilcane_lca.net_GWP_breakdown_alternative2
    return gwp_breakdown[4]/1e6
def get_cooling_GWP():
    gwp_breakdown = oilcane_lca.net_GWP_breakdown_alternative2
    return gwp_breakdown[5]/1e6
def get_non_cooling_GWP():
    gwp_breakdown = oilcane_lca.net_GWP_breakdown_alternative2
    return gwp_breakdown[6]/1e6

def verify_GWP():
    gwp1 = oilcane_lca.net_annual_GWP/1e6
    gwp2 = oilcane_lca.net_annual_GWP_alternative2/1e6
    if abs(gwp1 - gwp2) > 0.001:
        diff = gwp1 - gwp2
    else:
        diff = 0
    return diff

def get_hmf_GWP_energy_grave():
    return oilcane_lca_grave.HMF_GWP_energy
def get_hmf_GWP_economic_grave():
    return oilcane_lca_grave.HMF_GWP_economic
def get_hmf_GWP_displacement_grave():
    return oilcane_lca_grave.HMF_GWP_displacement
def get_biodiesel_GWP_energy_grave():
    return oilcane_lca_grave.Biodiesel_GWP_energy
def get_biodiesel_GWP_economic_grave():
    return oilcane_lca_grave.Biodiesel_GWP_economic

def get_net_GWP_grave(): #million kg CO2-eq/yr
    return oilcane_lca_grave.net_annual_GWP/1e6
# Breakdown of GWP cradle to grave
def get_feedstock_GWP_grave():
    gwp_breakdown = oilcane_lca_grave.net_GWP_breakdown_alternative2
    return gwp_breakdown[0]/1e6
def get_other_material_GWP_grave():
    gwp_breakdown = oilcane_lca_grave.net_GWP_breakdown_alternative2
    return gwp_breakdown[1]/1e6
def get_total_direct_nonBT_emissions_GWP_grave():
    gwp_breakdown = oilcane_lca_grave.net_GWP_breakdown_alternative2
    return gwp_breakdown[2]/1e6
def get_brf_heating_GWP_grave():
    gwp_breakdown = oilcane_lca_grave.net_GWP_breakdown_alternative2
    return gwp_breakdown[3]/1e6
def get_remaining_heating_GWP_grave():
    gwp_breakdown = oilcane_lca_grave.net_GWP_breakdown_alternative2
    return gwp_breakdown[4]/1e6
def get_cooling_GWP_grave():
    gwp_breakdown = oilcane_lca_grave.net_GWP_breakdown_alternative2
    return gwp_breakdown[5]/1e6
def get_non_cooling_GWP_grave():
    gwp_breakdown = oilcane_lca_grave.net_GWP_breakdown_alternative2
    return gwp_breakdown[6]/1e6

def verify_GWP_grave():
    gwp1 = oilcane_lca_grave.net_annual_GWP/1e6
    gwp2 = oilcane_lca_grave.net_annual_GWP_alternative2/1e6
    if abs(gwp1 - gwp2) > 0.001:
        diff = gwp1 - gwp2
    else:
        diff = 0
    return diff

# Technical metrics for verification of parameters ###################################################################################################
#Lipid recovery
# def verify_microbial_lipid_yield():
#     current_X = F.unit.R401.cofermentation[0].X
#     current_yield = F.unit.R401.cofermentation[0].product_yield("TAG", basis="wt")
#     F.unit.R401.cofermentation[0]._chemicals.index('TAG')
#     if abs(current_X - current_yield/)
def get_vegetative_lipid_recovery():
    vegetative_lipid_recovery = (
        F.stream.vegetative_lipid.imass["Vegetative_lipid"]
        / F.S201.outs[0].imass["Vegetative_lipid"]
    )
    return vegetative_lipid_recovery*100
def get_microbial_lipid_recovery():
    microbial_lipid_recovery = (
        F.stream.microbial_lipid.imass["Microbial_lipid"]
        / F.T501.outs[0].imass["Microbial_lipid"]
    )
    return microbial_lipid_recovery*100
def get_vegetative_lipid_extraction_efficiency():
    return F.M505.lipid_extraction_efficiency()*100
def get_microbial_lipid_extraction_efficiency():
    return F.M502.lipid_extraction_efficiency()*100

def verify_vegetative_lipid_extraction_efficiency():
    diff = get_vegetative_lipid_extraction_efficiency() - 90
    if abs(diff) > 2:
        return diff
    else:
        return 0
def verify_microbial_lipid_extraction_efficiency():
    diff = get_microbial_lipid_extraction_efficiency() - 97
    if abs(diff) > 1:
        return diff
    else:
        return 0
#Fermentation
def get_actual_titer():
    return F.R401.get_titer()

def get_actual_productivity():
    return F.R401.titer/F.R401.tau

#Keep track of each section ###################################################################################################
feedstock_processing_units = area_groups[0].units
pretreatment_units = area_groups[1].units
bioproduct_recovery_units = area_groups[2].units
microbial_lipid_units = area_groups[3].units
biodiesel_production_units = area_groups[4].units
wwt_units = area_groups[5].units
coheat_power_units = area_groups[6].units
facilities_units = area_groups[7].units
#Installed equipment costs
def get_fp_cost():  # [MM$]
    cost_area = get_cost_area()
    return cost_area['Installed equipment cost [MM$]']['Feedstock processing']
def get_pt_cost():  # [MM$]
    cost_area = get_cost_area()
    return cost_area['Installed equipment cost [MM$]']['Pretreatment']
def get_brf_cost():  # [MM$]
    cost_area = get_cost_area()
    return cost_area['Installed equipment cost [MM$]']['Bioproducts recovery and purification']
def get_mlp_cost():  # [MM$]
    cost_area = get_cost_area()
    return cost_area['Installed equipment cost [MM$]']['Microbial lipids production']
def get_bp_cost():  # [MM$]
    cost_area = get_cost_area()
    return cost_area['Installed equipment cost [MM$]']['Biodiesel production']
def get_wwt_cost():  # [MM$]
    cost_area = get_cost_area()
    return cost_area['Installed equipment cost [MM$]']['Wastewater treatment']
def get_chp_cost():  # [MM$]
    cost_area = get_cost_area()
    return cost_area['Installed equipment cost [MM$]']['Coheat and power']
def get_facilities_cost():  # [MM$]
    cost_area = get_cost_area()
    return cost_area['Installed equipment cost [MM$]']['Utilities']

#Material costs
def get_fp_material_cost():  # [$/hr]
    cost_area = get_cost_area()
    return cost_area['Material cost [USD/hr]']['Feedstock processing']
def get_pt_material_cost():  # [$/hr]
    cost_area = get_cost_area()
    return cost_area['Material cost [USD/hr]']['Pretreatment']
def get_brf_material_cost():  # [$/hr]
    cost_area = get_cost_area()
    return cost_area['Material cost [USD/hr]']['Bioproducts recovery and purification']
def get_mlp_material_cost():  # [$/hr]
    cost_area = get_cost_area()
    return cost_area['Material cost [USD/hr]']['Microbial lipids production']
def get_bp_material_cost():  # [$/hr]
    cost_area = get_cost_area()
    return cost_area['Material cost [USD/hr]']['Biodiesel production']
def get_wwt_material_cost():  # [$/hr]
    cost_area = get_cost_area()
    return cost_area['Material cost [USD/hr]']['Wastewater treatment']
def get_chp_material_cost():  # [$/hr]
    cost_area = get_cost_area()
    return cost_area['Material cost [USD/hr]']['Coheat and power']
def get_facilities_material_cost():  # [$/hr]
    cost_area = get_cost_area()
    return cost_area['Material cost [USD/hr]']['Utilities']

#Heating demand
def get_fp_heating_demand():  # [GJ/hr]
    cost_area = get_cost_area()
    return cost_area['Heating duty [GJ/hr]']['Feedstock processing']
def get_pt_heating_demand():  # [GJ/hr]
    cost_area = get_cost_area()
    return cost_area['Heating duty [GJ/hr]']['Pretreatment']
def get_brf_heating_demand():  # [GJ/hr]
    cost_area = get_cost_area()
    return cost_area['Heating duty [GJ/hr]']['Bioproducts recovery and purification']
def get_mlp_heating_demand():  # [GJ/hr]
    cost_area = get_cost_area()
    return cost_area['Heating duty [GJ/hr]']['Microbial lipids production']
def get_bp_heating_demand():  # [GJ/hr]
    cost_area = get_cost_area()
    return cost_area['Heating duty [GJ/hr]']['Biodiesel production']
def get_wwt_heating_demand():  # [GJ/hr]
    cost_area = get_cost_area()
    return cost_area['Heating duty [GJ/hr]']['Wastewater treatment']
def get_chp_heating_demand():  # [GJ/hr]
    cost_area = get_cost_area()
    return cost_area['Heating duty [GJ/hr]']['Coheat and power']
def get_facilities_heating_demand():  # [GJ/hr]
    cost_area = get_cost_area()
    return cost_area['Heating duty [GJ/hr]']['Utilities']

#Cooling demand
def get_fp_cooling_demand():  # [GJ/hr]
    cost_area = get_cost_area()
    return cost_area['Cooling duty [GJ/hr]']['Feedstock processing']
def get_pt_cooling_demand():  # [GJ/hr]
    cost_area = get_cost_area()
    return cost_area['Cooling duty [GJ/hr]']['Pretreatment']
def get_brf_cooling_demand():  # [GJ/hr]
    cost_area = get_cost_area()
    return cost_area['Cooling duty [GJ/hr]']['Bioproducts recovery and purification']
def get_mlp_cooling_demand():  # [GJ/hr]
    cost_area = get_cost_area()
    return cost_area['Cooling duty [GJ/hr]']['Microbial lipids production']
def get_bp_cooling_demand():  # [GJ/hr]
    cost_area = get_cost_area()
    return cost_area['Cooling duty [GJ/hr]']['Biodiesel production']
def get_wwt_cooling_demand():  # [GJ/hr]
    cost_area = get_cost_area()
    return cost_area['Cooling duty [GJ/hr]']['Wastewater treatment']
def get_chp_cooling_demand():  # [GJ/hr]
    cost_area = get_cost_area()
    return cost_area['Cooling duty [GJ/hr]']['Coheat and power']
def get_facilities_cooling_demand():  # [GJ/hr]
    cost_area = get_cost_area()
    return cost_area['Cooling duty [GJ/hr]']['Utilities']

#Electricity consumption
def get_fp_electricity_consumption():  # [MW]
    cost_area = get_cost_area()
    return cost_area['Electricity consumption [MW]']['Feedstock processing']
def get_pt_electricity_consumption():  # [MW]
    cost_area = get_cost_area()
    return cost_area['Electricity consumption [MW]']['Pretreatment']
def get_brf_electricity_consumption():  # [MW]
    cost_area = get_cost_area()
    return cost_area['Electricity consumption [MW]']['Bioproducts recovery and purification']
def get_mlp_electricity_consumption():  # [MW]
    cost_area = get_cost_area()
    return cost_area['Electricity consumption [MW]']['Microbial lipids production']
def get_bp_electricity_consumption():  # [MW]
    cost_area = get_cost_area()
    return cost_area['Electricity consumption [MW]']['Biodiesel production']
def get_wwt_electricity_consumption():  # [MW]
    cost_area = get_cost_area()
    return cost_area['Electricity consumption [MW]']['Wastewater treatment']
def get_chp_electricity_consumption():  # [MW]
    cost_area = get_cost_area()
    return cost_area['Electricity consumption [MW]']['Coheat and power']
def get_facilities_electricity_consumption():  # [MW]
    cost_area = get_cost_area()
    return cost_area['Electricity consumption [MW]']['Utilities']

# # unit_by_groups = [feedstock_processing_units, pretreatment_units, bioproduct_recovery_units, microbial_lipid_units,
# #                   biodiesel_production_units, wwt_units, coheat_power_units, facilities_units]
#
# # def get_group_cost(): # [MM$]
# #     installed_cost_by_groups = []
# #     for i in unit_by_groups:
# #         installed_cost_by_groups.append(sum([j.installed_cost for j in i])/1e6)
# #     return installed_cost_by_groups
# # cost_by_groups = get_group_cost()
# def get_pretreatment_cost():  # [MM$]
#     return sum([j.installed_cost for j in pretreatment_units])/1e6
# def get_mlp_cost(): # [MM$]
#     return sum([j.installed_cost for j in microbial_lipid_units])/1e6
# def get_wwt_cost():  # [MM$]
#     return sum([j.installed_cost for j in wwt_units])/1e6
# def get_chp_cost():  # [MM$]
#     return sum([j.installed_cost for j in coheat_power_units])/1e6
#
# # def get_group_electricity_consumption():  # [kWh/yr]
# #     electricity_consumption_by_groups = []
# #     for i in unit_by_groups:
# #         electricity_consumption_by_groups.append(sum([j.power_utility.consumption for j in i]) * oilcane_lca.annual_factor)
# #     return electricity_consumption_by_groups
# # electricity_consumption_by_groups = get_group_electricity_consumption()
# def get_mlp_electricity_consumption():  # [MW]
#     return sum([j.power_utility.consumption for j in microbial_lipid_units])/1000
# def get_facilities_electricity_consumption():  # [MW]
#     return sum([j.power_utility.consumption for j in facilities_units])/1000
#
# # def get_group_heating_duty():  # [KJ/yr]
# #     heating_demand_by_groups = []
# #     for i in unit_by_groups:
# #         heating_demand_by_groups.append(sum([sum([j.duty for j in unit.heat_utilities if j.flow > 0 and j.duty > 0]) for unit in i]) * oilcane_lca.annual_factor)
# #     return heating_demand_by_groups
# # heating_duty_by_groups = get_group_heating_duty()
# def get_pretreatment_heating_duty():  # [GJ/hr]
#     return sum([sum([j.duty for j in unit.heat_utilities if j.flow > 0 and j.duty > 0]) for unit in pretreatment_units])/1e6
# def get_brf_heating_duty():  # [GJ/hr]
#     return sum([sum([j.duty for j in unit.heat_utilities if j.flow > 0 and j.duty > 0]) for unit in
#                 bioproduct_recovery_units])/1e6
#
# # def get_group_cooling_duty():  # [KJ/yr]
# #     cooling_demand_by_groups = []
# #     for i in unit_by_groups:
# #         cooling_demand_by_groups.append(-sum([sum([j.duty for j in unit.heat_utilities if j.flow > 0 and j.duty < 0])
# #                                               for unit in i]) * oilcane_lca.annual_factor)
# #     return cooling_demand_by_groups
# # cooling_duty_by_groups = get_group_cooling_duty()
#
# def get_brf_cooling_duty():  # [GJ/hr]
#     return -sum([sum([j.duty for j in unit.heat_utilities if j.flow > 0 and j.duty < 0])
#             for unit in bioproduct_recovery_units]) /1e6
# def get_mlp_cooling_duty():  # [GJ/hr]
#     return -sum([sum([j.duty for j in unit.heat_utilities if j.flow > 0 and j.duty < 0])
#             for unit in microbial_lipid_units]) /1e6

# ######################################################################################################################
metrics = (
    #Main economic indicators
    Metric("MFPP", get_MFPP, "$/MT"),
    Metric("HMF MPSP", get_MPSP_hmf, "$/kg"),
    Metric("Biodiesel MPSP", get_MPSP_biodiesel, "$/kg"),

    #Main environmental indicators
    Metric("Energy HMF GWP", get_hmf_GWP_energy, "kg CO2-eq/kg HMF"),
    Metric("Energy HMF GWP cradle to grave", get_hmf_GWP_energy_grave, "kg CO2-eq/kg HMF"),

    #Other environmental indicators (cradle to gate)
    Metric("Economic HMF GWP", get_hmf_GWP_economic, "kg CO2-eq/kg HMF"),
    Metric("HMF GWP all displaced", get_hmf_GWP_displacement, "kg CO2-eq/kg HMF"),
    Metric("Energy Biodiesel GWP", get_biodiesel_GWP_energy, "kg CO2-eq/kg biodiesel"),
    Metric("Economic Biodiesel GWP", get_biodiesel_GWP_economic, "kg CO2-eq/kg biodiesel"),

    #Other environmental indicators (cradle to grave)
    Metric("Economic HMF GWP cradle to grave", get_hmf_GWP_economic_grave, "kg CO2-eq/kg HMF"),
    Metric("HMF GWP all displaced cradle to grave", get_hmf_GWP_displacement_grave, "kg CO2-eq/kg HMF"),
    Metric("Energy Biodiesel GWP cradle to grave", get_biodiesel_GWP_energy_grave, "kg CO2-eq/kg biodiesel"),
    Metric("Economic Biodiesel GWP cradle to grave", get_biodiesel_GWP_economic_grave, "kg CO2-eq/kg biodiesel"),

    #Annual production of products
    Metric("Annual furfural production", get_total_yield_furfural, "million kg/yr"),
    Metric("Annual HMF production", get_total_yield_hmf, "million kg/yr"),
    Metric("Annual biodiesel production", get_total_yield_biodiesel, "million kg/yr"),
    Metric(
        "Annual crude glycerol production", get_total_yield_crude_glycerol, "million kg/yr"
    ),
    Metric("Annual electricity production (sold as coproduct)", get_net_electricity_yield, "kWh/yr"),

    #Utilities usage
    Metric("Total heating demand", get_system_heating_demand, "GJ/hr"),
    Metric("Total cooling demand", get_system_cooling_demand, "GJ/hr"),
    Metric("Total electricity usage", get_total_electricity_demand, "MW"),
    Metric("Total purchased electricity", get_purchased_electricity, "MW"),
    Metric("Total natural gas usage", get_natural_gas_demand, "million kg/yr"),

    #Technical metrics to make sure the system is working properly:
    Metric("Furfural purity", get_purity_furfural, "%"),
    Metric("HMF purity", get_purity_hmf, "%"),
    Metric("Vegetative lipid recovery", get_vegetative_lipid_recovery, "%"),
    Metric('Vegetative lipid extraction efficiency', get_vegetative_lipid_extraction_efficiency, "%"),
    Metric("Vegetative lipid extraction efficiency verification", verify_vegetative_lipid_extraction_efficiency, "%"),
    Metric("Microbial lipid recovery", get_microbial_lipid_recovery, "%"),
    Metric('Microbial lipid extraction efficiency', get_microbial_lipid_extraction_efficiency, "%"),
    Metric("Microbial lipid extraction efficiency verification", verify_microbial_lipid_extraction_efficiency, "%"),
    Metric("Actual lipid titer of fermentation", get_actual_titer, "g/L"),
    Metric("Actual lipid productivity of fermentation", get_actual_productivity, "g/L/hr"),
    # Metric("Burned bagasse fraction", get_burned_bagasse_fraction, "%"),

    #Economic metrics
    Metric("Total capital investment", get_overall_TCI, "MM$"),
    Metric("Total direct cost", get_overall_TDC, "MM$"),
    Metric("Annual operating cost", get_operating_cost, "MM$/yr"),
    Metric("Annual material cost", get_material_cost, "MM$/yr"),
    Metric("NPV", get_NPV, "$"),
    Metric("Internal rate of return", solve_IRR, "%"),

    #Environmental total & breakdown (cradle to gate)
    Metric("Annual system GWP", get_net_GWP, "million kg CO2-eq/yr"),
    #GWP breakdown
    Metric("GWP from feedstock", get_feedstock_GWP, "million kg CO2-eq/yr"),
    Metric("GWP from other material", get_other_material_GWP, "million kg CO2-eq/yr"),
    Metric("GWP from bioproducts recovery heating demand", get_brf_heating_GWP, "million kg CO2-eq/yr"),
    Metric("GWP from remaining heating demand", get_remaining_heating_GWP, "million kg CO2-eq/yr"),
    Metric("GWP from cooling demand", get_cooling_GWP, "million kg CO2-eq/yr"),
    Metric("GWP from non-cooling demand", get_non_cooling_GWP, "million kg CO2-eq/yr"),
    #GWP verification
    Metric("GWP verification", verify_GWP, "million kg CO2-eq/yr"),

    #Environmental total & breakdown (cradle to grave)
    Metric("Annual system GWP cradle to grave", get_net_GWP_grave, "million kg CO2-eq/yr"),
    #GWP breakdown
    Metric("GWP from feedstock cradle to grave", get_feedstock_GWP_grave, "million kg CO2-eq/yr"),
    Metric("GWP from other material cradle to grave", get_other_material_GWP_grave, "million kg CO2-eq/yr"),
    Metric("GWP from bioproducts recovery heating demand cradle to grave", get_brf_heating_GWP_grave, "million kg CO2-eq/yr"),
    Metric("GWP from remaining heating demand cradle to grave", get_remaining_heating_GWP_grave, "million kg CO2-eq/yr"),
    Metric("GWP from cooling demand cradle to grave", get_cooling_GWP_grave, "million kg CO2-eq/yr"),
    Metric("GWP from non-cooling demand cradle to grave", get_non_cooling_GWP_grave, "million kg CO2-eq/yr"),
    #GWP verification
    Metric("GWP verification cradle to grave", verify_GWP_grave, "million kg CO2-eq/yr"),

    #Metrics of sections:
    Metric("Feedstock processing IEC", get_fp_cost, "MM$"),
    Metric("Feedstock processing MC", get_fp_material_cost, "$/hr"),
    Metric("Feedstock processing heating demand", get_fp_heating_demand, "GJ/hr"),
    Metric("Feedstock processing cooling demand", get_fp_cooling_demand, "GJ/hr"),
    Metric("Feedstock processing electricity consumption", get_fp_electricity_consumption, "MW"),

    Metric("Pretreatment IEC", get_pt_cost, "MM$"),
    Metric("Pretreatment MC", get_pt_material_cost, "$/hr"),
    Metric("Pretreatment heating demand", get_pt_heating_demand, "GJ/hr"),
    Metric("Pretreatment cooling demand", get_pt_cooling_demand, "GJ/hr"),
    Metric("Pretreatment electricity consumption", get_pt_electricity_consumption, "MW"),

    Metric("Brf IEC", get_brf_cost, "MM$"),
    Metric("Brf MC", get_brf_material_cost, "$/hr"),
    Metric("Brf heating demand", get_brf_heating_demand, "GJ/hr"),
    Metric("Brf cooling demand", get_brf_cooling_demand, "GJ/hr"),
    Metric("Brf electricity consumption", get_brf_electricity_consumption, "MW"),

    Metric("Mlp IEC", get_mlp_cost, "MM$"),
    Metric("Mlp MC", get_mlp_material_cost, "$/hr"),
    Metric("Mlp heating demand", get_mlp_heating_demand, "GJ/hr"),
    Metric("Mlp cooling demand", get_mlp_cooling_demand, "GJ/hr"),
    Metric("Mlp electricity consumption", get_mlp_electricity_consumption, "MW"),

    Metric("Bp IEC", get_bp_cost, "MM$"),
    Metric("Bp MC", get_bp_material_cost, "$/hr"),
    Metric("Bp heating demand", get_bp_heating_demand, "GJ/hr"),
    Metric("Bp cooling demand", get_bp_cooling_demand, "GJ/hr"),
    Metric("Bp electricity consumption", get_bp_electricity_consumption, "MW"),

    Metric("Wwt IEC", get_wwt_cost, "MM$"),
    Metric("Wwt MC", get_wwt_material_cost, "$/hr"),
    Metric("Wwt heating demand", get_wwt_heating_demand, "GJ/hr"),
    Metric("Wwt cooling demand", get_wwt_cooling_demand, "GJ/hr"),
    Metric("Wwt electricity consumption", get_wwt_electricity_consumption, "MW"),

    Metric("Chp IEC", get_chp_cost, "MM$"),
    Metric("Chp MC", get_chp_material_cost, "$/hr"),
    Metric("Chp heating demand", get_chp_heating_demand, "GJ/hr"),
    Metric("Chp cooling demand", get_chp_cooling_demand, "GJ/hr"),
    Metric("Chp electricity consumption", get_chp_electricity_consumption, "MW"),

    Metric("Facilities IEC", get_facilities_cost, "MM$"),
    Metric("Facilities MC", get_facilities_material_cost, "$/hr"),
    Metric("Facilities heating demand", get_facilities_heating_demand, "GJ/hr"),
    Metric("Facilities cooling demand", get_facilities_cooling_demand, "GJ/hr"),
    Metric("Facilities electricity consumption", get_facilities_electricity_consumption, "MW"),
)

model = Model(
    baseline_sys,
    metrics,
    # retry_evaluation=False,
    exception_hook="warn",
    # raise immediately stop once simulation fails, then need to restart
    # warn leaves the failed simulation blank and continues with other simulations
)

# parser = argparse.ArgumentParser()
# # parser.add_argument('job_number', type=int)
# parser.add_argument('data_path')
# parser.add_argument('filename')
#
# args = parser.parse_args()
#
# # N_samples = args.job_number # 2000 samples are enough
# data_path = args.data_path
# filename = args.filename


parameter_distribution = pd.read_excel('06212024_N2000_parameter_distribution.xlsx',
                                      usecols='C:G',sheet_name='Sheet1') #parameter name as index
parameter_distribution.index = parameter_distribution['Name']
# Uncertainty analysis needs at least five parameters
# Vegetative lipid related
# Lipid fraction in dry oilcane feedstock
flc_dist = parameter_distribution.loc["Feedstock lipid content"]
# lb_flc = 0.05
# ub_flc = 0.15
@model.parameter(
    name=flc_dist['Name'],
    element=F.stream.oilcane,
    kind="coupled",
    units="%",
    distribution=shape.Uniform(flc_dist['lower'], flc_dist['upper']),
    baseline=0.05,
)
def set_feedstock_lipid_content(lipid_content):
    feedstock_copy = bst.Stream("feedstock_copy")
    feedstock_copy.copy_like(feedstock)
    a = set_oil_fraction(
        lipid_fraction=lipid_content,
        stream=feedstock_copy,
        FFA_fraction=0.1,
        z_mass_carbs_baseline=0.149,
        PL_fraction=0.1,
    )
    feedstock.copy_like(a)


# Vegtative lipid recovery:
# 1. assume 5% loss of lipids at juicing (U104)
# 2. assume all lipids loss after juicing is after saccharification in separation of lignin and hydrolysate (lost part stay in hydrolystae)
# f.U402.isplit['Lipid'] = fraction of lipids that go to lignin for extraction
# lb_lr = 0.7*0.75
# up_lr = 0.7*1.25

vlr_dist = parameter_distribution.loc["Vegetative lipid recovery after pretreatment & saccharification"]
@model.parameter(
    name=vlr_dist['Name'],
    element=F.U402,
    kind="coupled",
    units="%",
    distribution=shape.Uniform(vlr_dist['lower'], vlr_dist['upper']),
    baseline=0.7,#70% of lipid goes to saccharified biomass for lipid extraction, 30% loss
)
def adjust_veg_lipid_recovery(lipid_recovery):
    F.U402.isplit["Lipid"] = lipid_recovery


# Nanofiltration related
# HMF retentions in both nanofiltrations
hmf_nf1_dist = parameter_distribution.loc["HMF retention in 1st nanofiltration"]
# lb_hmf_nf1 = 0.4072 * 0.5
# ub_hmf_nf1 = 0.4072

@model.parameter(
    name=hmf_nf1_dist['Name'],
    element=F.U301,
    kind="coupled",
    units="%",
    distribution=shape.Uniform(hmf_nf1_dist['lower'], hmf_nf1_dist['upper']),
    baseline=0.4072,
)
def set_rf_nf1_hmf(hmf_rf1):
    F.U301.rejection_factors["HMF"] = hmf_rf1

hmf_nf2_dist = parameter_distribution.loc["HMF retention in 2nd nanofiltration"]
lb_hmf_nf2 = 0.4755 * 0.5
ub_hmf_nf2 = 0.4755


@model.parameter(
    name=hmf_nf2_dist['Name'],
    element=F.U302,
    kind="coupled",
    units="%",
    distribution=shape.Uniform(hmf_nf2_dist['lower'], hmf_nf2_dist['upper']),
    baseline=0.4755,
)
def set_rf_nf2(hmf_rf2):
    F.U302.rejection_factors["HMF"] = hmf_rf2


# Furfural retentions in both nanofiltrations
furfural_nf1_dist = parameter_distribution.loc["Furfural retention in 1st nanofiltration"]
lb_furfural_nf1 = 0.5959* 0.5
ub_furfural_nf1 = 0.5959


@model.parameter(
    name=furfural_nf1_dist['Name'],
    element=F.U301,
    kind="coupled",
    units="%",
    distribution=shape.Uniform(furfural_nf1_dist['lower'], furfural_nf1_dist['upper']),
    baseline=0.5959,
)
def set_rf_nf1_furfural(furfural_rf):
    F.U301.rejection_factors["Furfural"] = furfural_rf

furfural_nf2_dist = parameter_distribution.loc["Furfural retention in 2nd nanofiltration"]
lb_furfural_nf2 = 0.7882 * 0.75
ub_furfural_nf2 = 0.7882 * 1.25


@model.parameter(
    name=furfural_nf2_dist['Name'],
    element=F.U302,
    kind="coupled",
    units="%",
    distribution=shape.Uniform(lb_furfural_nf2, ub_furfural_nf2),
    baseline=0.7882,
)
def set_rf_nf1_furfural(furfural_rf2):
    F.U302.rejection_factors["Furfural"] = furfural_rf2


# NF: membrane price, membrane lifetime
membrane_splitter = F.unit.S301
nml_dist = parameter_distribution.loc["Nanofiltration membrane lifetime"]
lb_nml = 0.5
ub_nml = 5.0

@model.parameter(name=nml_dist['Name'],
    element=(membrane_splitter),
    kind="coupled",
    units="yr",
    distribution=shape.Uniform(nml_dist['lower'], nml_dist['upper']),
    baseline=2.5,
)
def set_NFmembrane_lifetime(lifetime):
    membrane_splitter.nanofiltration_membrane_lifetime = lifetime
    N1_area = F.U301.design_results['Total filtration area'] #m2
    N2_area = F.U302.design_results['Total filtration area'] #m2
    membrane_thickness = 280/1e9 #nm, thickness of polyamide tfc membrane

    N1_vol = (N1_area*membrane_thickness)/(lifetime * get_annual_factor()) #m3/hr
    N2_vol = (N2_area*membrane_thickness)/(lifetime * get_annual_factor()) #m3/hr

    membrane_splitter.outs[0].ivol['polyamide'] = N1_vol
    membrane_splitter.outs[1].ivol['polyamide'] = N2_vol

# Membrane purchase cost
nmc_dist = parameter_distribution.loc["Nanofiltration membrane cost"]
membrane_cost_baseline = 20.8*803.2/394 #$/m2, 2023 price
lb_nmc = membrane_cost_baseline * 0.75
ub_nmc = membrane_cost_baseline * 1.25

@model.parameter(
    name=nmc_dist['Name'],
    element=(F.stream.nanofiltration_membrane),
    kind="cost",
    units="$/m2",
    distribution=shape.Uniform(nmc_dist['lower'], nmc_dist['upper']),
    baseline=membrane_cost_baseline,
)
def set_NFmembrane_cost(membrane_cost):
    F.stream.nanofiltration_membrane.price = membrane_cost/(280/1e9) #$/m3

# Microbial lipids related
# Microbial lipids fermentation: lipid production reaction conversion%
# lb_mly = 0.4  # 40% of theoretical yield
# # ub_mly = 0.18 / 0.33 *1.25
# ub_mly = 0.9  # 55% of theoretical yield
glu_lipid_yield_dist = parameter_distribution.loc["Microbial lipid yield from glucose"]
@model.parameter(
    name=glu_lipid_yield_dist['Name'],
    element=F.R401,
    units="%",
    kind="coupled",
    # distribution=shape.Uniform(0.18/0.33*0.75, 0.18/0.33*1.25),
    distribution=shape.Uniform(glu_lipid_yield_dist['lower'], glu_lipid_yield_dist['upper']),
    baseline=0.18 / 0.33,
)
def set_lipid_production_conversion_glucose(conversion_glu):
    F.unit.R401.cofermentation[0].product_yield(
        "TAG", basis="wt", product_yield=0.33 * conversion_glu
    )
    F.unit.R401.cofermentation[2].X = 0.999 - F.unit.R401.cofermentation[0].X

xyl_lipid_yield_dist = parameter_distribution.loc["Microbial lipid yield from xylose"]
@model.parameter(
    name=xyl_lipid_yield_dist['Name'],
    element=F.R401,
    units="%",
    kind="coupled",
    # distribution=shape.Uniform(0.18/0.33*0.75, 0.18/0.33*1.25),
    distribution=shape.Uniform(xyl_lipid_yield_dist['lower'], xyl_lipid_yield_dist['upper']),
    baseline=0.18 / 0.33,
)
def set_lipid_production_conversion_xylose(conversion_xyl):
    F.unit.R401.cofermentation[1].product_yield(
        "TAG", basis="wt", product_yield=0.34 * conversion_xyl
    )
    F.unit.R401.cofermentation[3].X = 0.999 - F.unit.R401.cofermentation[1].X

productivity_base = F.unit.R401.productivity
# lb_pro = 0.17 #from hydrolysate study
# lb_pro = productivity_base * 0.75
# ub_pro = productivity_base * 1.25
# try:
productivity_dist = parameter_distribution.loc["Microbial lipid productivity"]
@model.parameter(
    name=productivity_dist['Name'],
    # name='Microbial lipid productivity',
    element=F.R401,
    units='g/L/hr',
    kind='coupled',
    # distribution=shape.Uniform(productivity_base * 0.75, productivity_base * 1.25),
    distribution=shape.Uniform(productivity_dist['lower'], productivity_dist['upper']),
    baseline=productivity_base,
)
def set_lipid_production_productivity(prod):
    F.unit.R401.productivity = prod

# except:
#     pass

titer_base = F.unit.R401.titer
# lb_titer = titer_base *0.75
# ub_titer = titer_base * 1.25
titer_dist = parameter_distribution.loc["Microbial lipid titer"]
@model.parameter(
    name=titer_dist['Name'],
    element=F.R401,
    units="g/L",
    kind="coupled",
    # distribution=shape.Uniform(titer_base * 0.75, titer_base * 1.25),
    distribution=shape.Uniform(titer_dist['lower'], titer_dist['upper']),
    baseline=titer_base,
)
def set_lipid_production_titer(t):
    F.unit.R401.titer = t

# # Decided to exclude:
# lb_te = 0.8
# ub_te = 0.9
# @model.parameter(name='Turbogen efficiency',
#                      element=F.BT701,
#                      kind='coupled',
#                      distribution=shape.Uniform(lb_te, ub_te),
#                  baseline=0.85)
# def set_tubeff(X_tubeff):
#     F.BT701.turbogenerator_efficiency = X_tubeff
# #
# lb_beff = 0.75
# ub_beff = 0.88
# @model.parameter(name='Boiler efficiency',
#                      element=F.BT701,
#                      kind='coupled',
#                      distribution=shape.Uniform(lb_beff, ub_beff),
#                  baseline=0.8)
# def set_beff(X_beff):
#     F.BT701.boiler_efficiency = X_beff

# # Common parameters
# lb_cap = 3.33e+05 * 0.75
# ub_cap = 3.33e+05 * 1.25
# @model.parameter(name='Feedstock capacity',
#                  element=F.stream.oilcane,
#                  kind='coupled', units='kg/hr',
#                  distribution=shape.Uniform(lb_cap, ub_cap),
#                  baseline=3.33e+05,
#                  )
# def set_feedstock_input_flow(Cap):
#     feedstock.F_mass = Cap

# feedstock prices, sales prices
fp_dist = parameter_distribution.loc["Oilcane price"]
b_fp = feedstock.price
# lb_fp = feedstock.price * 0.75
# ub_fp = feedstock.price * 1.25  # Maximum price

@model.parameter(
    name=fp_dist['Name'],
    element=feedstock,
    kind="isolated",
    units="$/kg",
    distribution=shape.Uniform(fp_dist['lower'], fp_dist['upper']),
    baseline=b_fp,
)
def set_feedstock_price(feedstock_price):
    feedstock.price = feedstock_price


# sales prices
furprice_dist = parameter_distribution.loc["Furfural price"]
b_furprice = purified_furfural.price
# lb_furprice = purified_furfural.price * 0.75
# ub_furprice = purified_furfural.price * 1.25


@model.parameter(
    name=furprice_dist['Name'],
    element=purified_furfural,
    kind="isolated",
    units="$/kg",
    distribution=shape.Uniform(furprice_dist['lower'], furprice_dist['upper']),
    baseline=b_furprice,
)
def set_furfural_price(furfural_price):
    purified_furfural.price = furfural_price

hmfprice_dist = parameter_distribution.loc["HMF price"]
b_hmfprice = purified_hmf.price
# lb_hmfprice = purified_hmf.price * 0.9
# ub_hmfprice = purified_hmf.price * 1.1


@model.parameter(
    name=hmfprice_dist['Name'],
    element=purified_hmf,
    kind="isolated",
    units="$/kg",
    distribution=shape.Uniform(hmfprice_dist['lower'], hmfprice_dist['upper']),
    baseline=b_hmfprice,
)
def set_hmf_price(hmf_price):
    purified_hmf.price = hmf_price

bioprice_dist = parameter_distribution.loc["Biodiesel price"]
b_bioprice = biodiesel.price
# lb_bioprice = biodiesel.price * 0.75
# ub_bioprice = biodiesel.price * 1.25


@model.parameter(name=bioprice_dist['Name'],
                 element=biodiesel, kind='isolated', units='$/kg',
                 distribution=shape.Uniform(bioprice_dist['lower'], bioprice_dist['upper']),
                 baseline=b_bioprice, )
def set_biodiesel_price(biodiesel_price):
    biodiesel.price = biodiesel_price

cgp_dist = parameter_distribution.loc["Crude glycerol price"]
b_cgp = crude_glycerol.price
# lb_cgp = crude_glycerol.price * 0.75
# ub_cgp = crude_glycerol.price * 1.25


@model.parameter(
    name=cgp_dist['Name'],
    element=crude_glycerol,
    kind="isolated",
    units="$/kg",
    distribution=shape.Uniform(cgp_dist['lower'], cgp_dist['upper']),
    baseline=b_cgp,
)
def set_crude_glycerol_price(crude_glycerol_price):
    crude_glycerol.price = crude_glycerol_price


# raw material prices
caustic_dist = parameter_distribution.loc["Caustic price"]
b_caustic = F.stream.caustic.price
# lb_caustic = b_caustic*0.75
# ub_caustic= b_caustic*1.25
@model.parameter(name=caustic_dist['Name'],
                  element=F.stream.caustic, kind='isolated', units='$/kg',
                  distribution=shape.Uniform(caustic_dist['lower'], caustic_dist['upper']),
                 baseline=b_caustic,)
def set_caustic_price(caustic_price):
    F.stream.caustic.price = caustic_price

pg_dist = parameter_distribution.loc['Pure glycerine price']
b_pg = F.stream.pure_glycerine.price
# lb_pg = F.stream.pure_glycerine.price * 0.75
# ub_pg = F.stream.pure_glycerine.price * 1.25


@model.parameter(
    name=pg_dist['Name'],
    element=F.stream.pure_glycerine,
    kind="isolated",
    units="$/kg",
    distribution=shape.Uniform(pg_dist['lower'], pg_dist['upper']),
    baseline=b_pg,
)
def set_glycerine_price(glycerine_price):
    F.stream.pure_glycerine.price = glycerine_price

ng_dist = parameter_distribution.loc["Natural gas price"]
b_ng = F.stream.natural_gas.price
# lb_ng = F.stream.natural_gas.price * 0.75
# ub_ng = F.stream.natural_gas.price * 1.25


@model.parameter(
    name=ng_dist['Name'],
    element=F.stream.natural_gas,
    kind="isolated",
    units="$/kg",
    distribution=shape.Uniform(ng_dist['lower'], ng_dist['upper']),
    baseline=b_ng,
)
def set_natural_gas_price(natural_gas_price):
    F.stream.natural_gas.price = natural_gas_price

elecprice_dist = parameter_distribution.loc["Electricity price"]
b_elecprice = price["Electricity"]
# lb_elecprice = price["Electricity"] * 0.75
# ub_elecprice = price["Electricity"] * 1.25

#include since excss electricity is sold at this price too
@model.parameter(
    name=elecprice_dist['Name'],
    element=baseline_sys,
    kind="isolated",
    units="$/kWh",
    distribution=shape.Uniform(elecprice_dist['lower'], elecprice_dist['upper']),
    baseline=b_elecprice,
)
def set_electricity_price(electricity_price):
    bst.settings.electricity_price = electricity_price  # 2023 average price for industrial use (U.S. Energy Information Administration (EIA), 2023)

pre_dist = parameter_distribution.loc["Pretreatment reactor system cost"]
pre_base = F.unit.R201.cost_items["Pretreatment reactor system"].cost / 1000000
# lb_pre = pre_base * 0.75
# ub_pre = pre_base * 1.25


@model.parameter(
    name=pre_dist['Name'],
    element=F.unit.R201,
    kind="cost",
    units="MM$",
    distribution=shape.Uniform(pre_dist['lower'], pre_dist['upper']),
    baseline=pre_base,
)
def set_pre_reactor_cost(pre_reactor_cost):
    F.unit.R201.cost_items["Pretreatment reactor system"].cost = (
        pre_reactor_cost * 1000000
    )

#Always run baseline_metrics so that baseline_sys is at baseline parameters
baseline_metrics = model.metrics_at_baseline()
# baseline_metrics.to_excel('06132024_baseline_metrics_burnbagasse_satisfyTrue_N500.xlsx')
# model.get_distribution_summary()['Uniform'].to_excel('06172024_test2_parameter_distribution.xlsx')

#single point sensitvity
#baseline, lower, upper = model.single_point_sensitivity()
# mfpp,hmf_mpsp,biodiesel_mpsp,energy_hmf_gwp = model.metrics[0:4]
# metric_index = mfpp.index
# index = [i.describe(distribution=False) # Instead of displaying distribution, it displays lower, baseline, and upper values
#          for i in model.parameters]
# sp_plot = bst.plots.plot_single_point_sensitivity(100 * baseline[metric_index],
#                                         100 * lower[metric_index],
#                                         100 * upper[metric_index],
#                                         name='MFPP [$/MT]',
#                                         index=index)
# sp_plot[0].set_figheight(20)
# sp_plot[0].set_figwidth(50)
# sp_plot[0].show()

# data = pd.read_excel('05242024_samples_biomass_yield_5.xlsx',usecols='B:U')
# all = data.values
# model.load_samples(np.array([all[2]]))
# model.evaluate()
# baseline_metrics.to_excel('05242024_hxn_qmin_1e3.xlsx')
# baseline_sys.diagram(file='05242024_baseline_diagram_hxn',format='pdf')
# baseline_sys.save_report(file='05242024_baseline_sys_report_recycle_new_env.xlsx')
# #baseline_metrics.to_excel('01132024_baseline_metrics.xlsx') #y_yield = 0.8 for AnMBR
# N_samples = 2  # 1000 samples are enough
# rule = 'L'  # For Latin-Hypercube sampling
# np.random.seed(1234)  # For consistent results
# samples = model.sample(N_samples, rule)
# model.load_samples(samples)
# model.evaluate(
#     notify=2  # Also print elapsed Time after 100 simulations
# )
# # model.show()
# model.table.to_excel('05242024_results_biomass_yield_5.xlsx')

# N_samples = 5  # 1000 samples are enough
# rule = 'L'  # For Latin-Hypercube sampling
# np.random.seed(1234)  # For consistent results
# samples = model.sample(N_samples, rule)
# # pd.DataFrame(samples).to_excel('01192024_N1000_samples.xlsx')

# # data = pd.read_excel('01192024_failed_results_0.xlsx',usecols="B:U")
# # samples = data.values
# #
# # cols = []
# # for i in model.parameters:
# #     cols.append(i)
# # for i in model.metrics:
# #     cols.append(i)
# # results = pd.DataFrame(columns=cols)
# # failed_samples = pd.DataFrame(columns=model.parameters)
# #
# parser = argparse.ArgumentParser()
# parser.add_argument('job_number', type=int)
# JOB_ID = parser.parse_args().job_number
#
# RUNS_PER_JOB = 20
# #
# for i in tqdm(samples[RUNS_PER_JOB * JOB_ID : RUNS_PER_JOB * (JOB_ID + 1)]):
#     s = np.array(([i]))
#     # s = np.array(([samples[JOB_NUMBER]))
#     model.load_samples(s)
#     try:
#         model.evaluate()
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
#
# # results.to_excel('01132024_results_ybiomass8_R2.xlsx')
# # failed_samples.to_excel('01132024_failed_samples_ybiomass8_R2.xlsx')

# model.load_samples(samples)
# model.evaluate(
#     notify=2  # Also print elapsed Time after 100 simulations
# )
# # model.show()

 # a1 = model.table
# # model.table.to_excel('01152024_results_a_R1.xlsx')
# # s = a1[a1.columns[0:21]].values
# # failed_samples.loc[len(failed_samples)] = s[24]
# # model.load_samples(s[24+1:])
#
# # To combine data from different files
# # a = pd.read_excel('01182024_results_0.xlsx')
# # results = pd.DataFrame(columns=a.columns)
# # for i in range(50):
# #     ri = pd.read_excel(f"01182024_results_{i}.xlsx")
# #     r = pd.concat([results,ri])
# #     results = r

# model.table.to_excel('05222024_results_test.xlsx')
# # df_rho, df_p = model.spearman_r()
# # df_rho.to_excel('rho.xlsx')
# # df_p.to_excel('p.xlsx')
# # Spearman rho value: 0.5 - 1.0 or -0.5 to -1.0
# # no significant correlation at -0.5 to 0.5
# # but anything above abs(0.1) is vaguely impacting the results
# # p value smaller than 0.05 means the parameter is significant
# # A positive correlation means that as one variable increases, the other variable also tends to increase.
# # A negative correlation signifies that as one variable increases, the other tends to decrease.

# # #########################################################################################################
# # Generate table for supplementary materials: FOC, CAPEX, VOC, LCI, AF
# # cost_area = UnitGroup.df_from_groups(area_groups)
# # # # #fixed operating cost table
# foc = create_foc_table(baseline_sys.TEA)
# capex = create_capex_table(baseline_sys.TEA)
# foc.to_excel('N2000_current_06242024/baseline/08192024_foc.xlsx')
# capex.to_excel('N2000_current_06242024/baseline/08192024_capex.xlsx')
# # voc = bst.report.voc_table(systems = baseline_sys,
# # product_IDs = [F.stream.purified_hmf.ID,F.stream.purified_furfural.ID,F.stream.biodiesel.ID,
# #                F.stream.crude_glycerol.ID,F.stream.ash.ID])
# # voc.to_excel('01132024_voc.xlsx')
# # #life cycle inventory table
# lci = bst.report.lca_inventory_table(
#         systems=[baseline_sys],
#         key=GWP,items=[purified_hmf,purified_furfural,biodiesel,crude_glycerol])
# feeds = sorted({i.ID for i in sum([baseline_sys.feeds], []) if GWP in i.characterization_factors})
# cf = [F.stream[i].characterization_factors[GWP] for i in feeds]
# cf += [GWP_factors["electricity"],0,0,0,0]
# lci.insert(1, 'Characterization factor [kg CO2 eq./kg input]', cf)

# lci.to_excel('N2000_current_06242024/baseline/08192024_capex.xlsx.xlsx')

# #allocation factors
# af = pd.DataFrame(columns =oilcane_lca.allocation_factor_energy.keys())
# af.loc['Energy'] =  oilcane_lca.allocation_factor_energy.values()
# af.loc['Economic'] =  oilcane_lca.allocation_factor_economic.values()
# af.drop(columns=['Electricity'],inplace=True)
# af.to_excel('N2000_current_06242024/baseline/08192024_allocation_factor.xlsx')

# # VOC table
# def reformat(name):
#     name = name.replace('_', ' ')
#     if name.islower(): name = name.capitalize()
#     return name
#
# voc_ID = []
# voc_price=[]
# voc_cost=[]
# feeds = sorted({i.ID for i in baseline_sys.feeds if baseline_sys.has_market_value(i)})
# system_heat_utilities = bst.HeatUtility.sum_by_agent(baseline_sys.heat_utilities)
# product_IDs = ['purified_hmf', 'biodiesel', 'purified_furfural']
# coproducts = sorted({i.ID for i in baseline_sys.products if baseline_sys.has_market_value(i) and i.ID not in product_IDs})
# voc_ID.extend(feeds)
# system_heat_utilities = bst.HeatUtility.sum_by_agent(baseline_sys.heat_utilities)
# heating_agents = sorted(set(sum([[i.agent.ID for i in system_heat_utilities if i.cost and i.flow * i.duty > 0. and abs(i.flow) > 1e-6]], [])))
# cooling_agents = sorted(set(sum([[i.agent.ID for i in system_heat_utilities if i.cost and i.flow * i.duty < 0. and abs(i.flow) > 1e-6]], [])))
# voc_ID.extend(heating_agents)
# voc_ID.extend(cooling_agents)
# voc_ID.extend(coproducts)
# voc_ID.extend(product_IDs)
# voc_ID_reformated = [reformat(i) for i in voc_ID]
# for i in voc_ID:
#     if i not in cooling_agents:
#         s = F.stream[i]
#         voc_price.append(s.price*1000) #$/MT
#         voc_cost.append(s.price*s.F_mass*get_annual_factor()/1e6) #million $/yr
#     if i in cooling_agents:
#         voc_baseline_sys.heat_utilities[5].cost #$/hr
#
# voc_data = pd.DataFrame({'ID':voc_ID_reformated,'Price [$/MT]':voc_price,'Cost [MM$/yr]':voc_cost})
