# -*- coding: utf-8 -*-
"""
Created on Thu Nov  4 14:28:33 2021

@author: yrc2
"""

__all__ = (
    'GWP_characterization_factors',
    'set_GWPCF',
    'GWP',
)
lbsperMWh_to_KgsperKWh = 4.53592e-4
GWP = 'GWP100'
#try to find GWP factor for US, not the rest of the world (RoW)
#time doesn't really matter for GWP factor
# All values in cradle-to-gate except for CH4, which is in cradle-to-grave
GWP_characterization_factors = { # Material GWP cradle-to-gate [kg*CO2*eq / kg]
    # 'PEG': 3.7786, #RoW,estimated by modeling PEG using polyol library found in SimaPro software (Riofrio et al., 2022)
    'MTBE': 1.14, #Ecoinvent, market for methyl tert-butyl ether, GLO
    'sugarcane': 0.02931 * 0.30 / 0.25, # GREET, modified from moisture content of 0.75 to 0.70
    # 'sweet sorghum': 0.02821 * 0.30 / 0.25, # GREET, modified from moisture content of 0.75 to 0.70
    #'feedstock': 0.0607, # dry basis, Ecoinvent 2021,
    # 'protease': 8.07, # Assume same as cellulase
    'cellulase': 8.07, # GREET
    'H3PO4': 1.00, # GREET
    'lime': 1.164, # GREET
    # 'urea': 1.81,
    'hexane': 0.554, # Ecoinvent, market for hexane, GLO
    'pure-glycerol': 1.6678, # Ecoinvent, TRACI, market for glycerine, RoW; 
    'crude-glycerol': 0.36, # GREET
    # 'biodiesel': 1.13, # Soybean biodiesel GREET
    'HCl': 1.96, # GREET
    'NaOH': 2.01, # GREET
    'methanol': 0.45, # GREET, Natural gas to methanol
    'NaOCH3': 1.5871, # Ecoinvent, TRACI, sodium methoxide
    'pwc_makeup_water':0.00035559, #Ecoinvent
    'CH4': 0.33, # Natural gas from shell conventional recovery, GREET; includes non-biogenic emissions
    'electricity': 852.3*lbsperMWh_to_KgsperKWh, #Based on NAtional US average for electricty CO2 emissions [24]
    # 0.66 is the GWP from producing diesel from GREET; Conventional diesel from crude oil for US Refineries.
    # Downstream fuel emissions are added in. Accounts for how biodiesel has less energy than diesel.
    'biodiesel displacement': 0.92 * (0.66 +  (12 * 12.01 + 24 * 16) / (12 * 12.01 + 23 * 1.008)) ,
    'glycerine displacement': 3.73, #EcoInvent, 3.73 is GWP from producing glycerine from epicholohydrin
    'HMF displacement': 1.7, #EcoInvent, 1.7 is GWP from producing xylene from crude petroleum
    'furfural displacement':1.61,
    #0.831, #EcoInvent, 0.831 is GWP for formaldehyde
    #oil-based tolune: 1.61 #ecoinvent
    'electricity displacement': 1.86, #Ecoinvent, 1.86 is GWP for electricity from oil in US
    # Ecoinvent, market for sodium hypochlorite, without water, in 15% solution state, RoW,
    # converted to 12.5 wt% solution (15 vol%)
    'NaOCl': 2.4871 * 0.125,
    'Cooling_tower_chemicals': 0.2574,  # Based on sodium bicarbonate which is a common corrosion inhibitor [20].
    # Data for other cooling tower chemicals such as pH controllers and anti-scalants not found.
    'Boiler_chemicals': 1.5568,
    # Based on the production of sulfite which is commonly used as an oxygen scavenger in the boiler [20]
    # Data for other boiler chemicals (neutralising amines) not found.
    #High rate WWT only:
    'Citric_acid': 4.37, # Ecoinvent
    'Bisulfite':0.44,
    'NF membrane': 0.686,  # per m2 of membrane, Sen ́ an-Salinas et al. 2022
    'NF membrane maintenance': 0.02565,  # per m2 of membrane, Bonton et al. 2012
}

GWP_characterization_factors['methanol catalyst mixture'] = (
    GWP_characterization_factors['methanol'] * 0.75 + GWP_characterization_factors['NaOCH3'] * 0.25
)

def set_GWPCF(stream, name, dilution=1.):
    stream.characterization_factors[GWP] = GWP_characterization_factors[name] * dilution

# from thermosteam.units_of_measure import convert
# from thermosteam import Chemical
# CH4 = Chemical('CH4')
# CO2 = Chemical('CO2')
# electricty_produced_per_kg_CH4 = - convert(0.8 * 0.85 * CH4.LHV / CH4.MW, 'kJ', 'kWhr')
# GWP_per_kg_CH4 = 0.33 + CO2.MW / CH4.MW
# GWP_per_kWhr = GWP_per_kg_CH4 / electricty_produced_per_kg_CH4