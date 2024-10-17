#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd

biomass_based_diesel_price_2023 = pd.read_excel('biodiesel_price_2023.xlsx')['Biodiesel per kg']
biodiesel_avg_2023 = np.average(biomass_based_diesel_price_2023.values)

lb_to_kg = 0.453
cooling_tower_chem_density = 8.65  # lbs/gal [16]
boiler_chem_density = 9.10  # lb/ga[17]
cents_to_dollar = 100
tcf_to_cf = 1000  # thousand cubic feet to cubic feet
cf_to_m3 = 0.029  # cubic feet to cubic meter
kgal_to_L = 3.785412 * 1000

# Index was chosen for the general category of all chemicals and allied products
PPI_2016 = 270.400  # Dec 2016
PPI_2021 = 349.906  # Dec 2021
PPI_2022 = 357.404  # Dec 2022
PPI_2023 = 347.123  # Sept 2023
PPI_1997 = 143.500  # Dec 1997
PPI_2018 = 293.000  # Dec 2018
PPI_utility_2016 = 136.200  # Dec 2016
PPI_utility_2017 = 139.60  # Dec 2017
PPI_utility_2018 = 151.10  # Dec 2018
PPI_utility_2019 = 141.10  # Dec 2019
PPI_utility_2020 = 146.50  # Dec 2020
PPI_utility_2021 = 171.89  # Dec 2021
PPI_utility_2022 = 222.043  # Dec 2022
PPI_utility_2017_2022 = [PPI_utility_2017,
                         PPI_utility_2018,
                         PPI_utility_2019,
                         PPI_utility_2020,
                         PPI_utility_2021,
                         PPI_utility_2022]
PPI_utility_2023 = 197.546  # Sept 2023

ratio_2023from2021 = PPI_2023 / PPI_2021
ratio_2023from2022 = PPI_2023 / PPI_2022
ratio_2023from1997 = PPI_2023 / PPI_1997
ratio_2023from2016 = PPI_2023 / PPI_2016
ratio_2023from2018 = PPI_2023 / PPI_2018
ratio_utility_2023from2016 = PPI_utility_2023 / PPI_utility_2016
ratio_utility_2023from2017_2022_period_mean = PPI_utility_2023 / np.mean(PPI_utility_2017_2022)

natural_gas_prices_2017_to_2022 = [4.08, 4.19, 3.90, 3.32, 5.44,
                                   7.66]  # Provided by EIA, Dollars per Thousand Cubic Feet
average_electricity_prices_2017_to_2022 = [6.88, 6.92, 6.81, 6.67, 7.18, 8.45]  # Provided by EIA, cents/Kwh
average_water_rates_2008_2012_2016 = [2.44, 3.02, 3.38]  # $/kGal Provided by USDOE


def convert_lab_price_to_bulk_price(
        lab_quantity_in_kg=np.array([]),
        lab_price_per_kg=np.array([]),
        bulk_quantity_in_kg=30,  # used by (Hart and Sommerfeld, 1997)
        bulk_coeff=-0.75  # as mentioned in (Hart and Sommerfeld, 1997)
):
    bulk_price = np.array(
        [(lab_price_per_kg[i] * (bulk_quantity_in_kg / lab_quantity_in_kg[i]) ** (bulk_coeff)) for i in
         range(len(lab_quantity_in_kg))])
    return np.average(bulk_price)


HMF_bulk_price = convert_lab_price_to_bulk_price(lab_quantity_in_kg=np.array([
    1,
    1,
]),
    lab_price_per_kg=np.array([
        210,  # https://www.ambeed.com/products/67-47-0.html
        252,  # https://www.chemenu.com/products/67-47-0.html
    ]),
    bulk_quantity_in_kg=30, bulk_coeff=-0.75),

prices_per_Kg = {  # 2023 prices in USA dollars
    # Raw material prices:
    'Oilcane': 0.03455,
    # temporarily assumed to be the same price as sugarcane, maximum feedstock purchase price (MFPP) the biorefinery can afford will be computed.
    'H3PO4': 1.28 / 160.409 * 122.236,  # FRED converted from 2021 price on CatCost
    'Lime': 0.12 * 415.715 / 274.695,  # 2023 price, FRED & CatCost
    'HCl': 0.0718 * 112.111 / 95.300,  # 2023 price, FRED converted from Yoel 2019 price
    'Polymer': 0.0,
    'Urea': 0.415 / 88.623 * 183.565,  # 2023 price, FRED converted from Yoel 2019 price
    'Hexane': 1170 / 1000,  # 2023 price, https://www.chemanalyst.com/Pricing-data/n-hexane-1151
    'Acetone':0.98,#2023 price, https://businessanalytiq.com/procurementanalytics/index/acetone-price-index/
    'Pure_glycerine': 1.09,
    # 2023 price, https://businessanalytiq.com/procurementanalytics/index/glycerol-price-index/
    # 'PWC_makeup_water':0.0010,#makeup water only, 2023 price
    'PWC_makeup_water': np.mean(average_water_rates_2008_2012_2016) * ratio_utility_2023from2016 / kgal_to_L,
    # [13],[14]
    'Natural_gas_price': np.mean(natural_gas_prices_2017_to_2022) * ratio_utility_2023from2017_2022_period_mean / (
                tcf_to_cf * cf_to_m3),
    # Industrial natural gas price provided by EIA for 2017-2022 period[15]
    'Cellulase': 0.212 / 285.100 * 347.123,  # Fred converted from Yoel 2019 price
    'NaOH': 0.39,
    # 2023 price, (Mike, 2020) &https://businessanalytiq.com/procurementanalytics/index/sodium-hydroxide-price-index/
    'Methanol': 0.536,  # 2023 price, (Methanol Price Trend and Forecast, 2023)
    'Catalyst(NaOCH3)': 2.93 / 300.500 * 446.215,  # 2023 price, FRED converted from Yoel 2019 price
    # 'Cooling_tower_chemicals':3.933/278.300*337.082,#2023 price, FRED converted from Yoel 2019 price (inorganic chemicals other than alkalies and chlorine)
    'Cooling_tower_chemicals': 519 / (55 * cooling_tower_chem_density * lb_to_kg),
    # 2023 price,Based on [16], 55 gal drum costs 519$
    # 'Boiler_chemicals':6.563/175.000*224.159, #2023 price, FRED converted from Yoel 2019 price (water-treating compounds)
    'Boiler_chemicals': 1.216 / (55 * boiler_chem_density * lb_to_kg),
    # 2023 price,Based on [17],55 gal drum costs 609$
    'Natural_gas': 3.90 / (1000 * 0.68 * 0.0283),
    # assume natural gas density is 0.68 kg/m3, 2023 price, https://www.eia.gov/dnav/ng/ng_pri_sum_dcu_nus_m.htm
    'PEG': 5500 / 1000 / 111.785 * 100.038,
    # PEG600, 2023 price FRED (other synthetic organic chemicals) converted from its 2022 Sep price https://www.chemanalyst.com/Pricing-data/polyethylene-glycol-peg-1171
    'MTBE': 1.15,  # 2023 price, https://www.chemanalyst.com/Pricing-data/methyl-tertiary-butyl-ether-81
    'N2': 755 / 1000 / 285.992 * 295.224,
    # 2023 price, FRED converted from 2022 Sep price obtained from https://www.chemanalyst.com/Pricing-data/bulk-nitrogen-1097
    'Electricity': np.mean(
        average_electricity_prices_2017_to_2022) * ratio_utility_2023from2017_2022_period_mean / cents_to_dollar,
    # Industrial retail electricity price provided by EIA for 2017-2022 period[18]
    'Propane': 1.027458333 / 3.785 / 0.493,

    # High rate WWT only:
    'Citric_acid': 896 / 1000,
    # 2023 price, https://www.chemanalyst.com/Pricing-data/citric-acid-1438#:~:text=Between%20July%20and%20September%202023%2C%20the%20CFR%20New,levels%2C%20with%20an%20average%20quarterly%20increase%20of%200.17%25.
    'Bisulfite': 0.82 * 436.327 / 348.496,  # FRED converted from CatCost FOB price
    # By products prices:
    'Crude_glycerol': 0.20 / 243.100 * 308.325,  # 2023 price, FRED converted from Yoel 2019 price
    'Ash_disposal': -0.0318,  # default biosteam price
    # Product prices:
    # 'Biodiesel':biodiesel_avg_no_RIN,#B100 2023 July price (Bourbon & Science, 2023)
    # 'Cellulosic based diesel': cellulosic_diesel_avg,
    'Biomass based diesel': biodiesel_avg_2023,
    'Purified_HMF': HMF_bulk_price[0],
    'Purified_furfural': 1.501 / 115.961 * 146.661,
    # 2023 price, FRED converted from 2021 price (Grand View Research, 2022)
}

prices_per_Kg['methanol catalyst mixture'] = (
        prices_per_Kg['Methanol'] * 0.75 + prices_per_Kg['Catalyst(NaOCH3)'] * 0.25
)


def set_price(stream, name, dilution=1.):
    stream.price = prices_per_Kg[name] * dilution

# Uncertainity analysis
