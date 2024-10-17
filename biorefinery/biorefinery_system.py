import biosteam as bst
from system_setup.chemical import chem
from biosteam import main_flowsheet as F
from biorefinery.create_cellulosic_biodiesel_hmf import (
    create_oilcane_to_biodiesel_individual_oil_separation_microbial_vegetative)
from system_setup.prices_yj import prices_per_Kg as price
from system_setup import process_settings_yj
from system_setup.lca_characterization_factors_yj import GWP_characterization_factors as GWP_factors
from biosteam import preferences

bst.nbtutorial()
preferences.update(flow='kg/hr', T='degC', P='Pa', N=100, composition=True)
#sream volume = m3/hr, mass = kg/hr, mass/volume = kg/m3 = g/L
preferences.save()

F.set_flowsheet('oilcane_baseline')
bst.settings.set_thermo(chem, cache=True)  # takes all chemicals and let biosteam know the system_setup we want to use
process_settings_yj.load_process_settings()
GWP = 'GWP100'
bst.settings.set_electricity_CF(GWP, GWP_factors['electricity'], basis='kWhr', units='kg*CO2e')
bst.settings.electricity_price = price['Electricity']  # 2023 average price for industrial use (U.S. Energy Information Administration (EIA), 2023)

oilcane = bst.Stream('oilcane',
                          Ash=2000.042,
                          Cellulose=26986.69,
                          Glucose=2007.067,
                          Hemicellulose=15922.734,
                          Lignin=14459.241,
                          TAG_veg=8028.267200000001,
                          PL_veg=1003.5334000000001,
                          FFA_veg=1003.5334000000001,
                          Solids=5017.667,
                          Sucrose=22746.761,
                          Water=234157.798,
                          units='kg/hr',
                          price=0.03455,
                          )


oilcane_baseline_sys = create_oilcane_to_biodiesel_individual_oil_separation_microbial_vegetative(ins=oilcane,fed_batch=False,urea=False)
sys = F.create_system('baseline_sys')

# # #adjust bagasse to boiler to satisfy energy demand
# splitter = F.unit.S201
#
# minimum_fraction_burned = 0
# maximum_fraction_burned = 0.3
# natural_gas_streams = [F.stream['natural_gas']]
# @sys.add_bounded_numerical_specification(
#     x0=minimum_fraction_burned, x1=maximum_fraction_burned,
#     xtol=1e-4, ytol=100, args=(splitter,),
# )
# def adjust_bagasse_to_boiler(fraction_burned, splitter):
#     # Returns energy consumption at given fraction processed (not sent to boiler).
#     splitter.split[:] = 1 - fraction_burned
#     sys.simulate()
#     BT = F.unit.BT701
#     excess = BT._excess_electricity_without_natural_gas
#     if fraction_burned == minimum_fraction_burned and excess > 0:
#         splitter.neglect_natural_gas_streams = False  # No need to neglect
#         return 0  # No need to burn bagasse
#     elif fraction_burned == maximum_fraction_burned and excess < 0:
#         splitter.neglect_natural_gas_streams = False  # Cannot be neglected
#         return 0  # Cannot satisfy energy demand even at 30% sent to boiler (or minimum fraction processed)
#     else:
#         splitter.neglect_natural_gas_streams = True
#         return excess
# @sys.add_specification(args=(splitter,))
# def assume_negligible_natural_gas_streams(splitter):
#         if splitter.neglect_natural_gas_streams:
#             for i in natural_gas_streams: i.empty()

