# Description: Streams for the biorefinery system
from biorefineries.biodiesel.streams import *
from biorefineries.cane.streams import *
from biorefineries.cellulosic.streams import *

bagasse = stream_kwargs('bagasse')
bagasse_pellets = stream_kwargs('bagasse_pellets')
sugarcane = stream_kwargs('sugarcane',
    Water=0.7,
    Glucose=0.01208,
    Sucrose=0.1369,
    Ash=0.006,
    Cellulose=0.06115,
    Hemicellulose=0.03608,
    Lignin=0.03276,
    Solids=0.015,
    total_flow=333334.2,
    units='kg/hr',
    price=0.03455
)
lipidcane = stream_kwargs('lipidcane',
    Ash=2000.042,
    Cellulose=26986.69,
    Glucose=2007.067,
    Hemicellulose=15922.734,
    Lignin=14459.241,
    TAG=10035.334,
    Solids=5017.667,
    Sucrose=22746.761,
    Water=234157.798,
    units='kg/hr',
    price=0.03455,
)
oilcane = lipidcane.copy()
oilcane['ID'] = 'oilcane'
TAG = oilcane['TAG']
oilcane['TAG'] = 0.80 * TAG
oilcane['PL'] = 0.10 * TAG
oilcane['FFA'] = 0.10 * TAG
cane = stream_kwargs('cane')
condensate = stream_kwargs('condensate')
bagasse_to_boiler = stream_kwargs('bagasse_to_boiler')
shredded_cane = stream_kwargs('shredded_cane')
untreated_juice = stream_kwargs('untreated_juice')
H3PO4 = stream_kwargs('H3PO4',
    H3PO4=74.23,
    Water=13.1,
    units='kg/hr',
    price=0
)
lime = stream_kwargs('lime',
    CaO=333.0,
    Water=2200.0,
    units='kg/hr',
    price=0.077
)
polymer = stream_kwargs('polymer',
    Flocculant=0.83,
    units='kg/hr',
    price=0
)
clarified_juice = stream_kwargs('clarified_juice')
screened_juice = stream_kwargs('screened_juice',
    Glucose=3802,
    Sucrose=4.309e+04,
    Water=2.59e+05,
    H3PO4=83.33,
    units='kg/hr',
    T=372
)
stillage = stream_kwargs('stillage')
fiber_fines = stream_kwargs('fiber_fines')
evaporator_condensate = stream_kwargs('evaporator_condensate')
vent = stream_kwargs('vent')
vinasse = stream_kwargs('vinasse')
wastewater = stream_kwargs('wastewater')
emissions = stream_kwargs('emissions')
ash_disposal = stream_kwargs('ash_disposal')
molasses = stream_kwargs('molasses')
lipid = stream_kwargs('lipid')
cellmass = stream_kwargs('cellmass')
fermentation_effluent = stream_kwargs('fermentation_effluent')
spent_oil_wash_water = stream_kwargs('spent_oil_wash_water')
sugar = stream_kwargs('sugar', price=0.419) # https://markets.businessinsider.com/commodities/sugar-price?op=1 (3/18/2022)
acTAG = stream_kwargs('acTAG', price=1.633)
crude_oil = stream_kwargs('crude_oil', price=0.661) # 30 cts / lb (vegetable); 5yr average https://www.ams.usda.gov/mnreports/lswagenergy.pdf

biodiesel = stream_kwargs('biodiesel', price=1.38) # assuming density = 870 kg/m3
crude_glycerol = stream_kwargs('crude_glycerol',price=0.20)
wastewater = stream_kwargs('wastewater', price=-0.33)
acetone = stream_kwargs('acetone', price=0.80)
pure_glycerine = stream_kwargs('pure_glycerine', price=0.65)
degummed_oil = stream_kwargs('degummed_oil')
polar_lipids = stream_kwargs('polar_lipids')
wastewater = stream_kwargs('wastewater')
lipid_wash_water = stream_kwargs('lipid_wash_water')
washed_lipid = stream_kwargs('washed_lipid')
spent_wash_water = stream_kwargs('spent_wash_water')

#added by YJ
pretreatment_liquor = stream_kwargs('pretreatment_liquor')
pretreated_biomass = stream_kwargs('pretreated_biomass')
juice = stream_kwargs('juice')
hydrolysate = stream_kwargs('hydrolysate')
slurry = stream_kwargs('slurry')

purified_hmf = stream_kwargs('purified_hmf',price = 1.03)
purified_furfural = stream_kwargs('purified_furfural',price=1.8203)
concentrated_pretreatment_liqour = stream_kwargs('concentrated_pretreatment_liqour')
purification_wastewater = stream_kwargs('purification_wastewater')

oil_from_lignin = stream_kwargs('oil_from_lignin')
spent_lignin_to_boiler= stream_kwargs('spent_lignin_to_boiler')
extracted_cellmass_to_boiler = stream_kwargs('extracted_cellmass_to_boiler')
saccharified_biomass_solids = stream_kwargs('saccharified_biomass_solids')
water_condensate = stream_kwargs('water_condensate')
filter_cake = stream_kwargs('filter_cake')
nanofiltration_membrane = stream_kwargs('nanofiltration_membrane')