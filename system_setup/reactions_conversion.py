# Modified from:
#Humbird, D., Davis, R., Tao, L., Kinchin, C., Hsu, D., and Aden, A. (2011). Process Design and Economics for Biochemical Conversion of Lignocellulosic Biomass to Ethanol: Dilute-Acid Pretreatment and Enzymatic Hydrolysis of Corn Stover. Renewable Energy, 147.
#Cortés-Peña, Y. R., Kurambhatti, C., Eilts, K., Singh, V., and Guest, J. S. (2022). Economic and Environmental Sustainability of Vegetative Oil Extraction Strategies at Integrated Oilcane and Oil-Sorghum Biorefineries. ACS Sustainable Chemistry & Engineering, 10(42), 13980–13990. https://doi.org/10.1021/acssuschemeng.2c04204
#using experimental data in:
#Jia, Y., Maitra, S., and Singh, V. (2023). Chemical-free production of multiple high-value bioproducts from metabolically engineered transgenic sugarcane ‘oilcane’ bagasse and their recovery using nanofiltration. Bioresource Technology, 371, 128630. https://doi.org/10.1016/j.biortech.2023.128630

import thermosteam as tmo
import biosteam as bst
from system_setup.chemical import chem
from tqdm import trange
#see massbalance.xlsx in model folder for detailed calculations of conversion

bst.preferences.update(flow='kg/hr', T='degC', P='bar', N=100, composition=False)
#sream volume = m3/hr, mass = kg/hr, mass/volume = kg/m3 = g/L
bst.preferences.save()
f= bst.Flowsheet('test')
f.clear()
bst.main_flowsheet.set_flowsheet(f)
# System: move system to system.py later
bst.settings.set_thermo(chem, cache=True)  # takes all chemicals and let biosteam know the system_setup we want to use
Rxn = tmo.reaction.Reaction
ParallelRxn = tmo.reaction.ParallelReaction
# pretreatment_rxns_ori are based on mol basis rxns from Humbird 2011 (used by Yoel)
#Assum HMF is only produced from glucan and furfural is only produced from xylan
pretreatment_rxn_ori=ParallelRxn([
        #            Reaction definition                 Reactant    Conversion
        Rxn('Glucan + H2O -> Glucose',                   'Glucan', 0.0990, chem),#udpated
        Rxn('Glucan + H2O -> GlucoseOligomer',           'Glucan', 0.0030, chem),#udpated
        Rxn('Glucan -> HMF + 2 H2O',                     'Glucan', 0.0030, chem),#udpated
        Rxn('Galactan + H2O -> GalactoseOligomer',       'Galactan', 0.0240, chem),#same as Yoel
        # Rxn('Galactan -> HMF + 2 H2O',                   'Galactan', 0.0030, chem),
        Rxn('Xylan + H2O -> Xylose',                     'Xylan', 0.9000, chem),#udpated
        Rxn('Xylan + H2O -> XyloseOligomer',             'Xylan', 0.0240, chem),#same as Yoel
        Rxn('Xylan -> Furfural + 2 H2O',                 'Xylan', 0.0500, chem),#udpated
        Rxn('Arabinan + H2O -> Arabinose',               'Arabinan', 0.9000, chem),#same as Yoel
        Rxn('Arabinan + H2O -> ArabinoseOligomer',       'Arabinan', 0.0240, chem),#same as Yoel
        # Rxn('Arabinan -> Furfural + 2 H2O',              'Arabinan', 0.0050, chem),
        Rxn('Acetate -> AceticAcid',                     'Acetate', 1.0000, chem),#same as Yoel
        Rxn('Lignin -> SolubleLignin',                   'Lignin', 0.0500, chem),#udpated
        Rxn('Sucrose -> HMF + Glucose + 2H2O',           'Sucrose', 1.0000, chem),#same as Yoel
])

#convert mol basis rxns to wt basis
rxn_wt = []
for i in trange(12):
    rxn_wt.append(pretreatment_rxn_ori[i].copy(basis='wt'))
pretreatment_rxn_wt = ParallelRxn(rxn_wt)

glucan_to_glucose = pretreatment_rxn_wt[0]
glucan_to_glucoseoligomer = pretreatment_rxn_wt[1]
glucan_to_hmf = pretreatment_rxn_wt[2]
xylan_to_xylose = pretreatment_rxn_wt[4]
xylan_to_furfural = pretreatment_rxn_wt[6]
aceticacid_production = pretreatment_rxn_wt[9]
lignin_solubilization = pretreatment_rxn_wt[10]

#update conversion (X%) of reactions based on experimental data (Jia et al., 2023)
glucan_to_glucose.X = 0.2145
glucan_to_glucoseoligomer.X = 0.0794
glucan_to_hmf.X = 0.1440
xylan_to_xylose.X = 0.3123
xylan_to_furfural.X = 0.5110
lignin_solubilization.X = 0.1097

#Enzyamtic hydrolysis reactions
hydrolysis_rxn_ori = ParallelRxn([
    Rxn('Glucan -> GlucoseOligomer', 'Glucan', 0.0400, chem),
    Rxn('Glucan + 0.5 H2O -> 0.5 Cellobiose', 'Glucan', 0.0120, chem),
    Rxn('Glucan + H2O -> Glucose', 'Glucan', 0.9000, chem),
    Rxn('Cellobiose + H2O -> 2Glucose', 'Cellobiose', 1.0000, chem)]
)
#convert mol basis rxns to wt basis
hrxn_wt = []
for i in trange(4):
    hrxn_wt.append(hydrolysis_rxn_ori[i].copy(basis='wt'))
hydrolysis_rxn_wt = ParallelRxn(hrxn_wt)
glucan_to_glucose_hydrolysis = hydrolysis_rxn_wt[2]
glucan_to_glucose_hydrolysis.X = 0.9347 #based on experimental data (Jia et al., 2023),X is % of reactant converted (g converted reactant/g total reactant)
#Also works: glucan_to_glucose_hydrolysis.product_yield('Glucose',basis='wt',product_yield=1.0375),product_yield is g product/g reactant
