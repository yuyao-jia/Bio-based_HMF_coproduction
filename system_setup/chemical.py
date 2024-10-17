import thermosteam as tmo
from Unit_functions.Chemical_functions.cane_chemicals import (
    create_cellulosic_oilcane_chemicals,
    create_acyl_olein,)
import biosteam as bst
from chemicals import atoms_to_Hill
from thermosteam.utils import chemical_cache
from thermosteam import functional as fn

from system_setup.TAG_properties_yj import *
bst.nbtutorial()

# These defaults are used within system factories for pretreatment and fermentation.

#: Default liquid chemicals for saccharification solids-loading specification
default_nonsolids = ['Water', 'Ethanol', 'AceticAcid',
                     'Furfural', 'H2SO4', 'NH3', 'HMF']

#: Default insolible chemicals for saccharification solids-loading specification
default_insoluble_solids = ['Glucan', 'Mannan', 'Xylan',
                            'Arabinan', 'Galactan', 'Lignin']

#: Default ignored chemicals for saccharification solids-loading specification
default_ignored = ['TAG', 'DAG', 'MAG', 'FFA', 'PL']

#: Chemical groups for `get_grouped_chemicals` function
chemical_groups = dict(
        OtherSugars = ('Arabinose',
                       'Mannose',
                       'Galactose',
                       'Cellobiose',
                       'Sucrose'),
        SugarOligomers = ('GlucoseOligomer',
                          'XyloseOligomer',
                          'GalactoseOligomer',
                          'ArabinoseOligomer',
                          'MannoseOligomer'),
        OrganicSolubleSolids = ('AmmoniumAcetate',
                                'SolubleLignin',
                                'Extract',
                                'LacticAcid',
                                'Cellulase'),
        InorganicSolubleSolids = ('AmmoniumSulfate',
                                  'DAP',
                                  'NaOH',
                                  'HNO3',
                                  'NaNO3'),
        Furfurals = ('Furfural',
                     'HMF'),
        OtherOrganics = ('Glycerol',
                         'Denaturant',
                         'Oil',
                         'SuccinicAcid',
                         'Xylitol'),
        COxSOxNOxH2S = ('NO',
                        'NO2',
                        'SO2',
                        'CO',
                        'H2S'),
        Protein = ('Protein',
                   'Enzyme',
                   'DenaturedEnzyme'),
        CellMass = ('WWTsludge',
                    'Z_mobilis',
                    'T_reesei'),
        OtherInsolubleSolids = ('Tar',
                                'Ash',
                                'Lime'),
        OtherStructuralCarbohydrates = ('Arabinan',
                                        'Mannan',
                                        'Galactan')
)

chem = create_cellulosic_oilcane_chemicals(yeast_includes_nitrogen=None).copy()
#chem.show() can be used to check the chem list
# removed = {'HMF'}
# chem = tmo.Chemicals([
#         i for i in cellulosic_oilcane_chem.tuple if i.ID not in removed
#     ])
#adding chemicals existing in tmo.Chemical to my chemical list
#tmo.Chemical(search_ID="?") to search for compound in tmo.Chemical
FormicAcid=tmo.Chemical(ID="FormicAcid",search_ID="formic_acid")
chem.extend([FormicAcid]) #needs to be a list
PEG = tmo.Chemical(ID='PEG',search_ID='polyethylene_glycol')
chem.extend([PEG])
MTBE=tmo.Chemical(ID='MTBE',CAS='1634-04-4') #search_ID = Methyl tert-butyl ether
chem.extend([MTBE])
hexane=tmo.Chemical(ID='Hexane',search_ID='hexane')
chem.extend([hexane])
OleicAcid_veg = create_acyl_olein(0).copy('OleicAcid_veg')
MonoOlein_veg = create_acyl_olein(1).copy('MonoOlein_veg')
DiOlein_veg = create_acyl_olein(2).copy('DiOlein_veg')
TAG_veg = create_acyl_olein(3).copy('TriOlein_veg')
PL_veg = tmo.Chemical('Phosphatidylinositol', formula='C47H83O13P',
                                  search_db=False, CAS='383907-36-6', default=True,
                                  Hf=-1.779e6,  # Assume same as TAG on a weight basis
                                  aliases={'PL', 'PolarLipid'}, phase='l').copy('Phosphatidylinositol_veg')
chem.extend([OleicAcid_veg, MonoOlein_veg, DiOlein_veg, TAG_veg, PL_veg])
sodiumsulfate = tmo.Chemical(ID='SodiumSulfate', CAS='7757-82-6')
chem.extend([sodiumsulfate])
polyamide = tmo.Chemical(ID='Polyamide', CAS='63428-84-2')
chem.extend([polyamide])
for chemical in chem: chemical.default()#defaulting other missing properties if any to water
chem.compile()

soluble_solids = (chem.Xylose, chem.Arabinose, chem.XyloseOligomer,chem.ArabinoseOligomer,
                    chem.Cellobiose)
for chemical in soluble_solids:
    V = fn.rho_to_V(rho=1e5, MW=chemical.MW)
    chemical.V.add_model(V, top_priority=True)

chem.TriOlein.Tb = tmo.Chemical(ID='122-32-7').Tb
chem.TriOlein.Pc = tmo.Chemical(ID='122-32-7').Pc
chem.TriOlein.copy_models_from(tmo.Chemical(ID='122-32-7'),['Psat'])

chem.TriOlein_veg.Tb = tmo.Chemical(ID='122-32-7').Tb
chem.TriOlein_veg.Pc = tmo.Chemical(ID='122-32-7').Pc
chem.TriOlein_veg.copy_models_from(tmo.Chemical(ID='122-32-7'),['Psat'])

chem.set_synonym('OleicAcid_veg', 'FFA_veg')
chem.set_synonym('MonoOlein_veg', 'MAG_veg')
chem.set_synonym('DiOlein_veg', 'DAG_veg')
chem.set_synonym('TriOlein_veg', 'TAG_veg')
chem.set_synonym('Phosphatidylinositol_veg', 'PL_veg')

chem.define_group('Microbial_lipid', ['PL', 'FFA', 'MAG', 'DAG', 'TAG', 'AcTAG'])
chem.define_group('Oil', ['PL', 'FFA', 'MAG', 'DAG', 'TAG'])
chem.define_group('Vegetative_lipid', ['PL_veg','FFA_veg', 'MAG_veg', 'DAG_veg', 'TAG_veg'])
chem.define_group('Lipid',['PL', 'FFA', 'MAG', 'DAG', 'TAG', 'AcTAG','PL_veg','FFA_veg', 'MAG_veg', 'DAG_veg', 'TAG_veg'])
chem.define_group('Sugar',['Glucose','Sucrose'])
chem.define_group('Fiber',['Cellulose','Hemicellulose','Lignin'])
chem.define_group('Other',['Ash','Solids'])

#chem.show()

#System: move system to system.py later
# bst.settings.set_thermo(chem,cache=True)