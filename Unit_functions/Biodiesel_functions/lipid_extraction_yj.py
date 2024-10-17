# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import flexsolve as flx
import numpy as np
import biosteam as bst
from biosteam import SystemFactory
from system_setup import stream as s
import thermosteam as tmo
from system_setup import units_added as u

from biorefineries.biodiesel import (
    create_lipid_wash_system,
)
__all__ = (
    'create_lipid_exctraction_system',
    'create_post_fermentation_oil_separation_system',
    'create_lignin_oil_extraction_system',
)

@SystemFactory(
    ID='lipid_exctraction_sys',
    ins=[s.fermentation_effluent],
    outs=[s.lipid, s.cellmass, s.wastewater],
)
def create_lipid_exctraction_system(ins, outs):
    fermentation_effluent, = ins
    lipid, cellmass, wastewater, = outs
    U401 = bst.SolidsCentrifuge('U401', fermentation_effluent, ['', ''],
        split=dict(cellmass=0.98, lipid=0.98),
        moisture_content=0.5,
        solids=['cellmass'],
    )
    U402 = bst.DrumDryer('U402', 
        (U401-0, 'dryer_air', 'dryer_natural_gas'), 
        ('', 'dryer_outlet_air', 'dryer_emissions'),
        moisture_content=0.18, split=0.,
        utility_agent='Steam',
    )
    # X401 = bst.ThermalOxidizer('X401', (U403-1, 'oxidizer_air'), 'oxidizer_emissions')
    U403 = bst.ScrewPress('U403', U402-0, split=dict(cellmass=1, lipid=0.3, Water=0.8),)
    bst.ConveyingBelt('U405', U403-0, cellmass)
    lipid_wash_sys = create_lipid_wash_system(ins=U403-1, outs=lipid, mockup=True)
    washed_lipid, spent_wash_water = lipid_wash_sys.outs
    bst.Mixer(ins=[spent_wash_water, U401-1], outs=wastewater)

@SystemFactory(
    ID='post_fermentation_oil_separation_sys',
    ins=[s.stillage],
    outs=[s.lipid, s.extracted_cellmass_to_boiler,s.wastewater],
)
def create_post_fermentation_oil_separation_system(ins, outs,
                                                   # wastewater_concentration=None,
                                                   # target_oil_and_solids_content=60,
                                                   ):
                                                   # separate_cellmass=False):
    lipid, extracted_cellmass_to_boiler,wastewater = outs
    # if separate_cellmass:
    #     cellmass = bst.Stream('cellmass')
    #     outs.insert(1, cellmass)
    V605 = bst.MixTank('V605', ins)
    # if evaporation:
    #     evaporator_condensate = bst.Stream('evaporator_condensate')
    #     outs.insert(3, evaporator_condensate)
    #     P606 = bst.Pump('P606', V605 - 0)
    #     EvX = bst.MultiEffectEvaporator('Ev607',
    #         ins=P606-0,
    #         P=(101325, 69682, 47057, 30953),
    #         V=0.53, V_definition='Overall', #V is set to reach 6% of veg TAG in the total liquid, assuming
    #         thermo=lipid.thermo.ideal(),
    #         flash=False,
    #     ) #change V to adjust solid and lipid content
    #     EvX.target_oil_and_solids_content = target_oil_and_solids_content # kg / m3
    #     EvX.remove_evaporators = False
    #     P_original = tuple(EvX.P)
    #     Pstart = P_original[0]
    #     Plast = P_original[-1]
    #     N = len(P_original)
    #     P607 = bst.Pump('P607', EvX - 0, P=101325.)
    # # def x_oil(V): # Objective function for specification
    # #     EvX.V = V
    # #     EvX.run()
    # #     effluent = EvX.outs[0]
    # #     moisture = effluent.imass['Water']
    # #     total = effluent.F_mass
    # #     if total == 0:
    # #         return 0.
    # #     else:
    # #         return EvX.target_oil_and_solids_content - 1000. * (1. - moisture / total)
    # #
    # # @EvX.add_specification(run=False)
    # # def adjust_evaporation():
    # #     V_last = EvX.V
    # #     x0 = 0.
    # #     x1 = 0.5
    # #     EvX.P = P_original
    # #     EvX._reload_components = True
    # #     y0 = x_oil(x0)
    # #     if y0 <= 0.:
    # #         EvX.V = x0
    # #         return
    # #     else:
    # #         EvX._load_components()
    # #         for i in range(1, N):
    # #             if x_oil(1e-6) < 0.:
    # #                 EvX.P = np.linspace(Pstart, Plast, N - 1)
    # #                 EvX._reload_components = True
    # #             else:
    # #                 break
    # #         y1 = x_oil(x1)
    # #         EvX.V = flx.IQ_interpolation(x_oil, x0, x1, y0, y1, x=V_last, ytol=1e-5, xtol=1e-6)
    # else:
    P606 = bst.Pump('P606', V605 - 0)
    # C603_2 = bst.LiquidsSplitCentrifuge('C603_2', P607-0, (lipid, ''),
    #                                     split={'Lipid': 0.99,
    #                                            'Water': 0.0001})
    # C603 = bst.SolidLiquidsSplitCentrifuge('decanter',
    #                                            ins=P606-0,
    #                                            # outs=[lipid, '', ''],
    #                                            # outs=[backend_oil, aqueous_stream, cellmass],
    #                                            solids_split={'Cellmass':1.0, 'Lipid':0.99},
    #                                            #95% of the oil is in within cell, 1% is bagasse oil within broth
    #                                            aqueous_split={'Lipid':0.01,'Water':0.05},
    #                                            # oil_centrifuge.split = oil.split = 0.1, since most oil are within cells at this stage
    #                                            moisture_content=0.4,
    #                                            )
    C601=bst.SolidsCentrifuge('C601', P606-0,outs=['', wastewater],
                              split=dict(Microbial_lipid=0.99, Cellmass=1,Vegetative_lipid = 0.0), #vegetative lipid that remains in liquid phase is removed as waste
                              solids=('Cellmass','Microbial_lipid'),moisture_content=0.4)
    # U601 = bst.DrumDryer('U601',
    #     ins=(C601-0, 'dryer_air', 'dryer_natural_gas'),
    #     outs=('', 'dryer_outlet_air', 'dryer_emissions'),
    #     moisture_content=0.18, split=0.,
    #     utility_agent='Steam',
    #     #no natural gas is used since utility agent is steam
    #     #no emission since no natural gas is used
    # )
    # U602 = bst.ScrewPress('U602', ins=U601-0,outs=[lipid, ''],
    #                                    split=dict(Solids=0, Ash=0, Lignin=0, Glucan=0, Xylan=0,
    #                                               Arabinan=0, Galactan=0, Lipid=0.8, Water=0.2),)
    microbial_lipid_extraction_sys= create_lignin_oil_extraction_system(ins=C601-0,outs=(lipid,''),extractor_type='percolation',name=True,lipid_type = 'microbial')
     #20% lipid are lost in the cellmass
    #drying extracted cellmass to 40% moisture content before sending to boiler
    U601 = bst.DrumDryer('U601',ins=(microbial_lipid_extraction_sys.outs[1], 'dryer_air', ''),
                        outs=('', 'dryer_outlet_air', ''),
                        moisture_content=0.40, split=0.,utility_agent='Steam',)
    bst.ConveyingBelt('U602', U601 - 0, outs=extracted_cellmass_to_boiler)
    # if evaporation:
    #     S601 = bst.Splitter('S601', ins=EvX-1, outs=['', evaporator_condensate], split=0.5)
    #     M601 = bst.Mixer('M601', [S601-0, C601-0], wastewater)
    # M601.target_wastewater_concentration = 60. # kg / m3
    # @M601.add_specification(run=True, impacted_units=[S601])
    # def adjust_wastewater_concentration():
    #     concentrated_wastewater = M601.ins[1]
    #     waste = concentrated_wastewater.F_mass - concentrated_wastewater.imass['Water']
    #     current_concentration = waste / concentrated_wastewater.F_vol
    #     required_water = (1./M601.target_wastewater_concentration - 1./current_concentration) * waste * 1000.
    #     F_mass = S601.ins[0].F_mass
    #     if F_mass:
    #         split = required_water / F_mass
    #         if split < 0:
    #             split = 0.
    #         elif split > 1.:
    #             split = 1.
    #         S601.split[:] = split


@SystemFactory(
    ID='lignin_oil_extraction_sys',
    ins=[s.saccharified_biomass_solids],
    outs=[s.oil_from_lignin,s.spent_lignin_to_boiler]
)
def create_lignin_oil_extraction_system(ins, outs,lignin_extraction_area=None,extractor_type='percolation',name=None,lipid_type = 'microbial'):
    saccharified_biomass_solids = ins
    oil_from_lignin,spent_lignin_to_boiler = outs
    if lignin_extraction_area is None:lignin_extraction_area=500
    # Hexane oil extraction
    if name is None:
        hexane = bst.Stream('hexane1', Hexane=28834.83, units='kg/hr')
        hexane_recycled = bst.Stream('Recycled_Hexane1', Hexane=0.1, units='kg/hr')
    else:
        hexane = bst.Stream('hexane2', Hexane=23090.58, units='kg/hr')
        hexane_recycled = bst.Stream('Recycled_Hexane2', Hexane=0.1, units='kg/hr')

    mixer = bst.Mixer(lignin_extraction_area, ins=(hexane,hexane_recycled))
    T401 = bst.StorageTank(lignin_extraction_area, mixer - 0)
    P401 = bst.units.Pump(lignin_extraction_area, ins=T401 - 0)  # pump hexane to extractor

    partition_data_veg = {
        'K': np.array([0.0001,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2500,2500,5, 53306.45]),
        'IDs': ('Water','Glucose', 'Ash','Lignin','Solids','AceticAcid','Furfural','HMF','Xylose','SolubleLignin','GlucoseOligomer','GalactoseOligomer',
                'XyloseOligomer','ArabinoseOligomer','Glucan','Xylan','Cellobiose','Arabinan','Galactan','Cellulase','OleicAcid_veg','TriOlein_veg',
                'Phosphatidylinositol_veg','Hexane'),
    }
    partition_data_microbial = {'K': np.array([0.0001,0,200,53306.45]),
                                'IDs':('Water', 'Yeast','TriOlein','Hexane'),
                                }
    if extractor_type == 'immersion':
        U401 = bst.units.MixTank(lignin_extraction_area, ins=(P401 - 0, saccharified_biomass_solids),tau=12)  # extractor at room temp

    else:  # MultiStageMixerSettlers: counter-current extraction to simulate percolation (Cheng et al., 2017),separate solids after extraction by settler
        U401 = bst.MultiStageMixerSettlers(lignin_extraction_area, ins=(ins[0], P401 - 0), outs=(spent_lignin_to_boiler,''),N_stages=1,solvent_ID='Hexane')
    def calculate_lipid_extraction_efficiency():
        return U401.outs[1].imass['Lipid']/(U401.ins[0].imass['Lipid']+U401.ins[1].imass['Lipid'])  # temporary random guess
    U401.lipid_extraction_efficiency = calculate_lipid_extraction_efficiency # final lipid extraction efficiency
    if lipid_type == 'vegetative':
        U401.partition_data = partition_data_veg
    elif lipid_type == 'microbial':
        U401.partition_data = partition_data_microbial

    # @U401.add_specification(run=True)
    # def change_N_stages():
    #     stage_guess = 1
    #     stage_desired = flx.IQ_interpolation(get_lipid_extraction_efficiency(),x0=0.01,x1=0.99,x=stage_guess,ytol=0.01,maxiter=1000)

        #hexane typically extracts neutral lipids. PL are not typically extracted by hexane.
        #Ask shraddha: should I skip acetone washing?
    # #Composition of extracted cake:75% meal, 25% hexane (Le Clef & Kemper, 2015)
    # #Composition of miscella: 23% oil, 77% hexane (Le Clef & Kemper, 2015)
    # miscella = U401 - 1
    # cake = U401 - 0

    #determine temperature for recovery of solvent from miscella
    # bp=U401.outs[1].bubble_point_at_P()
    # U401.outs[1].T=bp.T
    H401 = bst.HXutility(lignin_extraction_area, ins=U401.outs[1], T=337.5)
    # #Desolvenization of miscella (oil&hexane): (14.6, Williams, 2013)
    if lipid_type == 'microbial':
        D401 = bst.units.BinaryDistillation(lignin_extraction_area, H401 - 0,
                                        outs=('distillate_hexane_water', 'bottom_oil'), LHK=('Hexane', 'TriOlein'),
                                        Lr=0.99, Hr=0.99, k=2, is_divided=True)
    elif lipid_type == 'vegetative':#vegetative lipid
        D401 = bst.units.BinaryDistillation(lignin_extraction_area, H401 - 0,
                                        outs=('distillate_hexane_water', 'bottom_oil'), LHK=('Hexane', 'TriOlein_veg'),
                                        Lr=0.99, Hr=0.99, k=2, is_divided=True)
    D401.check_LHK = False
    # miscella leaving the evaporator contains approximately 92–98% oil and 2–8% solvent (Le Clef & Kemper, 2015)
    # S401 = bst.Splitter(lignin_extraction_area,ins=D401-1, split=dict(Glucose=1,GlucoseOligomer=1,
    #                        GalactoseOligomer=1,))
    P402 = bst.Pump(lignin_extraction_area, ins = D401-1,P = 10000)
    # oil stripper/dryer to vacuum remove water and hexane (Le Clef & Kemper, 2015)
    F401 = bst.SplitFlash(lignin_extraction_area, T=273.15 + 40, P=10000.,  # 100 mbar
                          ins=P402 - 0,  # bottom_oil with minor solvent&water
                          outs=('',oil_from_lignin),
                          split=dict(Lipid=0.00001, Water=1.0, Hexane=1, AceticAcid=1, Furfural=1))
    M401 = bst.Mixer(lignin_extraction_area, ins=(D401-0, F401 - 0))
    H402 = bst.HXutility(lignin_extraction_area, ins=M401 - 0, V=0.0)
    # Hexane condenser: condense hexane & water mixture vapor from distillation & oil stripper
    # H402 = bst.HXutility(lignin_extraction_area, ins=M401 - 0, V=0.0,outs=hexane_recycled)# add recovered hexane to hexane storage tank
    P403 = bst.Pump(lignin_extraction_area, ins=H402-0, P=101325,outs=hexane_recycled) #adjust pressure before storage

    # #Desolvenization of solid biomass: toaster/dryer
    # U402 = bst.DrumDryer('U402',
    #         (U401-0, 'dryer_air', 'dryer_natural_gas'),
    #         ('', 'dryer_outlet_air', 'dryer_emissions'),
    #         moisture_content=0.40, split=dict(Hexane=1.0),
    #         utility_agent='Steam', T = 109+273.15
    #     )

    @mixer.add_specification(run=True)
    def adjust_hexane_to_extractor():
        if extractor_type == 'immersion':
            hexane_biomass_ratio = 5.0
        else:
            hexane_biomass_ratio = 1.0  # percolation (Pramparo et al., 2002)
        #mass of solids in saccarified biomass + mass of residual solidas in recycled hexane
        # target_hexane_biomass_ratio = hexane_biomass_ratio
        # current_hexane_biomass_ratio = ((ins[0].F_mass - ins[0].imass['Water']) + (mixer.ins[1].F_mass-mixer.ins[1].imass['Hexane']-mixer.ins[1].imass['Water']))/mixer.outs[0].imass['Hexane']
        target_hexane_kg = ((ins[0].F_mass - ins[0].imass['Water']) + (
                    mixer.ins[1].F_mass - mixer.ins[1].imass['Hexane'] - mixer.ins[1].imass[
                'Water'])) * hexane_biomass_ratio
        recycled_hexane_kg = mixer.ins[1].imass['Hexane']
        if abs(target_hexane_kg-recycled_hexane_kg)/target_hexane_kg > 1e-3: #0.1% error in flow rate
            if target_hexane_kg > recycled_hexane_kg:
                mixer.ins[0].imass['Hexane'] = target_hexane_kg - recycled_hexane_kg
            else:
                mixer.ins[0].imass['Hexane'] = 0
            # if missing_hexane_kg < 1e-3:  # negative value
            #     mixer.ins[1].imass['Hexane'] = 0
            # else:
            #     mixer.ins[0].imass['Hexane'] = missing_hexane_kg
        else:
            mixer.ins[0].imass['Hexane'] = 0
    #assume hexane:biomass = 1:1 using percolation extractor (Pramparo et al., 2002)