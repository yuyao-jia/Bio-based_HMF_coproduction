# -*- coding: utf-8 -*-
"""
"""

import biosteam as bst
from thermosteam import Stream
from system_setup import stream as s
from biorefineries.cellulosic import units


__all__ = (
    'create_hot_water_pretreatment_system',
)

@bst.SystemFactory(
    ID='hot_water_pretreatment_sys',
    ins=[s.bagasse],
    outs=[s.pretreated_biomass,
          s.pretreatment_liquor],
)
def create_hot_water_pretreatment_system(
        ins, outs,
        pretreatment_area=200,
        solids_loading=0.305,
        nonsolids=['Water'],
        milling=False,
        T_pretreatment_reactor=273.15 + 210,
        pretreatment_reactions = None,
        split_ratio=None,
    ):
    
    feedstock, = ins
    pretreated_biomass, pretreatment_liquor = outs
    
    warm_process_water = Stream('warm_process_water',
                              T=368.15,
                              P=4.7*101325,
                              Water=1)
    pretreatment_steam = Stream('pretreatment_steam',
                    phase='g',
                    T=268 + 273.15,
                    P=13 * 101325,
                    Water=24534+3490,
                    units='kg/hr')
    
    ### Pretreatment system
    n = pretreatment_area
    P = pretreatment_steam.chemicals['H2O'].Psat(T_pretreatment_reactor + 25)
    M203 = bst.SteamMixer(f'M{n+2}', (feedstock, pretreatment_steam, warm_process_water),
                          P=P, solids_loading=solids_loading)
    R201 = units.PretreatmentReactorSystem(f'R{n+1}', M203-0, T=T_pretreatment_reactor,reactions=pretreatment_reactions,tau=5/60)
    P201 = bst.Pump(f'P{n+1}', R201-1, P=101325)
    # F201 = units.PretreatmentFlash(f'F{n+1}', P201-0, P=101325, Q=0)
    # S201 = bst.SolidsSeparator(f'S{n+1}', P201-0,split=split_ratio,moisture_content=0.5) #outs[0]=retentate, outs[1]=permeate
    S201 = bst.PressureFilter(f'S{n+1}', P201-0, split=split_ratio, moisture_content=0.5)
    M204 = bst.Mixer(f'M{n+3}', (R201-0, S201-1),pretreatment_liquor)
    # units.WasteVaporCondenser(f'H{n+1}', M204-0, pretreatment_wastewater, V=0)
    P202 = units.HydrolyzatePump(f'P{n+2}', S201-0, None)
    if milling:
        bst.HammerMill(f'U{n+1}', P202-0, pretreated_biomass)
    else:
        P202.outs[0] = pretreated_biomass
