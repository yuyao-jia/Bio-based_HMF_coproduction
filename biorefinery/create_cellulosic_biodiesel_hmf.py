# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import thermosteam as tmo
import biosteam as bst
from biosteam import main_flowsheet as f
from biosteam import SystemFactory
from Unit_functions.Biodiesel_functions.lipid_prep_transesterification_funcs import (
    # create_lipid_pretreatment_system,
    create_lipid_pretreatment_system_with_degumming,
    create_transesterification_and_biodiesel_separation_system,
)
from Unit_functions.Fermentation_functions.fermentation_func import create_cane_to_combined_1_and_2g_fermentation
from Unit_functions.Biodiesel_functions.lipid_extraction_yj import (
    create_post_fermentation_oil_separation_system,
    create_lignin_oil_extraction_system
)
from system_setup import microbial_oil_baseline_yj as perf
from system_setup import units_added as units
from system_setup import stream as s
from Unit_functions.Facilities_functions.facilities_systems_yj import create_all_facilities
from system_setup.chemical import chem
from system_setup.reactions_conversion import pretreatment_rxn_wt, hydrolysis_rxn_wt

F = bst.Flowsheet('oilcane_baseline')
bst.main_flowsheet.set_flowsheet(F)
bst.settings.set_thermo(chem, cache=True)  # takes all chemicals and let biosteam know the system_setup we want to use

__all__ = (
    'create_oilcane_to_biodiesel_individual_oil_separation_microbial_vegetative',
)

'''

'''
@SystemFactory(
    ID='oilcane_biodiesel_sys',
    ins=[s.oilcane],
    outs=[s.biodiesel, s.crude_glycerol,s.purified_furfural,s.purified_hmf],
    fixed_outs_size=True,
)
def create_oilcane_to_biodiesel_individual_oil_separation_microbial_vegetative(
        ins, outs, fed_batch=True, WWT_kwargs=None,
        urea = True, #TODO:check whether urea is added in Yoel system
    ):
    oilcane, = ins
    biodiesel, crude_glycerol, purified_furfural,purified_hmf = outs
    if fed_batch:
        biomass_coeff = perf.fed_batch_biomass_growth_coefficient_mean
        lipid_yield = perf.fed_batch_lipid_yield_mean
        titer = perf.fed_batch_titer_mean
        productivity = perf.fed_batch_productivity_mean
    else:
        biomass_coeff = perf.batch_biomass_growth_coefficient_mean
        lipid_yield = perf.batch_lipid_yield_mean
        titer = perf.batch_titer_mean
        productivity = perf.batch_productivity_mean

    glucose_fermrxn = tmo.Rxn('O2 + Glucose -> H2O + TAG', 'Glucose', 1., correct_atomic_balance=True)
    glucose_fermrxn.product_yield('TAG', basis='wt', product_yield=lipid_yield)
    xylose_fermrxn = tmo.Rxn('O2 + Xylose -> H2O + TAG', 'Xylose', 1., correct_atomic_balance=True)
    xylose_fermrxn.product_yield('TAG', basis='wt', product_yield=lipid_yield)
    # arabinose_fermrxn = tmo.Rxn('O2 + Arabinose -> H2O + TAG', 'Arabinose', 1., correct_atomic_balance=True)
    # arabinose_fermrxn.product_yield('TAG', basis='wt', product_yield=lipid_yield)

    cellmass_rxn = tmo.Rxn(
        'Glucose + Urea -> Yeast + H2O + CO2', 'Glucose', 1., 
        correct_atomic_balance=True
    )
    cellmass_rxn.product_yield('Yeast', 'wt', biomass_coeff)
    combustion = tmo.Rxn('Glucose + O2 -> CO2 + H2O', 'Glucose', 1. - cellmass_rxn.X,
                         correct_atomic_balance=True)
    glucose_growrxn = cellmass_rxn + combustion
    glucose_growrxn.X = 0.999 - glucose_fermrxn.X
    
    cellmass_rxn = tmo.Rxn(
        'Xylose + Urea -> Yeast + H2O + CO2', 'Xylose', 1., 
        correct_atomic_balance=True
    )
    cellmass_rxn.product_yield('Yeast', 'wt', biomass_coeff)
    combustion = tmo.Rxn('Xylose + O2 -> CO2 + H2O', 'Xylose', 1. - cellmass_rxn.X,
                         correct_atomic_balance=True)
    xylose_growrxn = cellmass_rxn + combustion
    xylose_growrxn.X = 0.999 - xylose_fermrxn.X

    # cellmass_rxn = tmo.Rxn(
    #     'Arabinose + Urea -> Yeast + H2O + CO2', 'Arabinose', 1.,
    #     correct_atomic_balance=True)
    # cellmass_rxn.product_yield('Yeast', 'wt', biomass_coeff)
    # combustion = tmo.Rxn('Arabinose + O2 -> CO2 + H2O', 'Arabinose', 1. - cellmass_rxn.X,
    #                         correct_atomic_balance=True)
    # arabinose_growrxn = cellmass_rxn + combustion
    # arabinose_growrxn.X = 0.999 - arabinose_fermrxn.X

    cofermentation_rxn = tmo.PRxn(
        [glucose_fermrxn,
         xylose_fermrxn,
        # arabinose_fermrxn,
         glucose_growrxn,
         xylose_growrxn,
         # arabinose_growrxn
         ],
    )
    pretreatment_rxn_yj = pretreatment_rxn_wt.copy()  # wt basis
    hydrolysis_rxn_yj = hydrolysis_rxn_wt.copy()  # wt basis
    # Split values that are not 1 are from #2011 NREL report on cellulosic ethanol as given in [2]: Humbird?
    split_yj = {  # split pretreatment products to pretreatment liquor and solids
        'Glucose': 0.03647,
        'Xylose': 0.03766,
        'Ash': 1.0,
        'Lignin': 1.0,
        'Solids': 1.0,
        # assume all lipids stay in solid after pretreatment for easy modeling
        # Loss in pretreatment liquor is included after saccharification
        'TriOlein_veg': 1.0,
        'PL_veg': 1.0,
        'FFA_veg': 1.0,
        'AceticAcid': 0.03727,
        'Furfural': 0.03571,
        'HMF': 0.03571,
        'SolubleLignin': 0.03727,
        'GlucoseOligomer': 0.03722,
        'GalactoseOligomer': 0.03722,
        'XyloseOligomer': 0.03722,
        'ArabinoseOligomer': 0.03722,
        'Glucan': 1.0,
        'Xylan': 1.0,
        'Arabinan': 1.0,
        'Galactan': 1.0,
    }

    # based on experimental data (VRF10,pH adjusted)
    NF1_yj = {  # rejection factor for 1st Nanofiltration of pretreatment liquor (%mass balance based)
        'Glucose': 0.9839,
        'Xylose': 0.9818,
        'Arabinose': 0.9452,
        'AceticAcid': 0.8428,
        'Furfural': 0.5959,
        'HMF': 0.4072,
        # assume all sugar oligomers have the same rejection factor as cellobiose (1.0) because they are larger than cellobiose
        'GlucoseOligomer': 1.0,
        'GalactoseOligomer': 1.0,
        'XyloseOligomer': 1.0,
        'ArabinoseOligomer': 1.0,
        'SolubleLignin': 1.0,
    }

    # based on experimental data (VRF10,pH adjusted)
    NF2_yj = {
        # rejection factor for 2nd nanofiltration of permeate obtained from 1st nanofitlration (%mass balance based)
        'AceticAcid': 0.7754,
        'Furfural': 0.7882,
        # 'HMF': 0.3,  # lower number from Malmali 2014, can be used as lower bound of HMF rejection factor
        'HMF': 0.4755, #baseline
        'Glucose': 1.0,
        'Xylose': 1.0,
        'Arabinose': 1.0,
    }

    oilcane_to_fermentation_sys, ofs_dct = create_cane_to_combined_1_and_2g_fermentation('oilcane_to_fermentation_sys',
        ins=oilcane,
        outs=('','','','','',purified_furfural,purified_hmf,'','filter_cake'),
        # outs=('','','','','',purified_furfural,purified_hmf,'','','filter_cake'),
        product_group='Microbial_lipid',
        titer=titer,
        productivity=productivity,
        cofermentation_reactions=cofermentation_rxn,
        seed_train_reactions=bst.Rxn(None, 'Glucose', 1.), # Easier to simulate reactions only at cofermentation reactor
        #seed train's empty out stream is vent, but no reactions were simulated at seed_train, so no vent flow.
        CoFermentation=units.AeratedCoFermentation,
        SeedTrain=units.SeedTrain,
        include_scrubber=False,
        fed_batch=fed_batch,
        udct=True,
        mockup=False,
        add_urea = urea,
        pretreatment=pretreatment_rxn_yj,
        pretreatment_split=split_yj,
        NF1_rf = NF1_yj,
        NF2_rf = NF2_yj,
        saccharification_rxns = hydrolysis_rxn_yj,
    )
    beer, saccharified_biomass_solids, condensate, fiber_fines, bagasse_to_boiler, purified_furfural, purified_hmf, purification_wastewater,filter_cake = oilcane_to_fermentation_sys.outs
    # beer, saccharified_biomass_solids, condensate, fiber_fines, bagasse_to_boiler, purified_furfural, purified_hmf, purification_wastewater, water_condensate,filter_cake = oilcane_to_fermentation_sys.outs
    beer.ID = 'beer'
    # lignin_mixer = bst.Mixer(400,ins=lignin)
    hydrolysate_and_juice_mixer = bst.F.hydrolysate_and_juice_mixer
    post_fermentation_oil_separation_sys, pfls_dct = create_post_fermentation_oil_separation_system(
        ins=beer,
        mockup=False,
        area=500,
        udct=True,
    )
    microbial_lipid, cellmass_to_boiler, wastewater_micrbial_oil_extraction = post_fermentation_oil_separation_sys.outs
    microbial_lipid.ID = 'microbial_lipid'
    lignin_oil_extration_sys = create_lignin_oil_extraction_system(ID='veg_oil_extraction_sys',
                                                                        ins=saccharified_biomass_solids,
                                                                        outs=['vegetative_lipid','spent_lignin_to_boiler'],
                                                                        lignin_extraction_area=500,
                                                                        extractor_type='percolation',
                                                                        mockup=False,
                                                                        lipid_type = 'vegetative',)
    vegetative_lipid,spent_lignin_to_boiler = lignin_oil_extration_sys.outs
    oil_mixer = bst.Mixer(500, ins=[microbial_lipid, vegetative_lipid],outs='combined_extracted_lipid')
    @oil_mixer.add_specification(run=False)
    def convert_vegetative_lipid_back_to_normal():
        oil_mixer._run()
        combined = oil_mixer.outs[0]
        veg_lipids = ['PL_veg', 'FFA_veg', 'MAG_veg', 'DAG_veg', 'TAG_veg']
        microbial_lipid = ['PL', 'FFA', 'MAG', 'DAG', 'TAG']
        veg_imol = [oil_mixer.ins[1].imol[i] for i in veg_lipids]
        micro_imol = [oil_mixer.ins[0].imol[i] for i in microbial_lipid]
        lipid = [i+j for i,j in zip(veg_imol, micro_imol)]
        lipid_dict = dict(zip(microbial_lipid, lipid))
        for i in microbial_lipid:
            combined.imol[i] = lipid_dict[i]
        for i in veg_lipids: combined.imol[i] = 0

    wastewater_mixer = bst.Mixer(600,ins=(wastewater_micrbial_oil_extraction,purification_wastewater,))
    oil_pretreatment_sys, oil_pretreatment_dct = create_lipid_pretreatment_system_with_degumming(
        ins=oil_mixer-0,
        outs=['degummed_oil','polar_lipids'],
        mockup=False,
        area=500,
        udct=True
    )
    degummed_oil,polar_lipids = oil_pretreatment_sys.outs
    
    transesterification_and_biodiesel_separation_sys=create_transesterification_and_biodiesel_separation_system(
        ins=degummed_oil,
        outs=[biodiesel, crude_glycerol, 'wastewater_biodiesel'],
        mockup=False,
        area=500,
    )
    burner_mixer = bst.Mixer(700,ins=(bagasse_to_boiler,polar_lipids,filter_cake),outs='boiler_feed')
    s = f.stream
    u = f.unit
    if WWT_kwargs is None:
        WWT_kwargs = dict(area=600)
    else:
        WWT_kwargs['area'] = 600
    # WWT_kwargs['kind'] = 'highrate'
    WWT_kwargs['kind'] = 'conventional'
    # if WWT_kwargs['kind'] is 'highrate':
    #     AnMBR_kwargs = dict(Y_biomass=0.1) ##Biodegradability = 0.96
    #     # AnMBR_kwargs ['Y_biogas'] =0.9
    #     # AnMBR_kwargs = dict(Y_biomass=0.05) #used by Yalin Li
    #     # increased biodegradability of organic components due to low COD in inputs into AnMBR unit (Bokhary, Maleki, Hong, Hai, and Liao, 2020;Jensen et al., 2015)
    #     WWT_kwargs['AnMBR_kwargs'] = AnMBR_kwargs
    #     # WWT_kwargs['skip_AnMBR'] =True
    #     # WWT_kwargs['skip_IC'] = True
    #     # biodegradability = 1 #1 by Yalin Li
    #     # AnMBR_kwargs = dict(biodegradability=biodegradability)
    #     # IC_kwargs = dict(biodegradability=biodegradability)
    #     # WWT_kwargs['AnMBR_kwargs'] = AnMBR_kwargs
    #     # WWT_kwargs['IC_kwargs'] = IC_kwargs
    #     # WWT_kwargs['skip_AeF'] = True

    # HXN_kwargs = {}
    # HXN_kwargs['ID'] = 800
    #TODO: testing which units to ignore by areas
    #Progess so far:
    # works if all units are ignored, doesn't work if only boilers are included
    #works if only area100 is ignored:
        #works if only H101 or H102 is ignored, so they need to be both ignored?
    #works if only area300 is ignored:
        #works if H301 to H303 are ignored: checking which one is the problem?
            #doesn't work if only H302 is ignored
        #doesn't work if H304 to H307 are ignored: H304 to H307 doesn't need to be ignored
    #doesn't work if only area400 or 500 is ignored: doesn't need to be ignored
    # HXN_kwargs['ignored'] = [u.D301.boiler, u.D302.boiler,u.D303.boiler,u.D501.boiler,u.D502.boiler,u.D503.boiler,u.D504.boiler,
    #                          u.H101,
    #                          u.H102,
    #                          u.H301,
    #                          u.H302,
    #                          u.H303,
    #                          u.H304,u.H305,u.H306,u.H307,
    #                          # u.H401,u.H402,
    #                          # u.H501,u.H502,u.H503,u.H504,u.H505, #doesn't work if these are ignored
    #                          u.H506,u.H507,u.H508,u.H509,u.H510,u.H511
    #                          ]
    facilities = create_all_facilities(
        feedstock=s.oilcane,
        # recycle_process_water_streams=(condensate, water_condensate),
        recycle_process_water_streams=(condensate,),
        WWT_kwargs=WWT_kwargs,
        # HXN_kwargs=HXN_kwargs,
        HXN = None,
        area=700,
    )