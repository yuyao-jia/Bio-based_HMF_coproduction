from pickle import FALSE

import thermosteam as tmo
import biosteam as bst
from biosteam import SystemFactory
from system_setup import stream as s
from system_setup import units_added as units
import flexsolve as flx
__all__ = (
    'hmf_furfural_purification_sys',
)
@SystemFactory(
    ID='hmf_furfural_purification_sys',
    ins=[s.pretreatment_liquor],
    outs=[s.purified_furfural,s.purified_hmf,s.concentrated_pretreatment_liqour,s.purification_wastewater,
          # s.water_condensate
          ],
)
def create_hmf_furfural_purification_system(ins, outs,
        nanofiltration_paramete=None,#setup later based on NF units from Lavanya
        purification_area=None,
        NF1_rejection_factors=None,
        NF2_rejection_factors=None,
    ):
    """
    Create a system that separates and purifies HMF & furfural from cellulosic pretreatment liquor.

    """
    pretreatment_liquor = ins
    purified_furfural,purified_hmf,concentrated_pretreatment_liquor,purification_wastewater = outs
    # purified_furfural,purified_hmf,concentrated_pretreatment_liquor,purification_wastewater,water_condensate = outs
    if purification_area is None: purification_area=500
    if NF1_rejection_factors is None: print('Please input rejection factors for 1st nanofiltration')
    if NF2_rejection_factors is None: print('Please input rejection factors for 2nd nanofiltration')

    H1 = bst.units.HXutility(purification_area, ins=ins[0], T=30+273.15)#cool pretreatment liquor to 30C for NF
#1st step nanofiltration to separate sugars and HMF&furfural apart
    N1 = units.Nanofiltration_oilcane_unit(purification_area, ins=H1-0, outs=(concentrated_pretreatment_liquor,'N1_permeate'),
                                                                          volume_reduction_factor=10,
                                                                          rejection_factors=NF1_rejection_factors)

    pump_N2 = bst.Pump(purification_area, ins=N1-1, P=101325*15)
    N2 = units.Nanofiltration_oilcane_unit(purification_area, ins = pump_N2-0,outs=('N2_retentate_furfural','N2_permeate_hmf') ,volume_reduction_factor=10,
                                           rejection_factors=NF2_rejection_factors)

    nanofiltration_membrane = bst.Stream('nanofiltration_membrane',polyamide=1, units='kg/hr')
    membrane_splitter = bst.ReversedSplitter(purification_area,ins=nanofiltration_membrane,outs=('nf_membrane1','nf_membrane2'))
    membrane_splitter.nanofiltration_membrane_lifetime = 5 #yrs, typical lifetime of nanofiltration plants

#N2.outs = ('retentate_furfural', 'permeate_hmf')

# #recovery of HMF
#HMF distillation based on Ralf Böhling et al., 2013
    PEG_feed = bst.Stream('PEG', PEG=2389.47, units='kg/hr') #high-boiling solvent
    PEG_recycled = bst.Stream('PEG_recycled',PEG=0.1, units='kg/hr')
    # PEG_recycle_mixer = bst.Mixer(purification_area,ins=(PEG_feed))
    PEG_tank = bst.StorageTank(purification_area, ins=PEG_feed)

    # F4 = bst.SplitFlash(purification_area,ins=N2-1,T=30+273.5,P=4000,split=dict(Water=0.9))
    # @F4.add_specification(run=False)
    # def adjust_F4():
    #     F4.run()
    #     F4.P = 101325
    # #ignore the vacuum system for simulation
    # check flash,remove vacuum too if too much vapor: checked, no need to remove vacuum for other flash
    # try removing the flash and directly go to the distillation column: tried, too much burden on flashes after distillation
    # F4b = bst.SplitFlash(purification_area,ins=F4-1,T=30+273.5,P=4000,split=dict(Water=0.6))
    # F4 = bst.MultiEffectEvaporator(purification_area,ins=N2-1,P=(4000,3000),V=0.5)
    # M1 = bst.units.Mixer(purification_area, ins=(F4-1, PEG_recycled))

    #add a distillation column to remove water instead of flash
    water_distillation = units.ShortcutColumn(purification_area, ins = N2-1,LHK=('Water','HMF'),Lr=0.95,Hr= 0.99,k=1.1, is_divided=True,P=4000,
                                  partial_condenser = False)
    M1 = bst.units.Mixer(purification_area, ins=(water_distillation-1, PEG_recycled))
    PEG_mixer = bst.Mixer(purification_area, ins=(M1-0, PEG_tank-0))
    #k = 1.1 is preferred, reflux ratio over minimum reflux ratio, at least over 1, but lower is desired if operating cost is high
    #if capital cost too high, want higher k value
    #Try 2 distillation columns: 1 for water, 1 for PEG
    D1 = units.ShortcutColumn(purification_area, ins = PEG_mixer-0,
                              LHK=('HMF','PEG'),Lr=0.97,Hr= 0.99,k=1.1, is_divided=True,P=4000,
                                  partial_condenser = False) #do not use partial condenser when under vacuum
    #pressure & yield based on patent
    D1.check_LHK=False
    PEG_recycle_pump = bst.Pump(purification_area, ins=D1-1,outs=PEG_recycled, P=101325)
    target_ratio_PEG = 0.05 #mass percentage of PEG in total mass of solution
#P&recovery based on paper, D1.outs=('distillate_hmf','bottoms_D1')
    @PEG_mixer.add_specification(run=True)
    def adjust_PEG():#to reach 5% PEG in solution fed into distillation (Ralf Böhling et al., 2013)
        M1_PEG = M1.outs[0].imass['PEG']
        M1_PEG_ratio =  M1.outs[0].get_mass_composition('PEG')
        mass_non_PEG = M1.ins[0].F_mass - M1.ins[0].imass['PEG'] + M1.ins[1].F_mass - M1.ins[1].imass['PEG']
        # if M1_PEG_ratio < 0.29999:
        # if (M1_PEG_ratio  - target_ratio_PEG)/target_ratio_PEG > 1e-3: #0.1% error
        #     # PEG_needed = 3/7*mass_non_PEG
        #     PEG_needed = 1/19*mass_non_PEG
        #     if PEG_needed - M1_PEG < 1e-3: #negative value
        #         PEG_mixer.ins[1].imass['PEG'] = 0
        #     else:
        #         PEG_mixer.ins[1].imass['PEG'] = PEG_needed - M1_PEG
        # else:
        #     PEG_mixer.ins[1].imass['PEG'] = 0
        if abs(M1_PEG_ratio  - target_ratio_PEG) > 1e-3: #0.1% error
            # PEG_needed = 3/7*mass_non_PEG
            PEG_needed = 1/19*mass_non_PEG
            if PEG_needed - M1_PEG < 1e-3: #negative value
                PEG_mixer.ins[1].imass['PEG'] = 0
            else:
                PEG_mixer.ins[1].imass['PEG'] = PEG_needed - M1_PEG
        else:
            PEG_mixer.ins[1].imass['PEG'] = 0

#crytallization to purify HMF based on Rilsager et al., 2013
#vacuum dry crude HMF
    # H2 = bst.HXutility(purification_area, ins=D1-0, V=0,rigorous=True) #condense vapor
    crude_hmf_pump = bst.Pump(purification_area, ins=D1-0, P=101325)
    recycle_hmf_mixer = bst.Mixer(purification_area, ins=(crude_hmf_pump-0))
    # F1 = bst.units.SplitFlash(purification_area, ins=recycle_hmf_mixer-0,
    #                       T=273.15+50,P=10000,split=dict(Water=0.9))#dried crude HMF
    # condensate_mixer = bst.Mixer(purification_area, ins=(F1-0, F4-0))
    # condensate_mixer = bst.Mixer(purification_area, ins=(water_distillation-0))
    # condenser = bst.HXutility(purification_area, ins=condensate_mixer-0, V=0,outs=water_condensate)
#F1.outs=('water_condensate','crude_HMF')
    # P1= bst.units.Pump(purification_area, ins=F1-1, P=101325)
    MTBE_feed = bst.Stream('MTBE', MTBE=5057.9, units='kg/hr') #organic solvent for crystallization
    # MTBE_mixer = bst.Mixer(purification_area, ins=(MTBE_feed))
    # MTBE_tank =  bst.StorageTank(purification_area,ins=MTBE_mixer-0)
    # M2 = bst.units.Mixer(purification_area, ins = (P1-0,MTBE_tank-0))#crude HMF in organic solvent for crystallization
    MTBE_recycled = bst.Stream('MTBE_recycled',MTBE=0.1, units='kg/hr')
    MTBE_tank = bst.StorageTank(purification_area, ins=MTBE_feed)
    # MTBE_recycle_mixer = bst.Mixer(purification_area, ins=(P1-0,MTBE_recycled))
    MTBE_recycle_mixer = bst.Mixer(purification_area, ins=(recycle_hmf_mixer-0,MTBE_recycled))
    MTBE_mixer = bst.Mixer(purification_area, ins=(MTBE_recycle_mixer-0,MTBE_tank-0))
    target_ratio_MTBE = 3 #ratio of L of MTBE per kg of crude HMF

    # U1 = units.HMFCrystallizer(purification_area,ins=M2-0,T=-30+273.15,solid_purity=0.98,target_recovery=0.90)
    # U2=units.HMFCrystalSeparator(purification_area,ins=U1-0,moisture_content=0,MTBE_content=0.2)
    U1 = units.HMFCrystallizerSeparator(purification_area,ins=MTBE_mixer-0)
    #wet_crystal, remain_liquid = U1.outs
    # F2 = bst.SplitFlash(purification_area,ins=U1-0,T=273.15+25,P=30000,split=dict(MTBE=0.999))#dry separated HMF crystal(solid)
    # H3 = bst.HXutility(purification_area, ins=U1-0, T=273.15+40)
    F2 = bst.Flash(purification_area,ins=U1-0,T=273.15+25,P=2000)#dry separated HMF crystal(solid)
    # def hmf_purity(pressure):
    #     F2.P = pressure
    #     F2._run()
    #     conc = F2.outs[1].get_mass_composition('HMF')
    #     targetconc = 0.99
    #     return conc - targetconc
    # @F2.add_specification(run=True)
    # def adjust_purity():
    #     p_guess = 6000
    #     F2.P = flx.IQ_interpolation(hmf_purity,x0=1000,x1=10000,x=p_guess,ytol=0.01,maxiter=1000)

    P2 = bst.Pump(purification_area, ins=F2-1, P=101325)
    T1 = bst.StorageTank(purification_area, ins=P2-0,outs=purified_hmf)

#recover & recycle MTBE from remaining liquid obtained from crystallization
    #temp based on boiling point at P + 5
    F3 = bst.units.SplitFlash(purification_area,ins=U1-1,T=273.15+24,P=10000,split=dict(MTBE=0.99))
    @F3.add_specification(run=False)
    def adjust_F3():
        F3.run()
        F3.P = 101325

#F3-1 is wastewater
    M3 = bst.units.Mixer(purification_area,ins=(F2-0,F3-0))#mix recovered MTBE
    P3 = bst.units.Pump(purification_area,ins=M3-0,P=101325,outs=MTBE_recycled)#pumping recovered MTBE to atmospheric pressure
    @P3.add_specification(run=True)
    def adjust_recovered_MTBE_phase():
        outlet = P3.outs[0]
        outlet.copy_like(P3.ins[0])
        outlet.phase = 'l'

    #recycle and adjust MTBE to reach target ratio
    @MTBE_mixer.add_specification(run=True)
    def adjust_MTBE():
        recycle_MTBE_L = MTBE_recycle_mixer.ins[1].ivol['MTBE']*1000
        target_MTBE_L = (MTBE_recycle_mixer.ins[0].F_mass + MTBE_recycle_mixer.ins[1].F_mass - MTBE_recycle_mixer.ins[1].imass['MTBE']) * target_ratio_MTBE
        if (target_MTBE_L - recycle_MTBE_L)/target_MTBE_L > 1e-3: #0.1% error in flow rate
            missing_MTBE_L = target_MTBE_L - recycle_MTBE_L
            missing_MTBE_kg = missing_MTBE_L*0.7404 #density of MTBE = 0.7404 kg/L
            MTBE_mixer.ins[1].imass['MTBE'] = missing_MTBE_kg
        else:
            MTBE_mixer.ins[1].imass['MTBE'] = 0
        # if abs(recycle_MTBE_L - target_MTBE_L) > 1e-3:
        #     missing_MTBE_L = target_MTBE_L - recycle_MTBE_L
        #     missing_MTBE_kg = missing_MTBE_L*0.7404 #density of MTBE = 0.7404 kg/L
        #     MTBE_mixer.ins[1].imass['MTBE'] = missing_MTBE_kg
        # else:
        #     MTBE_mixer.ins[1].imass['MTBE'] = 0
#recovery of furfural
#distillation to purifiy furfural based on Zeitsch, 2000 (some specific data) and Nhien et al., 2016 (configuration/design)
#distillation tray column to azeotropic vapor of furfural and water at 370.9K: lower than Tb of acetic acid and HMF (Nhien et al., 2016)
    # D2 = bst.units.BinaryDistillation(purification_area, ins = N2-0,
    #                               LHK=('Water','Furfural'),Lr=0.5,Hr= 0.001,k=2, is_divided=False,P=101325)#azeotropic column
    preheat= bst.units.HXutility(purification_area, ins=N2-0, T=70+273.15)#temp based on paper
    D2 = units.BinaryDistillation(purification_area, ins = preheat-0,
                                  LHK=('Furfural','Glucose'),Lr=0.99,Hr= 0.99,k=2, is_divided=False,P=10500, partial_condenser = False)#azeotropic column
    D2.check_LHK=False

    # @D2.add_specification(run=True)
    # def adjut_D2_P():
    #     try:
    #         D2._run()
    #     except:
    #         try:
    #             print('D2 failed to converge at 101325, trying 50000')
    #             D2 = bst.units.BinaryDistillation(purification_area, ins=preheat - 0,
    #                                       LHK=('Furfural', 'Glucose'), Lr=0.99, Hr=0.99, k=2, is_divided=False,
    #                                       P=50000)  # azeotropic column
    #             D2.check_LHK = False
    #             D2._run()
    #         except:
    #             print('D2 failed to converge at 50000, trying 30000')
    #             D2 = bst.units.BinaryDistillation(purification_area, ins=preheat - 0,
    #                                   LHK=('Furfural', 'Glucose'), Lr=0.99, Hr=0.99, k=2, is_divided=False,
    #                                   P=30000)  # azeotropic column
    #             D2.check_LHK = False
    #             D2._run()

    # D2.check_LHK=False
    #Lr=recovery of light key in distillate
    #Hr=recovery of heavy key in bottoms
    #D2-1:wastewater

    #deal with overall losses of furfural ayhgbnbv       hnbbt last step
    #65% water in azeotropic vapor (Zeitch.pdf)
    # def get_azeotropic_vapor_D2_water(Lr_test):
    #     D2.Lr = Lr_test
    #     D2.run()
    #     concentration=D2.outs[0].get_mass_composition('Water')
    #     target_concentration=0.65
    #     return concentration - target_concentration
    # def get_azeotropic_vapor_D2_furfural(Lr_test):
    #     D2.Lr = Lr_test
    #     D2.run()
    #     concentration=D2.outs[0].get_mass_composition('Furfural')
    #     target_concentration=0.35
    #     return concentration - target_concentration
    # @D2.add_specification(run=True)
    # def adjust_Lr_for_azeotropic_vapor_D2():
    #     Lr_guess = 0.5
    #     try:
    #         D2.Lr = flx.IQ_interpolation(get_azeotropic_vapor_D2_water,x0=0.01,x1=0.99,x=Lr_guess,ytol=0.05)
    #     except:
    #         print('D2 failed to adjust water at 0.05')
    #         try:
    #             D2.Lr = flx.IQ_interpolation(get_azeotropic_vapor_D2_furfural(), x0=0.01, x1=0.99, x=Lr_guess, ytol=0.075)
    #             # D2.Lr = flx.IQ_interpolation(get_azeotropic_vapor_D2_water, x0=0.01, x1=0.99, x=Lr_guess, ytol=0.075)
    #         except:
    #             print('D2 failed to adjust furfural at 0.075')
    #             try:
    #                 D2.Lr = flx.IQ_interpolation(get_azeotropic_vapor_D2_furfural(), x0=0.01, x1=0.99, x=Lr_guess,
    #                                              ytol=0.1)
    #             except:
    #                 print('D2 failed to adjust furfural at 0.1')
    #                 D2.Lr = 0.5
    # D2.outs=('distillate_vapor_D2','bottom_with_impurities')
# D2-0: azeotropic vapor of furfural and water(~65% of water) (Zeitsch, 2000)

#S1-1: wastewater
#mix all azeotropic vapor of furfural and water
    H3 = bst.HXutility(purification_area,ins=D2-0,V=0,rigorous=True)#condense vapor for decanter (Zeitsch, 2000)
    P5 = bst.Pump(purification_area,ins=H3-0,P=101325)
    # decanter1 = bst.LLECentrifuge(purification_area,ins=H3-0,top_chemical=('Water'),#Furfural is heavier than water
    decanter1 = bst.LLECentrifuge(purification_area, ins=P5 - 0, top_chemical=('Water'),
                                  # Furfural is heavier than water
                                  forced_split_IDs=('Furfural','Glucose','AceticAcid','Water','HMF'),
                                  forced_split=(0.01,1.0,0.3,0.9,0.99))
                                  # forced_split_IDs=('Glucose','Furfural'),forced_split=(1.0,0.01))
    [crude_hmf_pump-0, decanter1-0] - recycle_hmf_mixer  # recycle waste HMF from furfural purification to HMF purification
    # M6 = bst.Mixer(purification_area,ins=(H3-0))#mix with light phase from decanter2
    # decanter1 = bst.LLECentrifuge(purification_area,ins=H3-0,top_chemical=('Water'),#Furfural is heavier than water
    #                               forced_split_IDs=('Furfural','Water','AceticAcid','HMF'),forced_split=(0.01,0.9,0.35,1.0))#forced split to low density fluid (outs[0])
    # decanter1.water_split = 0.95#random number
    # def furfural_concentration(water_sp):
    #     decanter1.water_split = water_sp
    #     decanter1.forced_split = (0.01,decanter1.water_split)
    #     decanter1.run()
    #     conc = decanter1.outs[1].get_mass_composition('Furfural')
    #     targetconc = 0.94
    #     return conc - targetconc
    # # @decanter1.add_specification(run=True)
    # # def adjust_split_for_94_furfural():#adjust split of water to get 94% furfural in decanter1-0
    # #     split_guess = 0.5
    # #     water_sp = flx.IQ_interpolation(furfural_concentration,x0=0.01,x1=0.99,x=split_guess,ytol=0.01,maxiter=1000)
#decanter-1: furfural rich phase
#decanter-0: light water phase is wastewater
    D3 = units.BinaryDistillation(purification_area,ins=decanter1-1,LHK=('AceticAcid','Furfural'),
                                      Lr=0.9,Hr=0.99,k=2,
                                  is_divided=False,P=101325)#further dehydtrate and purify furfural
    D3.check_LHK=False
    # D3-1:purified furfural; D3-0:wastewater
    H5 = bst.HXutility(purification_area,D3-1,T=25+273.15)#cool furfural to room temp
    T2 = bst.StorageTank(purification_area, ins=H5-0,outs=(purified_furfural))#liquid furfural
    M8 = bst.Mixer(purification_area, ins=(D2 - 1, D3 - 0))
    H6 = bst.HXutility(purification_area, M8-0, V=0)
    P4 = bst.Pump(purification_area,ins=F3-1,P=101325)
    P5 = bst.Pump(purification_area,ins=water_distillation-0,P=101325)
    M9 = bst.Mixer(purification_area, ins=(H6 - 0, P4 - 0, P5-0),outs=purification_wastewater)