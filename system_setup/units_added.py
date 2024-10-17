# -*- coding: utf-8 -*-
"""
.. contents:: :local:

Reactors
--------
.. autoclass:: biorefineries.cane.units.SeedTrain
.. autoclass:: biorefineries.cane.units.CoFermentation

Separations
-----------
.. autoclass:: biorefineries.cane.units.OleinCrystallizer
    

"""
import biosteam as bst
from biorefineries.cellulosic.units import (
    SeedTrain,
    CoFermentation
)
from biosteam.units.decorators import cost
from biosteam.units.design_tools.cost_index import CEPCI_by_year
from thermosteam import PRxn, Rxn

from system_setup.units_cellulosic_yoel import *
from system_setup.units_distillation_added import *

__all__ = ('SeedTrain', 'CoFermentation', 'HMFCrystallizer', 'Nanofiltration_oilcane_unit',
           'HMFCrystalSeparator', 'HMFSolidsDryer', 'HMFCrystallizerSeparator')


class SeedTrain(SeedTrain):

    def __init__(self, ID='', ins=None, outs=(), thermo=None, reactions=None, saccharification=False):
        bst.Unit.__init__(self, ID, ins, outs, thermo)
        self.saccharification = saccharification
        chemicals = self.chemicals
        self.reactions = reactions or PRxn([
            #   Reaction definition                   Reactant    Conversion
            Rxn('Glucose -> 2 Ethanol + 2 CO2', 'Glucose', 0.9000, chemicals),
            Rxn('3 Xylose -> 5 Ethanol + 5 CO2', 'Xylose', 0.8000, chemicals),
            Rxn('Glucose -> Cellmass', 'Glucose', 0.0473, chemicals, correct_mass_balance=True),
            Rxn('Xylose -> Cellmass', 'Xylose', 0.0421, chemicals, correct_mass_balance=True),
        ])

    def _setup(self):
        super()._setup()
        self.outs[0].phase = 'g'

    def _run(self):
        vent, effluent = self.outs
        effluent.mix_from(self.ins, energy_balance=False)
        self.reactions.force_reaction(effluent)
        effluent.mol.remove_negatives()
        effluent.T = self.T
        vent.empty()
        vent.copy_flow(effluent, ('CO2', 'O2', 'N2'), remove=True)


class CoFermentation(CoFermentation):

    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 tau=36, N=None, V=3785.4118, T=305.15, P=101325,
                 Nmin=2, Nmax=36, cofermentation=None):
        bst.BatchBioreactor.__init__(self, ID, ins, outs, thermo, tau, N, V, T, P, Nmin, Nmax)
        self.P = P
        chemicals = self.chemicals
        self.loss = None
        self.cofermentation = cofermentation or PRxn([
            #   Reaction definition                   Reactant    Conversion
            Rxn('Glucose -> 2 Ethanol + 2 CO2', 'Glucose', 0.9500, chemicals),
            Rxn('3 Xylose -> 5 Ethanol + 5 CO2', 'Xylose', 0.8500, chemicals),
            Rxn('Glucose -> Cellmass', 'Glucose', 0.05, chemicals, correct_mass_balance=True),
            Rxn('Xylose -> Cellmass', 'Xylose', 0.05, chemicals, correct_mass_balance=True),
        ])

        if 'CSL' in chemicals:
            self.CSL_to_constituents = Rxn(
                'CSL -> 0.5 H2O + 0.25 LacticAcid + 0.25 Protein', 'CSL', 1.0000, chemicals, basis='wt',
            )
            self.CSL_to_constituents.basis = 'mol'
        else:
            self.CSL_to_constituents = None
        if all([i in self.chemicals for i in ('FFA', 'DAG', 'TAG', 'Glycerol')]):
            self.lipid_reaction = self.oil_reaction = PRxn([
                Rxn('TAG + 3Water -> 3FFA + Glycerol', 'TAG', 0.23, chemicals),
                Rxn('TAG + Water -> FFA + DAG', 'TAG', 0.02, chemicals)
            ])
        else:
            self.lipid_reaction = self.oil_reaction = None


class AeratedCoFermentation(bst.AeratedBioreactor):  # For microbial oil production
    V_max_default = 500  # m3 #From Book: Industrial Microbiology and Biotechnology

    # Verma, P.(Ed.), 2022.Industrial Microbiology and Biotechnology.Springer. Singapore, Singapore.https: // doi.org / 10.1007 / 978 - 981 - 16 - 5214 - 1

    def __init__(
            self, ID='', ins=None, outs=(), thermo=None,
            *, cofermentation, theta_O2=0.5,
            dT_hx_loop=8,
            Q_O2_consumption=-460240,
            batch=True,
            # [kJ/kmol] equivalent to 110 kcal / mol as in https://www.academia.edu/19636928/Bioreactor_Design_for_Chemical_Engineers
            **kwargs,
    ):
        bst.StirredTankReactor.__init__(self, ID, ins, outs, thermo, batch=batch, dT_hx_loop=dT_hx_loop, **kwargs)
        chemicals = self.chemicals
        self.theta_O2 = theta_O2
        self.hydrolysis_reaction = Rxn('Sucrose + Water -> 2Glucose', 'Sucrose', 1.00, chemicals)
        self.cofermentation = cofermentation
        self.lipid_reaction = self.oil_reaction = PRxn([
            Rxn('TAG + 3Water -> 3FFA + Glycerol', 'TAG', 0.23, chemicals),
            Rxn('TAG + Water -> FFA + DAG', 'TAG', 0.02, chemicals)
        ])
        self.Q_O2_consumption = Q_O2_consumption
        self.optimize_power = True
        self.kLa_coefficients = "Van't Riet"

    def _run_vent(self, vent, effluent):
        vent.copy_flow(effluent, ('CO2', 'O2', 'N2'), remove=True)
        assert not effluent.imol['CO2', 'O2', 'N2'].any()

    def run_reactions(self, effluent):
        self.hydrolysis_reaction.force_reaction(effluent)
        self.lipid_reaction.force_reaction(effluent)
        if effluent.imol['H2O'] < 0.: effluent.imol['H2O'] = 0.
        self.cofermentation.force_reaction(effluent)


class AeratedFermentation(bst.AeratedBioreactor):  # For microbial oil production
    V_max_default = 500

    def __init__(
            self, ID='', ins=None, outs=(), thermo=None,
            *, fermentation_reaction, cell_growth_reaction, theta_O2=0.5,
            dT_hx_loop=8,
            Q_O2_consumption=-460240,
            # [kJ/kmol] equivalent to 110 kcal / mol as in https://www.academia.edu/19636928/Bioreactor_Design_for_Chemical_Engineers
            **kwargs,
    ):
        bst.StirredTankReactor.__init__(self, ID, ins, outs, thermo, dT_hx_loop=dT_hx_loop, **kwargs)
        chemicals = self.chemicals
        self.theta_O2 = theta_O2
        self.hydrolysis_reaction = Rxn('Sucrose + Water -> 2Glucose', 'Sucrose', 1.00, chemicals)
        self.fermentation_reaction = fermentation_reaction
        self.cell_growth_reaction = cell_growth_reaction
        self.lipid_reaction = self.oil_reaction = PRxn([
            Rxn('TAG + 3Water -> 3FFA + Glycerol', 'TAG', 0.23, chemicals),
            Rxn('TAG + Water -> FFA + DAG', 'TAG', 0.02, chemicals)
        ])
        self.Q_O2_consumption = Q_O2_consumption
        self.optimize_power = False

    def _run_vent(self, vent, effluent):
        vent.receive_vent(effluent, energy_balance=False, ideal=True)

    def run_reactions(self, effluent):
        self.hydrolysis_reaction.force_reaction(effluent)
        self.lipid_reaction.force_reaction(effluent)
        if effluent.imol['H2O'] < 0.: effluent.imol['H2O'] = 0.
        self.fermentation_reaction.force_reaction(effluent)
        self.cell_growth_reaction.force_reaction(effluent)


class HMFCrystallizer(bst.BatchCrystallizer):

    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                 T,
                 solid_purity=0.99,
                 target_recovery=0.90,  # (% recovery of HMF)
                 order=None):
        bst.BatchCrystallizer.__init__(self, ID, ins, outs, thermo,
                                       tau=12, V=1e6, T=T)
        self.solid_purity = solid_purity
        self.target_recovery = target_recovery

    def set_effluent_composition_from_recovery(self, recovery, solid_pure):
        in_stream = self.ins[0]
        other_mass = in_stream.imass['AceticAcid'] + in_stream.imass['H2O']
        acetic_acid_fraction = in_stream.imass['AceticAcid'] / other_mass
        water_fraction = in_stream.imass['H2O'] / other_mass
        self.outs[0].imass['s', 'HMF'] = recovery * in_stream.imass['HMF']
        self.outs[0].imass['l', 'HMF'] = (1 - recovery) * in_stream.imass['HMF']
        nonpurity_solid = self.outs[0].imass['s', 'HMF'] / solid_pure * (1 - solid_pure)
        self.outs[0].imass['s', 'AceticAcid'] = nonpurity_solid * acetic_acid_fraction
        self.outs[0].imass['s', 'H2O'] = nonpurity_solid * water_fraction

    # @property
    # def Hnet(self):
    #     feed = self.ins[0]
    #     effluent = self.outs[0]
    #     if 's' in feed.phases:
    #         solid = feed if feed.phase == 's' else feed['s']
    #         H_in = - sum([i.Hfus * j for i,j in zip(self.chemicals, solid.mol) if i.Hfus])
    #     else:
    #         H_in = 0.
    #     solids = effluent['s']
    #     H_out = - sum([i.Hfus * j for i,j in zip(self.chemicals, solids.mol) if i.Hfus])
    #     return H_out - H_in

    def _run(self):
        feed = self.ins[0]
        outlet = self.outs[0]
        target_recovery = self.target_recovery
        outlet.copy_like(feed)
        outlet.phases = ('s', 'l')
        outlet.T = self.T
        solid_purity = self.solid_purity
        self.set_effluent_composition_from_recovery(target_recovery, solid_purity)


class HMFCrystalSeparator(bst.PressureFilter):
    _N_ins = 1

    # _units = {'Solids flow rate': 'kg/hr'}

    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                 MTBE_content=0.1, moisture_content=None, split=None):
        bst.PressureFilter.__init__(self, ID, ins, outs, thermo,
                                    moisture_content=moisture_content,
                                    split=split)
        self.split = None
        self.MTBE_content = MTBE_content

    def _run(self):
        feed = self.ins[0]
        solids = self.outs[0]
        solids.phases = ('s', 'l')
        liquids = self.outs[1]
        MTBE_content = self.MTBE_content
        solids.imass['s', 'HMF'] = feed.imass['s', 'HMF']
        solids.imass['s', 'AceticAcid'] = feed.imass['s', 'AceticAcid']
        solids.imass['s', 'H2O'] = feed.imass['s', 'H2O']
        solids.imass['l', 'MTBE'] = solids.F_mass * MTBE_content
        # some liquid solvent were separated together with solids (stay liquid)
        liquids.imass['HMF'] = feed.imass['l', 'HMF']
        liquids.imass['AceticAcid'] = feed.imass['l', 'AceticAcid']
        liquids.imass['PEG'] = feed.imass['PEG']
        liquids.imass['Furfural'] = feed.imass['Furfural']
        liquids.imass['MTBE'] = feed.imass['MTBE'] - solids.imass['MTBE']
        liquids.imass['H2O'] = feed.imass['H2O'] - solids.imass['H2O']
        solids.T = liquids.T = feed.T


class HMFSolidsDryer(bst.Flash):
    line = 'Flash'

    def __init__(self, ID='', ins=None, outs=(), thermo=None, *, split,
                 order=None, T=None, P=None, Q=None,
                 vessel_material='Carbon steel',
                 vacuum_system_preference='Liquid-ring pump',
                 has_glycol_groups=False,
                 has_amine_groups=False,
                 vessel_type=None,
                 holdup_time=15,
                 surge_time=7.5,
                 has_mist_eliminator=False):
        bst.Splitter.__init__(self, ID, ins, outs, thermo, split=split, order=order)
        self._load_components()

        self.T = T  #: Operating temperature (K)
        self.P = P  #: Operating pressure (Pa)
        self.Q = Q  #: Duty (kJ/hr)

        #: [str] Vessel construction material
        self.vessel_material = vessel_material

        #: [str] If a vacuum system is needed, it will choose one according to this preference.
        self.vacuum_system_preference = vacuum_system_preference

        #: [bool] True if glycol groups are present in the mixture
        self.has_glycol_groups = has_glycol_groups

        #: [bool] True if amine groups are present in the mixture
        self.has_amine_groups = has_amine_groups

        #: [str] 'Horizontal', 'Vertical', or 'Default'
        self.vessel_type = vessel_type

        #: [float] Time it takes to raise liquid to half full (min)
        self.holdup_time = holdup_time

        #: [float] Time it takes to reach from normal to maximum liquied level (min)
        self.surge_time = surge_time

        #: [bool] True if using a mist eliminator pad
        self.has_mist_eliminator = has_mist_eliminator

    split = bst.Splitter.split
    V = None

    def _run(self):
        top, bot = self.outs
        feed, = self.ins
        feed_mol = feed.mol
        top.mol[:] = top_mol = feed_mol * self.split
        bot.mol[:] = feed_mol - top_mol
        top.phase = 'g'
        bot.phase = 's'
        bot.T = top.T = self.T
        bot.P = top.P = self.P

    def _design(self):
        self.heat_exchanger.simulate_as_auxiliary_exchanger(self.ins, self.outs, vle=False)
        super()._design()


class HMFCrystallizerSeparator(bst.Unit, isabstract=True):
    _N_ins = 1
    _N_outs = 2

    auxiliary_unit_names = ('crystallizer',
                            'separator',)

    # 'crystal_dryer')

    def __init__(self, ID='', ins=(), outs=(),
                 thermo=None, T=-30 + 273.15, target_recovery=0.90, moisture_content=None, MTBE_content=0.1):
        bst.Unit.__init__(self, ID, ins, outs, thermo)
        self.crystallizer = crystallizer = HMFCrystallizer(None, ins='crude_HMF_MTBE_mixture',
                                                           T=T, solid_purity=0.99, target_recovery=target_recovery)
        self.separator = separator = HMFCrystalSeparator(None, ins='crystallized_mixture',
                                                         moisture_content=moisture_content, MTBE_content=MTBE_content)
        # self.crystal_dryer = crystal_dryer = bst.SplitFlash(None,ins='separated_crystal',T=273.15+25,
        #                                                     P=30000,split=split)#dry separated HMF crystal(solid)

    def _run(self):
        feed = self.ins[0]
        wet_crystal, remain_liquid = self.outs
        self.crystallizer.ins[0].copy_like(feed)
        self.crystallizer._setup()
        self.crystallizer._run()

        self.separator.ins[0].copy_like(self.crystallizer.outs[0])
        self.separator._setup()
        self.separator._run()
        # crystal.copy_like(self.separator.outs[0])
        remain_liquid.copy_like(self.separator.outs[1])
        wet_crystal.copy_like(self.separator.outs[0])
        wet_crystal.phase = 'l'  # easier to simulate later, converts back to solid later
        # self.crystal_dryer.ins[0].copy_like(self.separator.outs[0])
        # self.crystal_dryer._setup()
        # self.crystal_dryer._run()
        # liquid_from_flash.copy_like(self.crystal_dryer.outs[0])
        # dryed_crystal.copy_like(self.crystal_dryer.outs[1])

    def _design(self):
        self.crystallizer._design()
        self.separator._design()
        # self.crystal_dryer._design()

    def _cost(self):
        self.crystallizer._cost()
        self.separator._cost()
        # self.crystal_dryer._cost()


# Lavanya Sep 26th 2023
# Ref: Rules of Thumb
# Molar mass cut off = 200 = 0.2 kDa. The usual range is 0.01–1 kDa.
# Pressure: (between UF and RO) = 0.3–1.4 MPa.
# @cost(
#     ID='Nanofiltration membrane',
#     basis='Total filtration area',
#     units='m^2',
#     cost=20.8,  # doi:10.1016/j.desal.2005.08.030
#     CE=CEPCI_by_year[2001],  # https://toweringskills.com/financial-analysis/cost-indices/
#     n=1,
#     S=1,
#     lifetime=5  # typical lifetime of nanofiltration plants
# )
# Below three are based on cost modelling of ultrafiltration membranes
# https://doi.org/10.1089/ees.2000.17.61
@cost(
    ID='Instruments and controls',
    basis='Total filtration area',
    units='m^2',
    cost=1445.50,
    CE=CEPCI_by_year[2000],
    n=0.66,
    S=1,
)
@cost(
    ID='Tanks and frames',
    basis='Total filtration area',
    units='m^2',
    cost=3047.21,
    CE=CEPCI_by_year[2000],
    n=0.53,
    S=1,
)
# The miscellaneous
# component includes process equipment building, electrical supply and distribution,
# disinfection facilities, treated water storage and pumping, and wash-water recovery system
@cost(
    ID='Miscellaneous nanofiltration expenses',
    basis='Total filtration area',
    units='m^2',
    cost=7865.02,
    CE=CEPCI_by_year[2000],
    n=0.57,
    S=1,
)
# varies between 0.5 to 1 for a proper cleaned system, can go up to 4$ per 1000 gallons
# cost per gal is 0.5/1000 gal
# 1m3 is 265 US gallons
# https://samcotech.com/cost-to-properly-maintain-membrane-filtration-systems/

@cost(
    ID='Membrane cleaning and regeneration',
    basis='Throughput',
    units='m^3',
    cost=0.5 / 1000, CE=CEPCI_by_year[2021],
    n=1,
    S=1,
)
@cost(
    ID='Antiscalants_and_filters',
    basis='Throughput',
    units='m^3',
    cost=0.01 / 1000, CE=CEPCI_by_year[2021],
    n=1,
    S=1,
)
class Nanofiltration_oilcane_unit(bst.Unit, isabstract=True):
    _N_ins = 1
    _N_outs = 2
    _units = {
        'Total filtration area': 'm^2',
        'Throughput': 'm^3'
    }

    def __init__(self, ID='', ins=(), outs=(), thermo=None,
                 rejection_factors=None,  # dict {'chemical_name',rejection_fraction}
                 volume_reduction_factor=None  # provide a number to account for reduction in amount of water
                 ):
        bst.Unit.__init__(self, ID, ins, outs, thermo)
        self.rejection_factors = rejection_factors
        self.volume_reduction_factor = volume_reduction_factor
        self.feed_rate_L_per_sec = 0.5  # based on L/s.m^2 provided in the book, range is between 0.25 to 0.5

    def _run(self):
        feed = self.ins[0]
        retenate, permeate, = self.outs
        if len(self.rejection_factors) == len(feed.available_chemicals) - 1:
            for i in self.rejection_factors:
                mass_water = feed.imass['Water']
                temp_chem_mass = feed.imass[i]
                retenate.imass[i] = temp_chem_mass * (self.rejection_factors[i])
                permeate.imass[i] = temp_chem_mass * (1 - self.rejection_factors[i])
                retenate.imass['Water'] = r_water_mass = mass_water / self.volume_reduction_factor
                permeate.imass['Water'] = mass_water - r_water_mass
        else:
            print('Please provide all the rejection factors except water')

    def _design(self):
        m3_per_h_to_L_per_s = 0.2777
        m3_per_h_to_gal_per_h = 265
        process_vol_rate = self.ins[0].F_vol * m3_per_h_to_L_per_s  # provides value in L/s
        self.design_results['Total filtration area'] = process_vol_rate / self.feed_rate_L_per_sec
        self.design_results['Throughput'] = self.ins[0].F_vol / m3_per_h_to_gal_per_h


# # Only modified to convert vegetative lipids (separated for easy calculation of lipid recvoery) back to normal lipids
# class GlycerolysisReactor_modified(bst.CSTR):
#     _ins_size_is_fixed = False
#     _N_ins = 2
#     _N_outs = 2
#     T_default = 273.15 + 230
#     P_default = 101325
#     tau_default = 2
#
#     def _setup(self):
#         super()._setup()
#         # End result is near 100% conversion of FFAs with selectivity
#         # of 42% MAG, 47% DAG, 11% TAG:
#         # [1] Erik Anderson. Superior Process Technologies.
#         # Glycerolysis for Lowering Free Fatty Acid Levels
#
#         # [2] Kapil Mamtani, KavehShahbaz Mohammed M.Farid.
#         # Glycerolysis of free fatty acids: A review.
#         # https://doi.org/10.1016/j.rser.2020.110501
#
#         self.glycerolysis_baseline = PRxn([
#             Rxn('FFA + Glycerol -> MAG + H2O', reactant='FFA', X=1.00),
#             Rxn('DAG + Glycerol -> 2MAG', reactant='DAG', X=1.00),
#             Rxn('TAG + 2Glycerol -> 3MAG', reactant='TAG', X=1.00),
#             Rxn('FFA_veg + Glycerol -> MAG + H2O', reactant='FFA_veg', X=1.00),
#             Rxn('DAG_veg + Glycerol -> 2MAG', reactant='DAG_veg', X=1.00),
#             Rxn('TAG_veg + 2Glycerol -> 3MAG', reactant='TAG_veg', X=1.00),
#         ])
#         self.glycerolysis = PRxn([
#             Rxn('2MAG -> DAG + Glycerol', reactant='MAG', X=0.49),
#             Rxn('3MAG -> TAG + 2Glycerol', reactant='MAG', X=0.13),
#         ])
#
#         # Alternative preliminary modeling for legacy purposes:
#         # self.glycerolysis = ParallelReaction([
#         #   Reaction('FFA + Glycerol -> MAG + H2O', reactant='FFA',  X=0.80),
#         #   Reaction('2FFA + Glycerol -> DAG + 2H2O', reactant='FFA',  X=0.15),
#         #   Reaction('3FFA + Glycerol -> TAG + 3H2O', reactant='FFA',  X=0.05),
#         # ])
#
#     def _run(self):
#         feed, *other, N2 = self.ins
#         vent, effluent = self.outs
#         vent.P = effluent.P = self.P
#         vent.T = effluent.T = self.T
#         vent.phase = 'g'
#         effluent.mix_from(self.ins, energy_balance=False)
#         self.glycerolysis_baseline.force_reaction(effluent)
#         self.glycerolysis(effluent)
#         effluent.mol.remove_negatives()  # The correct glycerol flow rate is taken care of in a unit specification
#         vent.copy_flow(effluent, ('N2', 'Water'), remove=True)



