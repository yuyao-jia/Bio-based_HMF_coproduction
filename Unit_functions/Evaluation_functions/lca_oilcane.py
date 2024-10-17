from copy import copy as copy_
import thermosteam as tmo
from system_setup.lca_characterization_factors_yj import (
    GWP_characterization_factors as GWP_factors,
)
import biosteam as bst
# from model_uncertainty_analysis import (baseline_sys, chem, feedstock, area_groups,
#                                         products,product_streams_chemical_IDs, product_stream_names)

# products = [purified_hmf, biodiesel, purified_furfural, crude_glycerol]
# product_streams_chemical_IDs = ['HMF', 'Biodiesel', 'Furfural', 'Glycerol']
# product_stream_names = ['HMF', 'Biodiesel', 'Furfural', 'Crude glycerol'] #for allocation factor & GWP calculation

class OilcaneLCA:
    def __init_subclass__(cls, isabstract=False):
        if isabstract: return
        for method in ('_DPI', '_TDC', '_FCI', '_FOC'):
            if not hasattr(cls, method):
                import pdb
                pdb.set_trace()
                raise NotImplementedError(
                    f"subclass must implement a '{method}' method unless the "
                     "'isabstract' keyword argument is True"
                )
    @staticmethod
    def like(system, other):
        """Create an LCA object from `system` with the same settings as `other`."""
        self = copy_(other)
        self.units = sorted([i for i in system.units if i._design or i._cost], key=lambda x: x.line)
        self.system = system
        self.feeds = system.feeds
        self.products = system.products
        system._TEA = self
        return self

    def __init__(self, system, system_chemicals, area_groups,
                 feedstock,
                 main_product, main_product_chemical_IDs, main_product_names,
                 add_EOL_GWP=False,
                 ):

        #: [System] System being evaluated.
        self.system = system
        self.chem = system_chemicals
        self.units = self.system.units
        self.flowsheet = self.system.flowsheet
        self.streams = self.flowsheet.stream
        self.tea = self.system.TEA
        self.GWP = 'GWP100'
        self.BT = self.flowsheet.BT701
        self.area_groups = area_groups

        self.add_EOL_GWP = add_EOL_GWP
        # self.feedstock_mass_FU_kind = feedstock_mass_FU_kind

        # self.feedstock_ID = feedstock_ID
        # self.priced_feeds = [i for i in feeds if i.price]
        # self.priced_emissions = [i for i in emissions if i.price]

        # self.CFs = CFs
        # self.demand_allocation_method = demand_allocation_method

        # self._chemical_IDs = [chemi.ID for chemi in self.chem]

        # self.GWP_CF_stream = CFs['GWP_CF_stream']

        self.feedstock = feedstock #main feedstock stream, ex. oilcane
        self.main_products = main_product #sreams
        self.main_product_chemical_ID = main_product_chemical_IDs #to calculate allocation factors
        self.main_product_names = main_product_names
        # self.by_products = by_products

        tmo.settings.set_thermo(self.chem)
        # self.LCA_stream = Stream('LCA_stream', units='kg/hr')
        # self.LCA_streams = [i for i in system.feeds if not i == feedstock]
        # self.chem_IDs = [i.ID for i in chem]

        # self.conc_CO2_sequestered_in_liquid_waste_streams = conc_CO2_sequestered_in_liquid_waste_streams
        system._LCA = self
    #Materials-related impacts
    @property
    def feeds(self):#all feed streams
        return self.system.feeds
    @property
    def annual_factor(self):
        return self.tea.operating_days * 24
    @property
    def feedstock_GWP(self): #for main feedstock, ex. oilcane
        return self.system.get_material_impact(self.feedstock, self.GWP)
    @property
    def natural_gas_stream(self):
        return self.streams.natural_gas
    @property
    def natural_gas_GWP(self):
        return self.system.get_material_impact(self.natural_gas_stream, self.GWP)
        # return self.streams.natural_gas.F_mass*self.streams.natural_gas.characterization_factors[self.GWP]*self.annual_factor
    @property
    def other_material_GWP(self):
        annual_factor = self.annual_factor

        #membrane related impacts
        membrane_splitter = self.flowsheet.unit.S301
        membrane_maintenance_CF = GWP_factors['NF membrane maintenance']
        membrane_CF = GWP_factors['NF membrane']
        total_membrane_CF = membrane_CF + membrane_maintenance_CF
        # membrane_density = 1.14*1000 #kg/m3, density of polyamide
        membrane_thickness = 280 / 1e9  # nm, thickness of polyamide tfc membrane
        #F_vol is already adjusted based on lifetime of membrane in model parameter
        total_membrane_area = membrane_splitter.ins[0].F_vol/membrane_thickness #m2/hr
        membrane_GWP = total_membrane_area * total_membrane_CF * annual_factor

        total_other_materials_impact = (self.system.get_total_feeds_impact(self.GWP) - self.feedstock_GWP - self.natural_gas_GWP # ng GWP is included in steam GWP
                    + membrane_GWP)
        return total_other_materials_impact
    @property
    def EOL_GWP(self):
        return (sum([stream.get_atomic_flow("C") * self.chem.CO2.MW for stream in self.main_products])) * self.annual_factor

    #Purchased electricity-related impacts
    # Electricity is produced using steam and natural gas/or is directly purchased from the grid
    # Net electricity use includes net electricity requirement of all the units including the facilities (production + consumption)
    # + power = only consumption, - power = only production
    # The excess electricity automatically sold back to the grid is included in TEA
    # Since we set the system not to purchase electricity from the grid, excessive electricity is always produced (-power),
    # and thus impact of purchased electricity is always zero; electricity-related impacts are included in direct BT emissions GWP
    @property
    def net_electricity_use(self): # consumption + production
        #power_utility.rate = consumption + production
        #negative power_utility.rate = production
        net_electricity_use = (sum(i.power_utility.rate for i in self.units)) * self.annual_factor #kWh/yr
        return net_electricity_use
    @property
    def purchase_electricity(self):
        net_electricity = self.net_electricity_use
        if net_electricity < 1: # net elec =  consumption - production
            # electricity is produced and sold as product, so not included within input LCA calculation
            return 0
        else:  # electricity is consumed from the grid
            return net_electricity
    @property
    def purchased_electricity_GWP(self):
        return self.purchase_electricity * GWP_factors['electricity']
    @property
    def electricity_yield(self):
        net_electricity = self.net_electricity_use# electricity provided by the grid,'kW' * hr/yr = kWh/yr
        if net_electricity < -1:  # electricity is produced and sold as product
            return -net_electricity
        else:  # electricity is consumed from the grid
            return 0

    #Emission related impacts
    @property
    def emissions(self):#all emission streams
        extra_streams = [self.flowsheet.nf_membrane1, self.flowsheet.nf_membrane2]
        return [i for i in self.system.products if i not in self.main_products and not i in extra_streams and i not in self.flowsheet.ADP801.outs]
    @property
    def total_emissions_GWP(self):
        return (sum([stream.get_atomic_flow("C") for stream in self.emissions]) * self.chem.CO2.MW) * self.annual_factor
    @property
    def total_direct_emissions_GWP(self):
        return self.total_emissions_GWP + int(self.add_EOL_GWP)*self.EOL_GWP
    @property
    def total_direct_BT_emissions_GWP(self):
        return  (sum([i.get_atomic_flow('C') for i in self.BT.outs]) * self.chem.CO2.MW)*self.annual_factor\
         / self.total_emissions_GWP * self.total_direct_emissions_GWP
    @property
    def total_direct_nonBT_emissions_GWP(self):
        return self.total_direct_emissions_GWP - self.total_direct_BT_emissions_GWP

    #Steam-related impacts & breakdown
    @property
    def total_steam_GWP(self):
        return self.total_direct_BT_emissions_GWP + self.natural_gas_GWP
    #total steam GWP is divided into heating and electricity-related (cooling + non-cooling) impacts
    @property #amount of steam for heating
    def BT_steam_kJph_heating (self):
        return sum([i.duty for i in self.BT.steam_utilities])
    # amount of stean for producing electricity  in turbogenerator
    # (This electricity can be for both cooling and non-cooling purposes)
    @property
    def BT_steam_kJph_turbogen(self):
        KJpersec_to_KJhr = 3600 #kJ/hr = kW
        #BT.electricity_demand = sum of electricity demand of all other units in biorefinery [kW]
        #BT.power_utility.production = electricity produced by turbogenerator [kW]
        #if no excessive electricity is produced, BT.power_utility.production = net electricity_demand
        return KJpersec_to_KJhr * self.BT.power_utility.production / self.BT.turbogenerator_efficiency #not all steam is converted to electricity
    @property
    def BT_steam_kJph_total(self):
        return self.BT_steam_kJph_heating + self.BT_steam_kJph_turbogen
    @property
    def steam_fraction_heating(self):
        return self.BT_steam_kJph_heating / self.BT_steam_kJph_total
    @property
    def steam_fraction_turbogen(self):
        return self.BT_steam_kJph_turbogen / self.BT_steam_kJph_total
    #Steam used by turbogenerator to produce electricity is divided into cooling and non-cooling purposes
    #Total electricity demand can be calculated from system electricity consumption
    @property
    def net_electricity_demand(self):#Only consumption of the biorefinery
        #same as summing up power_utility_consumption of all unites in the biorefinery
        return self.system.get_electricity_consumption() # kWh/yr
    @property
    def cooling_electricity_demand(self):
        return (self.flowsheet.CT801.power_utility.rate + self.flowsheet.CWP801.power_utility.rate)* self.annual_factor
    @property
    def non_cooling_electricity_demand(self):#used by pumps, agitators, etc.
        return self.net_electricity_demand - self.cooling_electricity_demand
    @property
    def elec_fraction_cooling(self):
        return self.cooling_electricity_demand / self.net_electricity_demand
    @property
    def elec_fraction_non_cooling(self):
        return self.non_cooling_electricity_demand / self.net_electricity_demand
    @property
    def steam_fraction_elec_cooling(self):
        return self.elec_fraction_cooling * self.steam_fraction_turbogen
    @property
    def steam_fraction_elec_non_cooling(self):
        return self.elec_fraction_non_cooling * self.steam_fraction_turbogen
    @property
    def total_heating_GWP(self): #all heating demand is satisfied by steam
        return self.steam_fraction_heating*self.total_steam_GWP
    @property
    def steam_GWP_cooling_elec(self): #steam fraction GWP for cooling
        return self.steam_fraction_elec_cooling*self.total_steam_GWP
    @property
    def steam_GWP_non_cooling_elec(self): #steam fraction GWP for non-cooling
        return self.steam_fraction_elec_non_cooling*self.total_steam_GWP
    # @property
    # def steam_GWP_excessive_elec(self):
    #     return self.steam_fraction_turbogen_elec_yield* self.steam_fraction_turbogen*self.total_steam_GWP
    #Both purchased electricity (if any) + steam-generated electricity are used for cooling and non-cooling purposes
    #GWP of total cooling and non-cooling electricity is calculated
    #by summing up each of their fraction in the GWP of purchased electricity and steam-generated electricity
    @property #purchased electricity GWP for cooling = fraction of cooling electricity in total electricity * GWP of purchased electricity
    def total_cooling_GWP(self): #steam fraction GWP for cooling + purchased electricity GWP for cooling (if purchased)
        return self.steam_GWP_cooling_elec + self.purchased_electricity_GWP*self.elec_fraction_cooling
    @property
    def total_non_cooling_GWP(self): #steam fraction GWP for non-cooling + purchased electricity GWP for non-cooling (if purchased)
        return self.steam_GWP_non_cooling_elec + self.purchased_electricity_GWP*self.elec_fraction_non_cooling
    @property
    def total_electricity_demand_GWP(self):
        return self.total_cooling_GWP + self.total_non_cooling_GWP

    @property
    def hmf_purification_heating_demand_kJpH(self):
        hmf_purification_units = self.area_groups[2].units
        return sum([sum([i.duty for i in unit.heat_utilities if i.flow > 0 and i.duty > 0])
            for unit in hmf_purification_units])
    @property
    def hmf_purification_heating_GWP(self):
        return self.hmf_purification_heating_demand_kJpH/self.BT_steam_kJph_heating*self.total_heating_GWP
    @property
    def non_hmf_purification_heating_demand_kJpH(self):
        return self.BT_steam_kJph_heating-self.hmf_purification_heating_demand_kJpH
    @property
    def non_hmf_purification_heating_GWP(self):
        return self.total_heating_GWP-self.hmf_purification_heating_GWP
    @property
    def net_annual_GWP(self):
        return self.feedstock_GWP + self.other_material_GWP + self.natural_gas_GWP +self.total_direct_emissions_GWP + self.purchased_electricity_GWP
    @property
    def net_annual_GWP_alternative(self):
        return self.feedstock_GWP + self.other_material_GWP + self.total_direct_nonBT_emissions_GWP\
            + self.total_heating_GWP+ self.total_cooling_GWP+self.total_non_cooling_GWP
    #include steam GWP (=natural gas + direct BT emissions) and purchased_electricity
    @property
    def net_GWP_breakdown_alternative2(self):
        return [self.feedstock_GWP, self.other_material_GWP, self.total_direct_nonBT_emissions_GWP,
                self.hmf_purification_heating_GWP, self.non_hmf_purification_heating_GWP, self.total_cooling_GWP, self.total_non_cooling_GWP]
    @property
    def net_annual_GWP_alternative2(self):
        return sum(self.net_GWP_breakdown_alternative2)
    #steam GWP + purchased electricity GWP = total heating demand + total electricity demand + excessive electricity (sold as coproduct if any)
    @property
    def allocation_factor_economic(self):
        # sequence = hmf, biodiesel, furfural,crude glycerol
        products = self.main_products
        product_stream_names = self.main_product_names
        product_chemical_ID = self.main_product_chemical_ID
        part_ID = product_stream_names.copy()
        part_ID.append('Electricity') #name of all products
        part = [i.price * i.F_mass for i in products] #$/hr
        if self.electricity_yield > 1:  # electricity is sold as coproduct
            electricity_part = self.electricity_yield * bst.settings.electricity_price/self.annual_factor #$/hr
            part.append(electricity_part)
        else:
            part.append(0)
        # total: total sale
        total = sum(part)
        # allocation factor
        factors = [i / total for i in part]
        allocation_factors_economic = dict(zip(part_ID, factors))
        return allocation_factors_economic
    @property
    def allocation_factor_energy(self):
        # sequence = hmf, biodiesel, furfural,crude glycerol
        products = self.main_products
        product_stream_names = self.main_product_names
        product_chemical_ID = self.main_product_chemical_ID
        part_ID = product_stream_names.copy()
        part_ID.append('Electricity') #name of all products

        chemical_for_LHV = [self.chem[i] for i in product_chemical_ID]
        LHV_kJ_per_kg = [i.LHV / i.MW for i in chemical_for_LHV] #kJ/kg
        # glycerol is 80% of the total mass in crude glycerol stream
        dilution_factor = [1 ,1 ,1 ,0.8]
        mass_for_LHV = [i.F_mass * j for i, j in zip(products, dilution_factor)]
        part = [i * j for i, j in zip(mass_for_LHV, LHV_kJ_per_kg)] #kJ/hr
        if self.electricity_yield > 1:  # electricity is sold as coproduct
            electricity_part = self.electricity_yield / self.annual_factor * 60 #kJ/hr
            part.append(electricity_part)
        else:
            part.append(0)
        # total: total sale for economic, total heating value for energy
        total = sum(part)
        # allocation factor
        factors = [i / total for i in part]
        allocation_factors_energy = dict(zip(part_ID, factors))
        return allocation_factors_energy

    @property
    def HMF_GWP_economic(self):
        factor = self.allocation_factor_economic['HMF']
        return self.net_annual_GWP * factor/(self.streams.purified_hmf.F_mass*self.annual_factor)
    @property
    def HMF_GWP_energy(self):
        factor = self.allocation_factor_energy['HMF']
        return self.net_annual_GWP * factor/(self.streams.purified_hmf.F_mass*self.annual_factor)
    @property
    def Biodiesel_GWP_economic(self):
        factor = self.allocation_factor_economic['Biodiesel']
        return self.net_annual_GWP * factor/(self.streams.biodiesel.F_mass*self.annual_factor)
    @property
    def Biodiesel_GWP_energy(self):
        factor = self.allocation_factor_energy['Biodiesel']
        return self.net_annual_GWP * factor/(self.streams.biodiesel.F_mass*self.annual_factor)
    @property
    def HMF_GWP_displacement(self):
        displaced_products = self.main_products[1:] #HMF is not included
        dispalced_products_ID = ["biodiesel displacement", "furfural displacement", "glycerine displacement"]
        dilution_factor = [1, 1, 0.8]
        displaced_GWP = [i.F_mass * GWP_factors[j] * k * self.annual_factor for i, j, k in zip(displaced_products, dispalced_products_ID, dilution_factor)]
        if self.electricity_yield > 1:  # electricity is sold as coproduct
            electricity_gwp = self.electricity_yield * GWP_factors['electricity'] # kJ/hr
        else:
            electricity_gwp = 0
        displaced_GWP.append(electricity_gwp)
        net_GWP_after_displacement = self.net_annual_GWP - sum(displaced_GWP)
        HMF_GWP_displacement = net_GWP_after_displacement/(self.streams.purified_hmf.F_mass*self.annual_factor)
        return HMF_GWP_displacement

    def __repr__(self):
        return f'{type(self).__name__}({self.system.ID}, ...)'

    # def show(self):
    #     """Prints information on unit."""
    #     print(self._info())
    # _ipython_display_ = show

def create_lca_oilcane(sys,sys_chemcials,group_units,
                       main_feedstock,
                       products_for_allocation,
                       product_chemical_IDs,
                       product_names,
                       EOL_GWP=None,
                       lca_method=None):
    # product_streams = [purified_hmf, biodiesel, purified_furfural, crude_glycerol]
    # product_streams_chemical_IDs = ['HMF', 'Biodiesel', 'Furfural', 'Glycerol']
    # product_stream_names = ['HMF', 'Biodiesel', 'Furfural', 'Crude glycerol']  # for allocation factor & GWP calculation
    if EOL_GWP is None: EOL_GWP = False
    if lca_method is None: lca_method = OilcaneLCA
    lca = OilcaneLCA(
        system=sys,
        system_chemicals=sys_chemcials,
        area_groups = group_units,
        feedstock=main_feedstock,
        main_product=products_for_allocation,
        main_product_chemical_IDs=product_chemical_IDs,
        main_product_names=product_names,
        add_EOL_GWP=EOL_GWP,
    )
    return lca

# #To test the class, keep as comment:
# lca_baseline = create_lca_oilcane(baseline_sys,
#                           sys_chemcials=chem, group_units=area_groups,
#                           main_feedstock=feedstock,
#                           products_for_allocation=products,
#                           product_chemical_IDs = product_streams_chemical_IDs,
#                           product_names = product_stream_names,
#                           )
# oilcane_lca =  baseline_sys.LCA