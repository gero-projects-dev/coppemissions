# -*- coding: utf-8 -*-

from .registry import ureg
from .basic_classes import *


class Crude_Oil_Storage:
    
    def __init__(self, identification, field, inlet_stream):
        
        loss_rate=field.data['storage_loss_rate']
        storage_gas_recovery_fraction = field.data['storage_gas_recovery_fraction']
        ch4_fraction_in_storage_gas = field.data['ch4_fraction_in_storage_gas']
        storage_gas_flare_fraction = field.data['storage_gas_flare_fraction']
        
        GWP = {'CO2':1, 'CH4':25, 'C2H6':3.1, 'C3H8':3.1, 'CO':1.6, 'N2O':298}
        storage_gas_composition = {"N2":0, "O2":0, "CO2":0.11, "H2O":0, "CH4":0.55, "C2H6":0.23, "C3H8":0.13, "CO":0, "H2":0}
        #storage_gas_recovery_fraction = 0.06
        #storage_gas_flare_fraction = 0.43
        #ch4_fraction_in_storage_gas = 0.55
        #loss_rate = 0.945 * ureg('kg/oil_bbl')
        
        oil_volume_storage = inlet_stream.oil_rate

        flashing_gas_mass = (loss_rate*oil_volume_storage/ch4_fraction_in_storage_gas).to('metric_ton/day')

        fugitive_emissions_from_tank = {}

        for i in list(storage_gas_composition.keys()):
            fugitive_emissions_from_tank[i] = (1-storage_gas_recovery_fraction-storage_gas_flare_fraction)*flashing_gas_mass*storage_gas_composition[i]

        storage_fugitives = 0

        for i in ['CO2', 'CH4', 'C2H6', 'C3H8', 'CO']:
            storage_fugitives += fugitive_emissions_from_tank[i]*GWP[i]
        
        recovered_mass =  storage_gas_recovery_fraction*flashing_gas_mass
        flared_mass = storage_gas_flare_fraction*flashing_gas_mass
        
        field.EMISSIONS['Flash_Fugitives'] = storage_fugitives
        
        flash_recovered_flow = recovered_mass/inlet_stream.gas.density()
        flash_flared_flow = flared_mass/inlet_stream.gas.density()
        
        self.recovered_stream = Stream(identification+"_recovered",field,inlet_stream.gas,inlet_stream.water_rate,inlet_stream.oil_rate, flash_recovered_flow,inlet_stream.temperature,field.data['flash_gas_initial_pressure'])
        self.flared_stream = Stream(identification+"_flared",field,inlet_stream.gas,inlet_stream.water_rate,inlet_stream.oil_rate,flash_flared_flow,inlet_stream.temperature,field.data['flash_gas_initial_pressure'])
        
        
