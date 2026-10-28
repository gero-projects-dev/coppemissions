# -*- coding: utf-8 -*-

from .registry import ureg
from .basic_classes import *
from .pump import *
from .compressor import *

import numpy as np



def gas_gathering(gathering_stream, fractional_loss_rate = 0.002554):
    """
    This function calculates the emissions from gas gathering.

    Parameters
    ----------
    gathering_stream : Stream object
        Post-vent stream to gas gathering.
    fractional_loss_rate : float
        Fractional Loss Rate [dimensionless]
    """

    GWP = {'CO2':1, 'CH4':25, 'C2H6':3.1, 'C3H8':3.1, 'CO':1.6, 'N2O':298}

    stream_mass_rate = gathering_stream.gas_mass_fraction_rate()

    fug_mass_rate = {}

    for i in stream_mass_rate.keys():
        fug_mass_rate[i] = fractional_loss_rate * stream_mass_rate[i]

    fugitive_emissions = np.zeros(len(gathering_stream.gas_mass_rate()))*ureg('metric_ton/day')

    for i in ['CO2', 'CH4', 'C2H6', 'C3H8', 'CO']:
        fugitive_emissions += fug_mass_rate[i] * GWP[i]

    
    outstream_mass_composition = {}

    for i in stream_mass_rate.keys():
        outstream_mass_composition[i] = (1-fractional_loss_rate) * gathering_stream.gas_mass_fraction()[i]

    molar = mass_fraction_to_molar(outstream_mass_composition)
    out_gas = Gas(molar)
    outstream = Stream('gas_gathering', gathering_stream.field, out_gas, gathering_stream.water_rate, gathering_stream.oil_rate, gathering_stream.gas_rate*(1-fractional_loss_rate), gathering_stream.temperature, gathering_stream.pressure)
    fugitive_emissions[np.isnan(fugitive_emissions)] = 0*ureg('metric_ton/day')

    return fugitive_emissions, outstream


def dehydration(field, inlet_stream, TEG_vol_rate=(689*ureg('gallons/day')).to('ft**3/day'), TEG_discharge_pressure = 786*ureg("psi"), TEG_eff=0.5, fractional_loss_rate = 0.000418):
    """
    This function calculates the emissions from triethylene glicol (TEG) gas dehydration.

    Parameters
    ----------
    inlet_stream : Stream object
        Inlet stream.
    TEG_vol_rate : Quantity object
        TEG volume rate [ft³/day] - Standard value = 689 gallons/day.
    TEG_discharge_pressure : Quantity object
        TEG discharge pressure [psi] - Standard value = 786 psi.
    TEG_eff : float
        TEG efficiency - Standard value = 0.5.
    fractional_loss_rate : array
        Fractional Loss Rate [dimensionless] - Standard value = 0.000418.
    """

    glycol_pump = Pump('glycol_pump',field, TEG_vol_rate*field.empty_field_array, 0, TEG_discharge_pressure*field.empty_field_array, TEG_eff*field.empty_field_array, save=0)
    electric_consumption = glycol_pump.power()*24

    GWP = {'CO2':1, 'CH4':25, 'C2H6':3.1, 'C3H8':3.1, 'CO':1.6, 'N2O':298}

    stream_mass_rate = inlet_stream.gas_mass_fraction_rate()

    fug_mass_rate = {}

    for i in stream_mass_rate.keys():
        fug_mass_rate[i] = fractional_loss_rate * stream_mass_rate[i]

    fugitive_emissions = np.zeros(len(inlet_stream.gas_mass_rate()))*ureg('metric_ton/day')

    for i in ['CO2', 'CH4', 'C2H6', 'C3H8', 'CO']:
        fugitive_emissions += fug_mass_rate[i] * GWP[i]

    
    outstream_mass_composition = {}

    for i in stream_mass_rate.keys():
        outstream_mass_composition[i] = (1-fractional_loss_rate) * inlet_stream.gas_mass_fraction()[i]

    molar = mass_fraction_to_molar(outstream_mass_composition)
    out_gas = Gas(molar)
    outstream = Stream('dehydration',inlet_stream.field, out_gas, inlet_stream.water_rate, inlet_stream.oil_rate, inlet_stream.gas_rate*(1-fractional_loss_rate), inlet_stream.temperature, inlet_stream.pressure)
    fugitive_emissions[np.isnan(fugitive_emissions)] = 0*ureg('metric_ton/day')

    return electric_consumption, fugitive_emissions, outstream


def amine_CO2_removal(field, inlet_stream, kvalue_amine_AGR = 0.000145 * ureg('gallons per minute')/ureg('ft**3/day'), contactor_pressure=500*ureg('psi'), eta_reboiler_AGR = 1.25, fractional_loss_rate = 0.000418):
    """
    This function calculates the emissions from amine CO2 removal.

    Parameters
    ----------
    inlet_stream : Stream object
        Inlet stream.
    kvalue_amine_AGR : Quantity object
        AGR amine K value [ft³/day] - Standard value = 145e-6 gallons/minute.
    contactor_pressure : Quantity object
        Pressure of the contactor [psi] - Standard value = 500 psi.  
    eta_reboiler_AGR : float 
        Reboiler Gas-fired heater fuel use [dimensionless] - Standard value = 1.25
    fractional_loss_rate : array
        Fractional Loss Rate [dimensionless]- Standard value = 0.000418
    """

    density_dict = {"N2":33.48, "O2":38.25, "CO2":52.61, "H2O":21.53, "CH4":19.17, "C2H6":35.94, "C3H8":52.70, "CO":33.48, "H2":2.41}

    reboiler_emission_factor = 0.59541e-7*ureg('metric_ton/BTU')


    GWP = {'CO2':1, 'CH4':25, 'C2H6':3.1, 'C3H8':3.1, 'CO':1.6, 'N2O':298}

    stream_mass_rate = inlet_stream.gas_mass_fraction_rate()
    stream_volumetric_rate = {}

    for i in ['N2', 'CO2', 'CH4', 'C2H6', 'C3H8', 'CO']:
        stream_volumetric_rate[i] = inlet_stream.gas_mass_fraction_rate()[i]/(density_dict[i] *(ureg.metric_ton/(1e6*ureg.ft**3)))



    amine_circulation_rate = kvalue_amine_AGR * stream_volumetric_rate['CO2']  

    circulator_pumps_hp = ((contactor_pressure+50*ureg('psi'))*amine_circulation_rate*0.00065).m*ureg('hp')
    boosters_reflux_pump_hp = 2 * amine_circulation_rate * 0.06 * ureg('hp')/ureg('gallons per minute')
    total_pumping_hp = circulator_pumps_hp + boosters_reflux_pump_hp
    c = 2967; r = -0.018
    total_pumping_kwh = total_pumping_hp * 24 * 2967*total_pumping_hp.m**r/3412.16 * ureg('kWh/(d * hp)')

    total_cooling_hp = 3* amine_circulation_rate * 0.36 * ureg('hp')/ureg('gallons per minute')
    total_cooling_kwh = total_cooling_hp * 24 * 2967*total_cooling_hp.m**r/3412.16 * ureg('kWh/(d * hp)')

    reboiler_heat_duty = (24*(amine_circulation_rate*72000)*1.15).m*ureg('BTU/d')
    reboiler_fuel_use = eta_reboiler_AGR*reboiler_heat_duty

    reboiler_emissions = reboiler_fuel_use*reboiler_emission_factor
    electric_consumption = total_cooling_kwh + total_pumping_kwh


    # FUGITIVES CALCULATION

    fug_mass_rate = {}

    for i in stream_mass_rate.keys():
        fug_mass_rate[i] = fractional_loss_rate * stream_mass_rate[i]

    fugitive_emissions = np.zeros(len(inlet_stream.gas_mass_rate()))*ureg('metric_ton/day')

    for i in ['CO2', 'CH4', 'C2H6', 'C3H8', 'CO']:
        fugitive_emissions += fug_mass_rate[i] * GWP[i]

    
    outstream_mass_composition = {}

    for i in stream_mass_rate.keys():
        outstream_mass_composition[i] = (1-fractional_loss_rate) * inlet_stream.gas_mass_fraction()[i]

    molar = mass_fraction_to_molar(outstream_mass_composition)
    out_gas = Gas(molar)
    outstream = Stream('amine_CO2_removal', inlet_stream.field, out_gas, inlet_stream.water_rate, inlet_stream.oil_rate, inlet_stream.gas_rate*(1-fractional_loss_rate), inlet_stream.temperature, inlet_stream.pressure)

    return reboiler_emissions, electric_consumption, fugitive_emissions, outstream


def chiller(inlet_stream, outlet_temperature = 274.817*ureg('K')):
    """
    This function calculates the energy consumption a gas chiller.
    Reference: National Energy Technology Laboratory (2012)

    Parameters
    ----------
    inlet_stream : Stream object
        Inlet stream.
    outlet_temperature : Quantity object
        Gas temperature exiting chiller [K] - Standard value = 35 F°
    """

    inlet_temperature = inlet_stream.temperature.to('rankine')
    outlet_temperature = outlet_temperature.to('rankine')    
    
    mass_flow = inlet_stream.gas_mass_rate()

    energy_consumption = 3.44*(mass_flow/6.111072)*(((inlet_temperature*5/9)-((outlet_temperature)*5/9))/56)*24*ureg('kilowatt_hour/(metric_ton*rankine)')
    inlet_stream.temperature = outlet_temperature.to('K')

    return energy_consumption.to('kW'), inlet_stream



def membrane_CO2_removal(field, inlet_stream, outlet_pressure = 500*ureg('psi'), eta_compressor_membrane = 0.75):
    """
    This function calculates the emissions from membrane CO2 removal.

    Parameters
    ----------
    inlet_stream : Stream object
        Inlet stream.
    inlet_temperature : Quantity object
        Gas temperature exiting chiller [K] - Standard value = 35 F°
    fractional_loss_rate : array
        Fractional Loss Rate [dimensionless]
    """

    density_dict = {"N2":33.48, "O2":38.25, "CO2":52.61, "H2O":21.53, "CH4":19.17, "C2H6":35.94, "C3H8":52.70, "CO":33.48, "H2":2.41}

    reboiler_emission_factor = 0.59541e-6*ureg('metric_ton/BTU')

    HC_sep_effectiveness = {"N2":0.58, 
                            "O2":1.0, 
                            "CO2":0.33, 
                            "H2O":1.0, 
                            "CH4":0.98, 
                            "C2H6":0.97, 
                            "C3H8":0.96,
                            "C4H10":0.94, 
                            "CO":1.0, 
                            "H2":1.0,
                            "H2S":0.49,
                            "SO2":1.0}
    
    CO2_sep_effectiveness = {"N2":0.42, 
                            "O2":0.0, 
                            "CO2":0.67, 
                            "H2O":0.0, 
                            "CH4":0.02, 
                            "C2H6":0.03, 
                            "C3H8":0.04,
                            "C4H10":0.06, 
                            "CO":0.0, 
                            "H2":0.0,
                            "H2S":0.51,
                            "SO2":0.0}
    
    HC_exit_P = inlet_stream.pressure - 150*ureg('psi')
    HC_exit_T = inlet_stream.temperature


    GWP = {'CO2':1, 'CH4':25, 'C2H6':3.1, 'C3H8':3.1, 'CO':1.6, 'N2O':298}

    stream_mass_rate = inlet_stream.gas_mass_fraction_rate()
    stream_volumetric_rate = {}

    for i in ['N2', 'CO2', 'CH4', 'C2H6', 'C3H8', 'CO']:
        stream_volumetric_rate[i] = inlet_stream.gas_mass_fraction_rate()[i]/(density_dict[i] *(ureg.metric_ton/(1e6*ureg.ft**3)))


    # FUGITIVES CALCULATION

    HC_mass_rate = {}
    CO2_mass_rate = {}

    for i in stream_mass_rate.keys():
        HC_mass_rate[i] = HC_sep_effectiveness[i] * stream_mass_rate[i]
        CO2_mass_rate[i] = CO2_sep_effectiveness[i] * stream_mass_rate[i]
    
    HC_out_mass_composition = {}
    CO2_out_mass_composition = {}

    for i in stream_mass_rate.keys():
        HC_out_mass_composition[i] = HC_sep_effectiveness[i] * inlet_stream.gas_mass_fraction()[i]
        CO2_out_mass_composition[i] = CO2_sep_effectiveness[i] * inlet_stream.gas_mass_fraction()[i]

    HCmolar = mass_fraction_to_molar(HC_out_mass_composition)
    HCout_gas = Gas(HCmolar)

    
    HCtotal_mass = sum(HC_mass_rate.values())
    
    HC_gas_rate = HCtotal_mass/HCout_gas.density()*inlet_stream.field.empty_field_array

    HCoutstream = Stream('HC_membrane_CO2_removal', inlet_stream.field, HCout_gas, inlet_stream.water_rate, inlet_stream.oil_rate, HC_gas_rate, ureg("296 K"), outlet_pressure)

    CO2molar = mass_fraction_to_molar(CO2_out_mass_composition)
    CO2out_gas = Gas(CO2molar)
    CO2outstream = Stream('CO2_membrane_CO2_removal',inlet_stream.field, CO2out_gas, inlet_stream.water_rate, inlet_stream.oil_rate, inlet_stream.gas_rate-HC_gas_rate, ureg("296 K"), HC_exit_P)


    compression_ratio = outlet_pressure / HC_exit_P

    Comp = Compressor('membrane_CO2_removal', field)
    c_compressor = Comp.get_compressor_energy_consumption(eta_compressor_membrane, compression_ratio,inlet_stream.gas_rate, HC_exit_T, HC_exit_P, inlet_stream.gas, plot=0)[0]

    return c_compressor, HCoutstream, CO2outstream



class GasProcessing():

    standard_parameters = {
        "gas_gathering_fractional_loss_rate":0.002554,
        "dehydration_TEG_vol_rate":(689*ureg('gallons/day')).to('ft**3/day'),
        "dehydration_TEG_discharge_pressure":786*ureg("psi"),
        "dehydration_TEG_efficiency":0.5,
        "dehydration_fractional_loss_rate":0.000418,
        "membrane_chiller_outlet_temperature":274.817*ureg('K'),
        "membrane_outlet_pressure": 500*ureg('psi'), 
        "membrane_eta_compressor": 0.75,
        "amine_kvalue_amine_AGR": 0.000145 * ureg('gallons per minute')/ureg('ft**3/day'), 
        "amine_contactor_pressure": 500*ureg('psi'), 
        "amine_eta_reboiler_AGR": 1.25, 
        "amine_fractional_loss_rate": 0.000418
    }

    def __init__(self, field, inlet_stream):
        self.standard_parameters = field.data
        self.field = field
        self.inlet_stream = inlet_stream
        self.arq = field.data['gas_processing_arq']
        
        if "Gas Gathering" not in self.arq:
            raise Exception("Gas Gathering is mandatory for the Gas Processing Architecture")
        if ("Amine" in self.arq) and ("Membrane" in self.arq):
            raise Exception("The Gas Processing Architecture can't include the Amine and Membrane separation methods simultaneously") 
        


        fugitive_emissions_dict = {}
        energy_consumption_dict = {}
        other_emissions_dict = {}
        fugitive_emissions = 0*self.field.empty_field_array*ureg('metric_ton/day')
        energy_consumption = 0*self.field.empty_field_array*ureg('kW')
        other_emissions = 0*self.field.empty_field_array*ureg('metric_ton/day')

        GG = gas_gathering(self.inlet_stream, self.standard_parameters["gas_gathering_fractional_loss_rate"])
        energy_consumption_dict['gas_gathering'] = GG[0]; fugitive_emissions += GG[0]
        outstream = GG[1]

        #print('gas gatheing'); print(outstream.gas_mass_fraction_rate())

        if "Dehydration" in self.arq:
            DEH = dehydration(self.field, outstream, self.standard_parameters["dehydration_TEG_vol_rate"],  TEG_discharge_pressure = self.standard_parameters["dehydration_TEG_discharge_pressure"], TEG_eff=self.standard_parameters["dehydration_TEG_efficiency"])
            energy_consumption_dict['dehydration'] = DEH[0]; energy_consumption += DEH[0]
            fugitive_emissions_dict['dehydration'] = DEH[1]; fugitive_emissions += DEH[1]
            outstream = DEH[2]

            #print('deh'); print(outstream.gas_mass_fraction_rate())
        
        if "Membrane" in self.arq:
            CHL = chiller(outstream, self.standard_parameters["membrane_chiller_outlet_temperature"])
            energy_consumption_dict['chiller'] = CHL[0]; energy_consumption += CHL[0]
            outstream = CHL[1]

            #print('chiller');print(outstream.gas_mass_fraction_rate())

            MEM = membrane_CO2_removal(self.field, outstream, outlet_pressure = self.standard_parameters["membrane_outlet_pressure"], eta_compressor_membrane = self.standard_parameters["membrane_eta_compressor"])
            energy_consumption_dict["membrane"] = MEM[0]; energy_consumption += MEM[0]
            outstream = MEM[1]
            CO2_outstream = MEM[2]

            #print('mem'); print(outstream.gas_mass_fraction_rate())

        if "Amine" in self.arq:
            AMN = amine_CO2_removal(self.field, outstream, self.standard_parameters["amine_kvalue_amine_AGR"], self.standard_parameters["amine_contactor_pressure"], self.standard_parameters["amine_eta_reboiler_AGR"], self.standard_parameters["amine_fractional_loss_rate"])
            other_emissions_dict['amine'] = AMN[0]; other_emissions += AMN[0]
            energy_consumption_dict['amine'] = AMN[1]; energy_consumption += AMN[1]
            fugitive_emissions_dict['amine'] = AMN[2]; fugitive_emissions += AMN[2]
            outstream = AMN[3]
            CO2_outstream = Stream('Gas_Processing_Outstream', outstream.field, outstream.gas, 0*outstream.water_rate, 0*outstream.oil_rate, 0*outstream.gas_rate, outstream.temperature, outstream.pressure)

        self.outstream = outstream
        self.CO2_outstream = CO2_outstream
        self.fugitive_emissions = fugitive_emissions
        self.energy_consumption = energy_consumption
        self.other_emissions = other_emissions
        self.fugitive_emissions_dict = fugitive_emissions_dict
        self.energy_consumption_dict = energy_consumption_dict
        self.other_emissions_dict = other_emissions_dict
        
        
        self.field.EMISSIONS['Gas Processing fugitives'] = self.fugitive_emissions + self.other_emissions
        self.field.POWER_DICT['Gas Processing'] = self.energy_consumption
        
        