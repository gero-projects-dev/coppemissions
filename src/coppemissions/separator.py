# -*- coding: utf-8 -*-

from .registry import ureg
from .basic_classes import *
from .compressor import *

import math
import numpy as np
import matplotlib.pyplot as plt
import copy


class Separator():
    """
    Class that represents a triphasic separator.
    """

    stand_dict = {
    'inlet_pressure_criterion':[1000*ureg('psi'), 500*ureg('psi')],
    'stage_pressures':[500*ureg('psi'), 250*ureg('psi'), 100*ureg('psi')],
    'liquid_separator_pressure': 100*ureg('psi'),
    'gas_separator_temperature': 300*ureg('K'),
    'separator_outlet_pressure': 500*ureg('psi'),
    'fractional_loss_rate': 0.0002
    }

    def __init__(self, field, stream):    
        """
        Parameters
        ----------
        rate : array (same length as time_array)
            Pump discharge rate [bbl/d]
        suction_pressure : Quantity object or array (same length as time_array)
            Pump suction pressure [psi]
        discharge_pressure : Quantity object or array (same length as time_array)
            Pump discharge pressure [psi]
        efficiency : float or array (same length as time_array)
            Pump efficiency [ratio].
        """
        self.field = field
        dict = field.data['separator_dict']
        
        if self.field.wellhead_p >= dict['inlet_pressure_criterion'][0]:
            self.separator_stages = 3
        elif self.field.wellhead_p >= dict['inlet_pressure_criterion'][1]:
            self.separator_stages = 2
        else: 
            self.separator_stages = 1

        self.stream = stream
        self.liquid_separator_pressure = dict['liquid_separator_pressure']
        self.gas_separator_temperature = dict['gas_separator_temperature']
        self.outlet_pressure = dict['separator_outlet_pressure']
        self.separator_stages_pressures = dict['stage_pressures']
        self.fractional_loss_rate = dict['fractional_loss_rate']
        self.compressor_efficiency = dict['separation_compressor_efficiency']
        
        #CHOKE
        
        self.stream.pressure = dict['stage_pressures'][0]
        
        field.POWER_DICT['Separator'] = self.calculate_consumption()[0]
        field.EMISSIONS['Separator fugitives'] = self.fugitives()[0]


    def solution_gor(self, T, P, oil_sg, gas_sg ):
        a1 = 5.527215
        a2 = 0.783716
        a3 = 1.841408
        
        R1 = P.to('psi')**(1/a2)*oil_sg**(-a1/a2)
        R2 = (math.e**((a3/a2)*gas_sg*oil_sg))
        R3 = (gas_sg * T.to('rankine'))
        return (R1*R2/R3).m*ureg('ft³/oil_bbl')

    def calculate_consumption(self,plot=0):

        T = self.stream.temperature
        P = self.stream.pressure
        c = Compressor('Separator', self.field)
        
        gor1 = self.solution_gor(T, P, self.stream.oil.specific_gravity(), self.stream.gas.specific_gravity)
        removed_free_gas1 = self.field.GOR - gor1

        rate = removed_free_gas1*self.stream.oil_rate

        compression_ratio1 = self.outlet_pressure / P
        
        C1, T1, P1 = c.get_compressor_energy_consumption(self.compressor_efficiency,compression_ratio1,rate, T, P, self.stream.gas, plot=0)
        
        C2, T2, P2, C3, T3, P3 = 0,0,0,0,0,0


        if self.separator_stages > 1:

            T = (T + self.gas_separator_temperature)/2
            P = self.separator_stages_pressures[1]

            gor2 = self.solution_gor(T, P, self.stream.oil.specific_gravity(), self.stream.gas.specific_gravity)
            removed_free_gas2 = gor1-gor2

            rate = removed_free_gas2*self.stream.oil_rate

            compression_ratio2 = self.outlet_pressure / P
            

            C2, T2, P2 = c.get_compressor_energy_consumption(self.compressor_efficiency,compression_ratio2,rate, T, P, self.stream.gas, plot=0)
            C3, T3, P3 = 0,0,0

            
        if self.separator_stages > 2:

            T = self.gas_separator_temperature
            P = self.separator_stages_pressures[2]

            gor3 = self.solution_gor(T, P, self.stream.oil.specific_gravity(), self.stream.gas.specific_gravity)
            removed_free_gas3 = gor2 - gor3
            
            rate = removed_free_gas3*self.stream.oil_rate

            compression_ratio3 = self.outlet_pressure / P
            
            C3, T3, P3 = c.get_compressor_energy_consumption(self.compressor_efficiency,compression_ratio3,rate, T, P, self.stream.gas, plot=0)

        if plot:
            plt.figure(); plt.plot(time_array,np.array(C1+C2+C3)); plt.title('Separator power consumption [kW]'); plt.grid()

        return (C1+C2+C3).to('kW'), [C1, C2, C3]
    
    def gas_outstream(self):
        #outstream = copy.copy(self.stream)
        outstream = Stream('SEP_out_gas', self.field, self.stream.gas, 0,0,0,0,0)
        outstream.water_rate = 0*self.field.empty_field_array*ureg(('oil_bbl/d'))
        outstream.oil_rate = 0*self.field.empty_field_array*ureg(('oil_bbl/d'))
        outstream.gas_rate = self.stream.gas_rate
        outstream.temperature = self.gas_separator_temperature
        outstream.pressure = self.outlet_pressure
        return outstream
    
    def water_outstream(self):
        #outstream = copy.copy(self.stream)
        outstream = Stream('SEP_out_wtr', self.field, self.stream.gas, 0,0,0,0,0)
        outstream.water_rate = self.stream.water_rate
        outstream.gas_rate = 0*self.field.empty_field_array*ureg("ft**3/d")
        outstream.oil_rate = 0*self.field.empty_field_array*ureg(('oil_bbl/d'))
        outstream.temperature = self.gas_separator_temperature
        outstream.pressure = self.liquid_separator_pressure
        return outstream
    
    def oil_outstream(self):
        #outstream = copy.copy(self.stream)
        outstream = Stream('SEP_out_oil', self.field, self.stream.gas, 0,0,0,0,0)
        outstream.water_rate = 0*self.field.empty_field_array*ureg(('oil_bbl/d'))
        outstream.gas_rate = 0*self.field.empty_field_array*ureg("ft**3/d")
        outstream.oil_rate = self.stream.oil_rate
        outstream.temperature = self.gas_separator_temperature
        outstream.pressure = self.liquid_separator_pressure
        return outstream
    
    def fugitives(self):

        fug_water = 0*self.field.empty_field_array*ureg(('oil_bbl/d'))
        fug_oil = 0*self.field.empty_field_array*ureg(('oil_bbl/d'))
        

        GWP = {'CO2':1, 'CH4':25, 'C2H6':3.1, 'C3H8':3.1, 'CO':1.6, 'N2O':298}

        stream_mass_rate = self.stream.gas_mass_fraction_rate()
        fug_mass_rate = {}

        for i in stream_mass_rate.keys():
            fug_mass_rate[i] = self.fractional_loss_rate * stream_mass_rate[i]

        fugitive_emissions = np.zeros(len(self.stream.gas_mass_rate()))*ureg('metric_ton/day')

        for i in ['CO2', 'CH4', 'C2H6', 'C3H8', 'CO']:
            fugitive_emissions += fug_mass_rate[i] * GWP[i]

        outstream_mass_composition = {}

        for i in stream_mass_rate.keys():
            outstream_mass_composition[i] = (self.fractional_loss_rate) * self.stream.gas_mass_fraction()[i]

        molar = mass_fraction_to_molar(outstream_mass_composition)
        out_gas = Gas(molar)
        fug_stream = Stream('Sep_fug', self.stream.field, out_gas, fug_water, fug_oil, self.stream.gas_rate*(self.fractional_loss_rate), self.gas_separator_temperature , self.outlet_pressure)

        return fugitive_emissions, fug_stream
    
