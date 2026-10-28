# -*- coding: utf-8 -*-

from .registry import ureg

import numpy as np
from .basic_classes import *

class Vent():

    def __init__(self, identification, field, vent_stream):
        """
        This function calculates the emissions from purposeful venting.
    
        Parameters
        ----------
        vent_tream : Stream object
            Post-flare stream to vent.
        surface_piping_leakage : array
            Surface Piping CH4 Leakage Array [metric_ton CH4 / day]
        fractional_loss_rate : float
            Fractional Loss Rate [dimensionless]
        """
        
        surface_piping_leakage = field.data["venting_surface_piping_leakage"]
        fractional_loss_rate = field.data["venting_fraction"]
    
        if surface_piping_leakage.m == 0.0:
            surface_piping_leakage = np.zeros(len(vent_stream.gas_mass_rate()))*ureg('metric_ton/day')
    
        if len(vent_stream.gas_mass_rate()) != len(surface_piping_leakage):
            raise ValueError("The Stream gas rate and surface piping leakage arrays must have the same length.")
        
        
        GWP = {'CO2':1, 'CH4':25, 'C2H6':3.1, 'C3H8':3.1, 'CO':1.6, 'N2O':298}
    
        stream_mass_rate = vent_stream.gas_mass_fraction_rate()
    
        vent_mass_rate = {}
    
        for i in stream_mass_rate.keys():
            if i not in ['CH4', 'N2', 'CO2']:
                vent_mass_rate[i] = 0.0*stream_mass_rate[i]
            else:
                vent_mass_rate[i] = fractional_loss_rate * stream_mass_rate[i]
    
        misc_mass_rate = {}
        misc_mass_rate['CH4'] = surface_piping_leakage
    
        for i in stream_mass_rate.keys():
            misc_mass_rate[i] = (misc_mass_rate['CH4'] / ((1-fractional_loss_rate) * stream_mass_rate['CH4'])) * stream_mass_rate[i]
    
        venting_emissions = np.zeros(len(vent_stream.gas_mass_rate()))*ureg('metric_ton/day')
    
        for i in ['CO2', 'CH4', 'C2H6', 'C3H8', 'CO']:
            venting_emissions += (vent_mass_rate[i] + misc_mass_rate[i]) * GWP[i]
    
        outstream_mass_composition = {}
    
        for i in stream_mass_rate.keys():
            outstream_mass_composition[i] = (1-fractional_loss_rate) * vent_stream.gas_mass_fraction()[i]
    
        molar = mass_fraction_to_molar(outstream_mass_composition)
        out_gas = Gas(molar)
        outstream = Stream(identification+"_outstream", vent_stream.field, out_gas, vent_stream.water_rate, vent_stream.oil_rate, vent_stream.gas_rate*(1-fractional_loss_rate), vent_stream.temperature, vent_stream.pressure)
    
        venting_emissions[np.isnan(venting_emissions)] = 0*ureg('metric_ton/day')
        
        self.emissions = venting_emissions
        field.EMISSIONS[identification] = venting_emissions
        
        self.outstream = outstream


