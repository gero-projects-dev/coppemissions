# -*- coding: utf-8 -*-

from .registry import ureg
import numpy as np

ngturbine_data = {
    'ngt_power_nominal': 32.94,           # Capacity (MW)
    'trl': 9,                             # Technology Readiness Level (1-9)
    'intermittence': 0,                   # Power Supply Intermittence (0-infinity)
    'ngt_type': '?'              # Type of Turbine (ex. Aeroderivative 15°C 60% RH)
    }

def ngt_curveload_emission_factor(ngt_type, load):
    x = load*100
    
    if ngt_type == 'AERO-25-30':
        y = 834815.6746*(x**0.25)+117401.739*(x**0.5)-763094.1303*(x**(1/3))-666.0073938*(x**1)-204868.4436
    if ngt_type == 'AERO-35-30':
        y = -0.001600320705*(x**3)+0.4110311468*(x**2)-36.49971768*(x**1)+1656.117979
    if ngt_type == 'AERO-15-60':
        y = 841158.504*(x**0.25)+118462.7765*(x**0.5)-769214.412*(x**(1/3))-674.6007922*(x**1)-206281.649
    if ngt_type == 'W31SG':
        y = 2095.94633 + 14.71266785*x -0.05290732899*(x**2)-562.1742963*np.log(x)
    if ngt_type == '?':
       y = -171892.2691*(x**(1/8))+86967.85478*(x**(1/4))-7424.013243*(x**(1/2))+64.92909963*(x**1)+99020.45587
        
    return y*1E-3*ureg('metric_ton/(MW*h)') #metric_ton/MWh


def ngt_emissions(demand_power,ngturbine_data=ngturbine_data):
    
    nominal = ngturbine_data['ngt_power_nominal']*ureg('MW')                                    # Nominal power of each turbine in MWh
    ngt_number_on = np.floor((demand_power/nominal).m)                                              # number of turbines in total load to feed the demand power
    total_load_power = ngt_number_on*nominal
    
    partial_load_power = demand_power - ngt_number_on*nominal                                   # power from partial load turbine to feed demand power in MWh
    load_ratio = partial_load_power / nominal                                                   # load ratio of partial load turbine per day
    
    part_load_factor = ngt_curveload_emission_factor(ngturbine_data['ngt_type'], load_ratio)    # emission factor of partial load turbine per day
    ngt_part_load_emissions = part_load_factor*partial_load_power                               # total emissions from partial load turbine in metric_ton
    
    total_load = np.ones(len(demand_power))
    ngt_total_load_factor = ngt_curveload_emission_factor(ngturbine_data['ngt_type'], total_load)
    ngt_total_load_emissions = ngt_total_load_factor*total_load_power                           # total emissions from total load turbines in metric_ton
    
    total_emissions_year = ngt_part_load_emissions+ngt_total_load_emissions
    
    return total_emissions_year.to('metric_ton/day')



def ngt_emissions2(demand_power,ngturbine_data=ngturbine_data):
    
    nominal = ngturbine_data['ngt_power_nominal']*ureg('MW')                                    # Nominal power of each turbine in MWh
    ngt_number_on = np.ceil(demand_power.to('kW')/nominal.to('kW'))                             # number of turbines to feed the demand power
    power_each = demand_power/ngt_number_on                                                     # Power of each turbine to feed the demand power
    load_ratio = (power_each / nominal).to('dimensionless').magnitude                           # Load ratio of each turbine to feed the demand power
    emission_factor = ngt_curveload_emission_factor(ngturbine_data['ngt_type'], load_ratio)     # emission factor of each turbine [metric_ton/MWh]

    emissions_each = emission_factor*power_each                                                 # total emissions from partial load turbine in metric_ton

    total_emissions = emissions_each*ngt_number_on
    total_emissions = np.nan_to_num(total_emissions)
    
    return total_emissions.to('metric_ton/day'), ngt_number_on

