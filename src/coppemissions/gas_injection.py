# -*- coding: utf-8 -*-

from .registry import ureg
from .compressor import *
import matplotlib.pyplot as plt


def Gas_Injection(stream, field, res_press, efficiency, plot=0):
    
    discharge_pressure = res_press + ureg('1 atm').to('psi')+ ureg('500 psi')
    c = Compressor()
    inlet_pressure = stream.pressure
    ratio = discharge_pressure/inlet_pressure

    power,T_out,P_out = c.get_compressor_energy_consumption(efficiency,ratio, stream.gas_rate, stream.temperature, inlet_pressure, stream.gas, plot=0)

    if plot:
            plt.figure(); plt.plot(field.time_array,power); plt.title('Gas injection power [kW]')
            
    return power.to('kW'),T_out,P_out