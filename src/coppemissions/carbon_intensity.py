# -*- coding: utf-8 -*-

from .registry import ureg

import matplotlib.pyplot as plt


def carbon_intensity(field,total_power, percentage, plot=0):
    
    power = total_power*percentage
    
    if plot:
            plt.figure(); plt.plot(field.time_array,power); plt.title('Heating requirement power [kW]')
            
    return power.to('kW')