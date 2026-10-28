# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt

class Pump():
    """
    Class that represents a pump.
    """
    def __init__(self,identification,field,rate,suction_pressure,discharge_pressure,efficiency,save=1):    
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
        self.suction_pressure=suction_pressure      #[psi]
        self.discharge_pressure=discharge_pressure   #[psi]
        self.efficiency=efficiency                  #[ratio]
        self.rate=rate.to('oil_bbl/day')             #[bbl/d]
        
        if save:
            field.POWER_DICT[identification] = self.power()

    def bhp(self, plot=0):
        '''
        Return the brake horsepower of the pump [hp]
        '''
        bhp_array = (self.rate*(self.discharge_pressure-self.suction_pressure)/self.efficiency).to('hp')

        if plot:
            plt.figure(); plt.plot(self.field.time_array,bhp_array); plt.title('Pump BHP [hp]')

        return bhp_array #[hp]
        
    def power(self, plot=0):
        '''
        Return the consuming power of the pump [kWh]
        '''
        power_array = (self.rate*(self.discharge_pressure-self.suction_pressure)/self.efficiency).to('kW') #[kWh]

        if plot:
            plt.figure(); plt.plot(self.field.time_array,power_array); plt.title('Pump power [kW]')
        
        return power_array

