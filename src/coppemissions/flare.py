# -*- coding: utf-8 -*-

from .registry import ureg
from .basic_classes import *

import copy
import matplotlib.pyplot as plt


class Flare():
    """
    Class used to represent flare emission sources.
    """
    
    def __init__(self, field,stream,volume_total=[0],FOR=0,volume_tanks=0,volume_mine=0):
        """
        Parameters
        ----------
        gas : fluid object
            Fluid object that is burned in the flare.
        volume_total : float
            Total volume of gas that is burned in the flare [ft³/d]. If not inserted, the next parameters must be inserted.
        FOR : float
            Flaring-to-oil ratio [ft³/bbl of oil].
        volume_tanks : float
            Volume of gas routed from oil storage tanks to flares [ft³/d].
        volume_mine : float
            Volume of gas routed to the flares from mine off-gas collection [ft³/d].
        """
        self.efficiency = (1-0.03)*0.9611 #IMPLEMENTAR REF: OPGEE documentation
        self.stream = stream
        self.gas = copy.copy(self.stream.gas)
        self.FOR = FOR
        self.field = field
        self.volume_tanks = volume_tanks
        self.volume_mine = volume_mine
        self.volume_total = volume_total
        if len(volume_total)==1:
            if volume_total == [0]:
                if FOR==0:
                    raise ValueError("If the total flared volume isn't inserted, the FOR is needed.")
                self.volume = (self.FOR*field.oil_prod + self.volume_tanks + self.volume_mine) #Volume of flared gas [ft³/d]
            else:
                self.volume = self.field.empty_field_array*volume_total[0]
        else:
            self.volume = copy.copy(volume_total)
        self.fraction_flared = self.volume/self.stream.gas_rate
        
        field.EMISSIONS['Flaring'] = self.emissions_eng()
        
        

    # def emissions_EF(self, plot=0):
    #     """
    #     This method returns the flare's emission calculated by an emission factor method [metric_tonCO2/d]
    #     """
    #     EM = emission_factors['flare'] * (ureg.g/(ureg.ft**3))
    #     if plot:
    #         plt.figure(); plt.plot(self.field.time_array,(EM*self.volume).to('metric_ton/day')); plt.title('Flare - Emissions calculated by Emission Factor Method')
    #     return (EM*self.volume).to('metric_ton/day') #[metric_tonCO2/d]
    
    def emissions_eng(self, plot=0):
        """
        This method returns the flare's emission calculated by an engineering method [metric_tonCO2/d]
        """
        Emission_noncombustion = 0
        Emission_combustion = 0

        GWP = {'CO2':1, 'CH4':25, 'C2H6':3.1, 'C3H8':3.1, 'CO':1.6, 'N2O':298}

        #print('/nVolume: '+str(self.volume)+'/n')
        
        for i in list(set(list(self.gas.composition.keys())) & set(["CO2","CH4","C2H6","C3H8"])):
            
            Emission_noncombustion += self.volume*(1-self.efficiency)*self.gas.composition[i]*density(i)*GWP[i]
            Emission_combustion += self.volume*self.efficiency*self.gas.composition[i]*density(i)*tau(i)

            #print('\n'+str(i)+' Comb: '+ str(self.volume*self.efficiency*self.gas.composition[i]*density(i)*tau(i))+'\n')
            #print('\n'+str(i)+' N_Comb: '+ str(self.volume*(1-self.efficiency)*self.gas.composition[i]*density(i)*gwp(i))+'\n')

        if plot:
            plt.figure(); plt.plot(self.field.time_array,(Emission_noncombustion +  Emission_combustion).to('metric_ton/day')); plt.title('Flare - Emissions calculated by Engineering Method')

        return (Emission_noncombustion +  Emission_combustion).to('metric_ton/day') #[metric_tonCO2/d]
    

    def methane_slip(self):
        pass
    
    def gas_stream_to_vent(self):
        outstream = copy.copy(self.stream)
        #print(self.fraction_flared)
        outstream.gas_rate = outstream.gas_rate * (1.0-self.fraction_flared)
        return outstream