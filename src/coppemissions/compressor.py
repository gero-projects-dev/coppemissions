# -*- coding: utf-8 -*-

from .registry import ureg

import math
import numpy as np
import gascompressibility as gascomp
from gascompressibility.pseudocritical import Sutton
import matplotlib.pyplot as plt

_power = [1, 1 / 2, 1 / 3, 1 / 4, 1 / 5]

class Compressor():
    """
    Class that represents a Compressor.
    """
    
    def __init__(self,identification, field, save=1):
        self.identification = identification
        self.field = field
        self.save = save
        
    def Z_factor(self,reduced_temperature, reduced_pressure):
        """
        Calculate the compressibility factor (Z) for a given reduced temperature and reduced pressure
        using the Redlich-Kwong equation of state. The Redlich-Kwong equation of state is a cubic equation.
        The function calculates the three possible roots (real or complex) and returns the highest real root as
        the compressibility factor. The input quantities are first converted to dimensionless fractions.

        Args:
            reduced_temperature (pint.Quantity): Reduced temperature as a Pint Quantity object, defined as T/Tc (ratio of temperature to critical temperature).
            reduced_pressure (pint.Quantity): Reduced pressure as a Pint Quantity object, defined as P/Pc (ratio of pressure to critical pressure).

        Returns:
            pint.Quantity: Compressibility factor (Z) as a Pint Quantity object (dimensionless fraction) for the given reduced temperature and reduced pressure.

        """
        tr = reduced_temperature.magnitude
        pr = reduced_pressure.magnitude

        a = 0.42748 * (pr / tr ** 2.5)
        b = 0.08664 * (pr / tr)

        alpha = (1 / 3) * (3 * (a - b - b ** 2) - 1)
        beta = (1 / 27) * (-2 + (9 * (a - b - b ** 2)) - (27 * a * b))
        d = (beta ** 2 / 4) + (alpha ** 3 / 27)

        if d < 0:
            theta = math.acos(-np.sign(beta) * (math.sqrt((beta ** 2 / 4) / (-alpha ** 3 / 27))))
            z_roots = [2 * math.sqrt(- alpha / 3) * math.cos((theta / 3) + (i * ((math.pi * 2) / 3))) + (1 / 3) for i in
                       range(3)]
        else:
            a_star = np.cbrt((-beta / 2) + np.sqrt(d))
            b_star = np.cbrt((-beta / 2) - np.sqrt(d))
            z_roots = [a_star + b_star + 1 / 3] if d > 0 else [a_star + b_star + 1 / 3] + [
                -(1 / 2) * (a_star + b_star) + (1 / 3) * i for i in [2, 3]]

        z = max(z_roots)

        return z*ureg('')

    def get_compressor_work_temp(self, inlet_temp, inlet_press, fluid, compression_ratio, num_of_compression):
        """

        :param field:
        :param inlet_temp:
        :param inlet_press:
        :param gas_stream:
        :param compression_ratio:
        :param num_of_compression:
        :return:(float) overall work from compressor which has maximum five stages (unit = hp*day/scf)
        """

        gas = fluid
        #corrected_temp = gas.corrected_pseudocritical_temperature(gas_stream)
        #corrected_press = gas.corrected_pseudocritical_pressure(gas_stream)
        #ratio_of_specific_heat = gas.ratio_of_specific_heat(gas_stream)

        z_factor = gascomp.calc_z(sg=gas.mix.mw()[0]/28.967, T=inlet_temp.to('degF').magnitude, P=inlet_press.to('psi').magnitude, CO2=gas.composition['CO2'], H2S=0, N2=gas.composition['N2'])
        ratio_of_specific_heat = gas.mix.gam(T=inlet_temp.to('kelvin'))

        work_temp1 = 3.027 * 14.7 / (60 + 460) * ratio_of_specific_heat / (ratio_of_specific_heat - 1)
        ratio = (ratio_of_specific_heat - 1) / ratio_of_specific_heat
        work = 0.0
        for j in range(num_of_compression):
            inlet_reduced_temp = Sutton().calc_Tr(sg=gas.mix.mw()[0]/28.967, T=inlet_temp.to('degF').magnitude, CO2 = gas.composition['CO2'],  H2S=0)
            inlet_reduced_press = Sutton().calc_Tr(sg=gas.mix.mw()[0]/28.967, T=inlet_press.to('psi').magnitude, CO2 = gas.composition['CO2'],  H2S=0)
            z_factor = self.Z_factor(inlet_reduced_temp*ureg(''), inlet_reduced_press*ureg(''))

            work_temp2 = (compression_ratio ** z_factor) ** ratio - 1
            work += work_temp1 * work_temp2 * inlet_temp.to("rankine")

            delta_temp = (inlet_temp.to("rankine") * compression_ratio ** (z_factor * ratio) - inlet_temp) * 0.2
            inlet_temp = ureg.Quantity(inlet_temp.m + delta_temp.m, "degF")
            inlet_press = compression_ratio * inlet_press

        work_sum = ureg.Quantity(work.m/1e6, "hp*day/ft³")
        return work_sum, inlet_temp, inlet_press

    @staticmethod
    def get_compression_ratio_stages(overall_compression_ratio_stages):
        max_stages = len(_power)
        compression_ratio_per_stages = []

        for compression_ratio in overall_compression_ratio_stages:
            for pow in _power:
                if compression_ratio ** pow < max_stages:
                    compression_ratio_per_stages.append(compression_ratio ** pow)
                    break

        return compression_ratio_per_stages

    @staticmethod
    def get_compression_ratio(overall_compression_ratio):
        max_stages = len(_power)
        result = 0

        for pow in _power:
            if overall_compression_ratio ** pow < max_stages:
                result = overall_compression_ratio ** pow
                return result

    @staticmethod
    def get_num_of_compression_stages(overall_compression_ratio_stages):
        num_of_compression_stages = []
        compression_ratio_per_stages = Compressor.get_compression_ratio_stages(overall_compression_ratio_stages)

        for overall_compression_ratio, compression_ratio in \
                zip(overall_compression_ratio_stages, compression_ratio_per_stages):
            for pow in _power:
                if overall_compression_ratio ** pow == compression_ratio:
                    num_of_compression_stages.append(int(1 / pow))
                    break

        return num_of_compression_stages

    @staticmethod
    def get_num_of_compression(overall_compression_ratio):
        result = 0
        compression_ratio = Compressor.get_compression_ratio(overall_compression_ratio)

        for pow in _power:
            if overall_compression_ratio ** pow == compression_ratio:
                result = int(1 / pow)
                return result

    def compressor_energy(self,
                                          eta_compressor,
                                          overall_compression_ratio,
                                          stream,
                                          plot=0):
        """
        Calculate compressor energy consumption

        :param field: (Field)
        :param prime_mover_type:
        :param eta_compressor:
        :param overall_compression_ratio:
        :param inlet_stream: (Stream)
        :param inlet_tp: (TemperaturePressure) the T and P at the inlet to override
           that in the inlet_stream
        :return:
        """
        #energy_consumption = ureg.Quantity(0, "mmbtu/day")
        #tp = inlet_tp or inlet_stream.tp
        #inlet_temp, inlet_press = tp.get()
        
        rate = stream.gas_rate
        inlet_temp = stream.temperature
        inlet_press = stream.pressure
        fluid = stream.gas
        
        if overall_compression_ratio < 1 or max(rate) == 0:
            return 0*np.ones(len(rate))*ureg('kW'), inlet_temp, inlet_press

        compression_ratio = Compressor.get_compression_ratio(overall_compression_ratio)
        num_stages = Compressor.get_num_of_compression(overall_compression_ratio)
        total_work, outlet_temp, outlet_press = self.get_compressor_work_temp(inlet_temp,inlet_press,fluid,compression_ratio,num_stages)
        #print(str(inlet_temp)+' '+str(inlet_press)+' '+str(compression_ratio)+' '+str(num_stages))
        #volume_flow_rate_STP = field.gas.volume_flow_rate_STP(inlet_stream)
        #total_energy = total_work * volume_flow_rate_STP
        total_energy = total_work * rate
        
        brake_horse_power = total_energy / eta_compressor
        power = brake_horse_power.to('kW')
        #energy_consumption = get_energy_consumption(prime_mover_type, brake_horse_power)
        
        if self.save:
            self.field.POWER_DICT[self.identification] = power

        return power, outlet_temp, outlet_press
    
    
    
    def get_compressor_energy_consumption(self,
                                          eta_compressor,
                                          overall_compression_ratio,
                                          rate,
                                          inlet_temp,
                                          inlet_press,
                                          fluid,plot=0):
        """
        Calculate compressor energy consumption

        :param field: (Field)
        :param prime_mover_type:
        :param eta_compressor:
        :param overall_compression_ratio:
        :param inlet_stream: (Stream)
        :param inlet_tp: (TemperaturePressure) the T and P at the inlet to override
           that in the inlet_stream
        :return:
        """
        #energy_consumption = ureg.Quantity(0, "mmbtu/day")
        #tp = inlet_tp or inlet_stream.tp
        #inlet_temp, inlet_press = tp.get()

        if overall_compression_ratio < 1 or max(rate) == 0:
            return 0*np.ones(len(rate))*ureg('kW'), inlet_temp, inlet_press

        compression_ratio = Compressor.get_compression_ratio(overall_compression_ratio)
        num_stages = Compressor.get_num_of_compression(overall_compression_ratio)
        total_work, outlet_temp, outlet_press = self.get_compressor_work_temp(inlet_temp,inlet_press,fluid,compression_ratio,num_stages)
        #print(str(inlet_temp)+' '+str(inlet_press)+' '+str(compression_ratio)+' '+str(num_stages))
        #volume_flow_rate_STP = field.gas.volume_flow_rate_STP(inlet_stream)
        #total_energy = total_work * volume_flow_rate_STP
        total_energy = total_work * rate
        
        brake_horse_power = total_energy / eta_compressor
        power = brake_horse_power.to('kW')
        #energy_consumption = get_energy_consumption(prime_mover_type, brake_horse_power)
        
        #self.field.POWER_DICT[self.identification] = power

        return power, outlet_temp, outlet_press
    
    
    

