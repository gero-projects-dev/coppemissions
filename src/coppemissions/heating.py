# -*- coding: utf-8 -*-

from .registry import ureg

class Heating():
    
    def __init__(self, identification, field):
        self.power = field.data['electrical_heating?']*ureg('MW')*(1/field.data['heater_efficiency'])*field.data['heating_demand']
        field.POWER_DICT[identification] = self.power
