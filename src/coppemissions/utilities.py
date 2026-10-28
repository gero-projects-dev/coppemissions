# -*- coding: utf-8 -*-

class Utilities():
    
    def __init__(self, identification, field):
        self.power = field.data['utilities_total_power']*field.empty_field_array
        field.POWER_DICT[identification] = self.power