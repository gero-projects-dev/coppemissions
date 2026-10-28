# -*- coding: utf-8 -*-

from .registry import ureg
import copy
import numpy as np
import math
from .pump import *

class Water_Injection():
    """
    This class represents the Water Injection process.

    Attributes:
    ----------
    field: Field object
    stream : Stream object
    data_dict : Dictionary

    data_dict keys
    ----------
    Q_w : Quantity object
        Total water flow to be injected [bbl/d].
    N_inj : int
        Number of injection wells.
    field_depth : Quantity object
        Field depth [ft].
    diameter : Quantity object
        Tube diameter [ft].
    field : Field object
        Field object.
    PI : Quantity object
        Productivity Index [bbl/d*psi]. 
    P_res : Quantity object
        Reservoir pressure [psi] - Standard value = 0 (when is zero, the function estimates)

    Methods:
    ----------
    power : Returns the power consumption of the process [kW].

    """

    stand_dict = {
        "injection_fields":3,
        "field_depth": 7240*ureg('ft'),
        "diameter": 0.23*ureg.ft,
        "PI": 6.8*ureg('oil_bbl/(day * psi)'),
        "P_res": 0,
        "pump_efficiency":0.55
    }

    def __init__(self, field, stream, data_dict = stand_dict):
        self.field = field
        self.stream = stream
        self.N_inj = data_dict["injection_fields"]
        self.field_depth = data_dict["field_depth"]
        self.diameter = data_dict["field_depth"]
        self.PI = data_dict["PI"]
        self.P_res = data_dict["P_res"]
        self.pump_efficiency = data_dict["pump_efficiency"]

    def power(self):

        Q_w = copy.copy(self.stream.water_rate)
        field = copy.copy(self.field)
        N_inj = self.N_inj
        field_depth = copy.copy(self.field_depth)
        diameter = copy.copy(self.diameter)
        PI = copy.copy(self.PI)
        P_res = copy.copy(self.P_res)

        field_age = []

        for i in field.time_array:
            field_age.append(float(i-field.time_array[0])*1e-9/(60*60*24*365))
        field_age = np.array(field_age)
        
        if P_res==0:
            P_res = (4.9e-8*field.age**2 + 0.95)*field_depth / 2.31
            P_res = (P_res.m)*ureg('psi') #resevoir default pressure
        
        Q_ww = Q_w / N_inj
        P_required = P_res + Q_ww/PI #required pressure to ensure that the water will flow into the formation 

        P_hs = field_depth * 0.433 * ureg('psi/ft') * field.empty_field_array # hydrostatic pressure of water at the depth of the reservoir

        f = 0.02 # Moody friction factor (default value)
        flow_velocity = (Q_ww * 4) / (math.pi * diameter**2)
        flow_velocity = flow_velocity.to('ft/s').m
        
        P_fric = 0.433 * (f * field_depth *flow_velocity**2) / (2 * diameter * 32.2) * ureg('psi')
    
        P_pump = P_required + P_fric - P_hs

        for i in range(len(P_pump)):
            if P_pump[i] < 0:
                P_pump[i] =0*ureg('psi') #the weight of the column of water in the injection wellbore is sufficient to drive the water into the formation on its own

        #print(P_pump)
        PUMP = Pump(field, Q_w/N_inj, self.stream.pressure, P_pump, self.pump_efficiency)

        return N_inj * PUMP.power()

