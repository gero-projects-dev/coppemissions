# -*- coding: utf-8 -*-

from .registry import ureg

import numpy as np
import cantera as ct
from pyromat import igtools as igt
import math
import pandas as pd
import matplotlib.pyplot as plt
import copy
from .gas_turbine import *

from pathlib import Path

fluid_data = pd.read_excel(
    Path(__file__).resolve().parents[2] / "data" / "fluid_data.xlsx"
)



def density(fluid): 
    """
    This function returns the density of a selected gas [g/ft³].

    Parameters
    ----------
    gas : str
        The name of the gas (CO2, CH4, N2O, C2H6, C3H8 or C4H10).
    """

    density_dict = {"N2":33.48, "O2":38.25, "CO2":52.61, "H2O":21.53, "CH4":19.17, "C2H6":35.94, "C3H8":52.70, "CO":33.48, "H2":2.41}

    return density_dict[fluid]*(ureg.g/ureg.ft**3)


def gwp(fluid):
    """
    This function returns the GWP of a selected gas [g CO2/g gas].

    Parameters
    ----------
    gas : str
        The name of the gas (CO2, CH4, N2O, C2H6, C3H8 or C4H10).
    """
    return float(fluid_data[fluid_data['fluid']==fluid]['gwp'].values[0])

def tau(fluid):
    """
    This function returns the stoichiometric relationship between a selected gas and product CO2 for complete combustion [g CO2/g gas].

    Parameters
    ----------
    gas : str
        The name of the gas (CO2, CH4, N2O, C2H6, C3H8 or C4H10).
    """
    return float(fluid_data[fluid_data['fluid']==fluid]['tau'].values[0])


def molar_mass(composition):
    molar_mass_dict = {"N2":35.70, "O2":31.25, "CO2":22.72, "H2O":55.52, "CH4":62.34, "C2H6":33.26, "C3H8":22.68, "CO":35.70, "H2":496.03}

    MM = 0
    for i in list(composition.keys()):
        if composition[i]>0:
            MM += composition[i]*molar_mass_dict[i]*(ureg.mol/ureg.kg)
    return MM


def mass_fraction_to_molar(mass_dict):

    base = {"N2":0, "O2":0, "CO2":0, "H2O":0, "CH4":0, "C2H6":0, "C3H8":0, "CO":0, "H2":0}
    gas = ct.Solution('gri30.yaml')
    gas.Y = mass_dict
    comp = gas.mole_fraction_dict()

    for i in list(comp.keys()):
        base[i] = comp[i]
        
    return base


def mass_fraction_to_molar2(mass_dict):

    molar_mass_dict = {"N2":35.70, "O2":31.25, "CO2":22.72, "H2O":55.52, "CH4":62.34, "C2H6":33.26, "C3H8":22.68, "C4H10":17.21, "CO":35.70, "H2":496.03}
    molar_mass_dict = {"N2":35.70, "O2":31.25, "CO2":22.72, "H2O":55.52, "CH4":62.34, "C2H6":33.26, "C3H8":22.68, "CO":35.70, "H2":496.03}
    molar = {"N2":0, "O2":0, "CO2":0, "H2O":0, "CH4":0, "C2H6":0, "C3H8":0, "CO":0, "H2":0}
    molar2 = {"N2":0, "O2":0, "CO2":0, "H2O":0, "CH4":0, "C2H6":0, "C3H8":0, "CO":0, "H2":0}

    for i in molar_mass_dict.keys():
        molar[i] = mass_dict[i]*molar_mass_dict[i]

    for i in molar_mass_dict.keys():
        molar2[i] = molar[i]/np.sum(list(molar.values()))

    return molar2
    

class Gas:
    """
    Class used to represent a gas model.
    """ 
    def __init__(self,composition):

        if abs(1 - sum(composition.values()))>0.01 :
            raise ValueError("Sum of composition must be 1. Sum equal "+ str(sum(composition.values())) )
        
        self.composition = composition

        self.mixture = " ".join(str(value) + " " + key + " +" for key, value in self.composition.items())
        self.mixture = self.mixture.rstrip(' +')

        self.mix = igt.IgtMix(self.mixture, units='mol')

        self.specific_gravity = (self.mix.mw()/28.96443)[0] 
        
        
    def density(self):

        density_dict = {"N2":33.48, "O2":38.25, "CO2":52.61, "H2O":21.53, "CH4":19.17, "C2H6":35.94, "C3H8":52.70, "CO":33.48, "H2":2.41}

        D = 0
        for i in list(self.composition.keys()):
            if self.composition[i]>0:
                D+=self.composition[i]*density_dict[i]*(ureg.metric_ton/(1e6*ureg.ft**3))
        return D




class Oil():
    def __init__(self, field):
        self.field = field
        self.a1 = 5.527215
        self.a2 = 0.783716
        self.a3 = 1.841408
    
    def specific_gravity(self):
        SG = 141.5 / (131.5 + self.field.API) * ureg("dimensionless")
        return SG

    



class Stream():
    """
    Class used to represent a stream.
    """
    def __init__(self, identificator, field, water_rate, oil_rate, gas_rate, T, P):
        self.gas = field.gas
        self.water_rate = water_rate
        self.oil_rate = oil_rate
        self.gas_rate = gas_rate
        self.pressure = P
        self.temperature = T
        self.oil = Oil(field)
        self.field = field

        self.a1 = 5.527215
        self.a2 = 0.783716
        self.a3 = 1.841408
        
        field.STREAMS[identificator] = self

    def solution_gor(self):
        R1 = self.pressure.to('psi')**(1/self.a2)*self.oil.specific_gravity()**(-self.a1/self.a2)
        R2 = (math.e**((self.a3/self.a2)*self.gas.specific_gravity*self.oil.specific_gravity()))
        R3 = (self.gas.specific_gravity * self.temperature.to('rankine'))
        return (R1*R2/R3).m*ureg('ft³/oil_bbl')
    
    def gas_mass_rate(self):
        mass_rate = (self.gas_rate * self.gas.density()).to('metric_ton/day')
        return mass_rate
    
    def gas_mass_fraction(self):
        base = {"N2":0, "O2":0, "CO2":0, "H2O":0, "CH4":0, "C2H6":0, "C3H8":0, "CO":0, "H2":0}
        gas = ct.Solution('gri30.yaml')
        gas.X = self.gas.composition
        comp = gas.mass_fraction_dict()

        for i in list(comp.keys()):
            base[i] = comp[i] 
        
        return base
    
    def gas_mass_fraction_rate(self):
        base = {"N2":0, "O2":0, "CO2":0, "H2O":0, "CH4":0, "C2H6":0, "C3H8":0, "CO":0, "H2":0}
        gas = ct.Solution('gri30.yaml')
        gas.X = self.gas.composition
        comp = gas.mass_fraction_dict()

        for i in list(base.keys()):
            if i in comp.keys():
                base[i] = comp[i] 
            base[i] = base[i] * self.gas_mass_rate()
        
        return base
    
    
    
    

class Field():
    """
    Class used to represent a Field.
    """
    def __init__(self,model_data_path):
        #self.data = pd.read_excel('dados campo presal.xlsx')
        
        model_data_path = str(Path(__file__).resolve().parents[2] / "data" / model_data_path)
        
        with open(model_data_path, 'r') as file:
                data_dict = file.read()
                
        self.data = eval(data_dict)

        field_data_path = Path(self.data['field_data_path'])

        if not field_data_path.is_absolute():
            field_data_path = (
                Path(__file__).resolve().parents[2]
                / "data"
                / field_data_path.name
            )
        
        field_data = pd.read_excel(field_data_path)
        
        self.field_start = field_data['Time'][0]
        self.oil_prod = (field_data['Oil Rate'].to_numpy()*ureg(self.data['field_data_units']['oil'])).to(ureg.oil_bbl/ureg.day)
        self.gas_prod = (field_data['Gas Rate'].to_numpy()*ureg(self.data['field_data_units']['gas'])).to(ureg.ft**3/ureg.day)
        self.water_prod = (field_data['Water Rate'].to_numpy()*ureg(self.data['field_data_units']['water'])).to(ureg.oil_bbl/ureg.day)
        self.water_inj = (field_data['Water Inj'].to_numpy()*ureg(self.data['field_data_units']['water_inj'])).to(ureg.oil_bbl/ureg.day)
        self.gas_inj = (field_data['Gas Inj'].to_numpy()*ureg(self.data['field_data_units']['gas_inj'])).to(ureg.ft**3/ureg.day)
        self.gas_lift = (field_data['Gas Lift'].to_numpy()*ureg(self.data['field_data_units']['gas_lift'])).to(ureg.ft**3/ureg.day)

        a = self.gas_lift
        b = a*0
        
        for i in range(len(a)):
            if i == 0:
                b[i] = a[0]
            elif (a[i]-a[i-1]) < 0:
                b[i] = a[0]
            else:
                b[i] = a[i]-a[i-1]
             
        self.gas_lift = b

        self.time_array = field_data['Time'].to_numpy()
        self.empty_field_array = np.ones(len(field_data))
        
        self.API = self.data['oil_API']
        self.Gas = Gas(self.data['gas_molar_fraction'])
        self.Oil = Oil(self)
        self.wellhead_t = self.data['wellhead_temperature']#*self.empty_field_array
        self.wellhead_p = self.data['wellhead_pressure']#*self.empty_field_array
        self.GOR = (self.gas_prod / self.oil_prod); self.GOR[np.isnan(self.GOR)] = 0
        self.WOR = self.water_prod / self.oil_prod; self.WOR[np.isnan(self.WOR)] = 0
        
        self.STREAMS = {}
        self.POWER_DICT = {}
        self.EMISSIONS = {}
        
        
    @property
    def VF_EMISSIONS(self):
        emiss1 = copy.deepcopy(self.EMISSIONS)
        del emiss1['Flaring']
        return emiss1
    
    @property
    def TOTAL_VF_EMISSIONS(self):
        emiss2 = copy.deepcopy(self.EMISSIONS)
        del emiss2['Flaring']
        vf_emiss2 = emiss2
        return sum(list(vf_emiss2.values()))
    
    @property
    def TOTAL_POWER(self):
        soma = sum(list(self.POWER_DICT.values()))
        return np.nan_to_num(soma, nan=0.0)
    
    @property
    def TURBINES_EMISSIONS(self):
        return ngt_emissions2(self.TOTAL_POWER*(1-self.data['electrification_percentage']))[0].to('metric_ton/day')
    
    @property
    def ELEC_EMISSIONS(self):
        return self.TOTAL_POWER*(self.data['electrification_percentage']*self.data['power_from_shore_footprint'])
    
    @property
    def TOTAL_POWER_EMISSIONS(self):
        return (self.TOTAL_POWER*(self.data['electrification_percentage']*self.data['power_from_shore_footprint']) + ngt_emissions2(self.TOTAL_POWER*(1-self.data['electrification_percentage']))[0]).to('metric_ton/day')
    
    @property
    def TOTAL_FVF_EMISSIONS(self):
        return sum(list(self.EMISSIONS.values()))
    
    @property
    def TOTAL_EMISSIONS(self):
        return self.TOTAL_POWER_EMISSIONS + self.TOTAL_FVF_EMISSIONS
    
    @property
    def BOE_ANNUAL(self):
        gas_exp = self.gas_prod
        return (gas_exp.to('ft**3/d')*self.data['BOE_gas'] + self.oil_prod)
        
    @property
    def BOE(self):
        return sum(self.BOE_ANNUAL)
        
    @property
    def IGEE_ANNUAL(self):
        return (self.TOTAL_EMISSIONS/self.BOE_ANNUAL).to('kg/oil_bbl')
    
    @property
    def IGEE(self):
        return (sum(self.TOTAL_EMISSIONS)/self.BOE).to('kg/oil_bbl')
    

    def graph_field_data(self):
        
        fig, ax1 = plt.subplots()
        scale = 'linear'
        
        ureg.define("million = 10**6")
        ureg.define("M = 10**3")
        
        x = self.time_array
        x = np.array(range(1,len(self.time_array)+1))
        y = self.oil_prod.to('M oil_bbl/day')
        y2 = self.water_prod.to('M oil_bbl/day')
        y3 = self.water_inj.to('M oil_bbl/day')
        
        y4 = self.gas_prod.to('million m**3/day')
        y5 = self.gas_inj.to('million m**3/day')
        y6 = self.gas_lift.to('million m**3/day')
        
        
        lns1 = ax1.plot(x,y,'k',label="Oil production [Mbbl/day]")
        lns2 = ax1.plot(x,y2,'g',label="Water production [Mbbl/day]")
        lns3 = ax1.plot(x,y3,'y',label="Water injection [Mbbl/day]")
        
        ax1.xaxis.grid(color='k', linestyle='-', linewidth=0.1)
        ax1.set_ylabel('Mbbl/day',color='black')
        ax1.set_xlabel('Year',color='black')
        #ax1.set_xticks(x, x, rotation='vertical')
        ax2 = ax1.twinx() 
        
        lns4 = ax2.plot(x,y4,'r',label="Gas production [MMSm³/day]")
        lns5 = ax2.plot(x,y5,'c',label="Gas injection [MMSm³/day]")
        lns6 = ax2.plot(x,y6,'b',label="Gas lift [MMSm³/day]")
        
        ax2.tick_params(axis='y', labelcolor='tab:red')
        ax2.set_ylabel('MMSm³/day',color='tab:red')
        lns = lns1+lns2+lns3+lns4+lns5+lns6
        labs = [l.get_label() for l in lns]
        ax1.legend(lns, labs, loc='upper center',bbox_to_anchor=(0.5, -0.15),ncol=2 )
        ax1.set_title('FIELD PRODUCTION CURVES')
        fig.tight_layout() 
        plt.show()
    
    

    def graph(self, attr, title, unit=0, total=1):
        
        if unit == 0:
            if type(attr) == dict:
                unit=list(attr.values())[0].units
            else:
                unit=attr.units
            
        fig, ax = plt.subplots()

        scale = 'linear'

        colours = ['y', 'g', 'r', 'b', 'm', 'c', 'brown', 'orange', 'gray', 'pink','olive','navy', 'limegreen']
        index = 0
        
        if type(attr) == dict:

            for process in attr.keys():
                lns1 = ax.plot(self.time_array,attr[process].to(unit),colours[index],label=process)
                index += 1
            
            if total:
                lns2 = ax.plot(self.time_array,sum(list(attr.values())).to(unit),'k',label="TOTAL")
                
        else:
            
            lns3 = ax.plot(self.time_array,attr.to(unit),colours[index],label=title)
        
        plt.grid()
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=1, fancybox=True, shadow=True)
        ax.set_title(title)
        #ax.set_xticks(self.field.time_array, self.field.time_array, rotation='vertical')
        ax.set_yscale(scale)
        ax.grid(color='k', linestyle='-', linewidth=0.1)
        plt.show()
        

    def graph_all(self):
        self.graph_field_data()
        
        #self.graph(self.STREAMS_GAS_RATES, 'GAS STREAMS', unit='foot ** 3 /year', total=0)
        self.graph(self.EMISSIONS['Flaring'], 'FLARING', unit='metric_ton/year')
        self.graph(self.POWER_DICT, 'POWER CONSUMPTION', unit='megawatt')
        #self.graph(self.POWER_EMISSIONS, 'POWER PRODUCING EMISSIONS', unit='metric_ton/year', total=1)
        self.graph(self.EMISSIONS, 'FLARING, VENTING, FUGITIVES', unit='metric_ton/year')

        
        fig, ax = plt.subplots()
        scale = 'linear'
        
        x = self.time_array
        y = self.TOTAL_EMISSIONS.to('metric_ton/year')
        y2 = self.TOTAL_POWER_EMISSIONS.to('metric_ton/year')
        y3 = self.TOTAL_FVF_EMISSIONS.to('metric_ton/year')
        
        lns1 = ax.plot(x,y,'k',label="TOTAL CO2eq EMISSIONS [metric_ton/year]")
        lns1 = ax.plot(x,y2,'g',label="TOTAL POWER EMISSIONS [metric_ton/year]")
        lns1 = ax.plot(x,y3,'y',label="TOTAL OTHER EMISSIONS [metric_ton/year]")
        
        # # Annotate each point
        # for i, txt in enumerate(y):
        #     plt.annotate(f'({x[i]}, {y[i]})', (x[i], y[i]), textcoords="offset points", xytext=(0,10), ha='center')
        
        plt.grid()
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=1, fancybox=True, shadow=True)
        ax.set_title('TOTAL EMISSIONS')
        #ax.set_xticks(self.field.time_array, self.field.time_array, rotation='vertical')
        ax.set_yscale(scale)
        ax.grid(color='k', linestyle='-', linewidth=0.1)
        plt.plot()
        
        
    def compare_all(self, comparation_file_path):
        
        questor = pd.read_excel(comparation_file_path)
  
    
  
        fig, ax = plt.subplots()
        
        x = self.time_array
        y = self.POWER_DICT['UTILITIES'].to('MW')
        z1 = np.array(questor['General demand'])*ureg('MW')
        
        lns1 = ax.plot(x,y,'k',label="General demand - COPPE [MW]")
        lns1 = ax.plot(x,z1,'k--', label="General demand- QUE$TOR [MW]")
        
        plt.grid()
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=1, fancybox=True, shadow=True)
        ax.set_title('GENERAL ELECTRIC DEMAND')
        ax.grid(color='k', linestyle='-', linewidth=0.1)
        plt.show()
        
        
        
        
        fig, ax = plt.subplots()
        
        x = self.time_array
        y = self.POWER_DICT['Oil_Export_Pump'].to('MW')
        z2 = np.array(questor['Oil export demand'])*ureg('MW')
        
        lns1 = ax.plot(x,y,'k',label="Oil export pump - COPPE [MW]")
        lns1 = ax.plot(x,z2,'k--', label="Oil export pump - QUE$TOR [MW]")
        
        plt.grid()
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=1, fancybox=True, shadow=True)
        ax.set_title('OIL EXPORT PUMP DEMAND')
        ax.grid(color='k', linestyle='-', linewidth=0.1)
        plt.show()
        
     
        

        fig, ax = plt.subplots()
        
        x = self.time_array
        y = self.POWER_DICT['Water_Injection'].to('MW')
        z3 = np.array(questor['Water injection demand'])*ureg('MW')
        
        lns1 = ax.plot(x,y,'k',label="Water Injection Pump - COPPE [MW]")
        lns1 = ax.plot(x,z3,'k--', label="Water Injection Pump - QUE$TOR [MW]")
        
        plt.grid()
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=1, fancybox=True, shadow=True)
        ax.set_title('WATER INJECTION DEMAND')
        ax.grid(color='k', linestyle='-', linewidth=0.1)
        plt.show()
        
        
        
        
        fig, ax = plt.subplots()
        
        x = self.time_array
        y3 = self.POWER_DICT['GAS_INJ'].to('MW') + self.POWER_DICT['GAS_LIFT'].to('MW') + self.POWER_DICT['FLASH_COMP'].to('MW') + self.POWER_DICT['Separator'].to('MW')
        z4 = np.array(questor['Compressor demand'])*ureg('MW')
        
       
        lns1 = ax.plot(x,y3,'b',label="Compressor demand - COPPE [MW]")
        lns1 = ax.plot(x,z4,'b--', label="Compressor demand - QUE$TOR [MW]")
        
       
        plt.grid()
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=1, fancybox=True, shadow=True)
        ax.set_title('COMPRESSOR DEMAND')
        ax.grid(color='k', linestyle='-', linewidth=0.1)
        plt.show()
        
        
       
        
        
        fig, ax = plt.subplots()
        
        x = self.time_array

        y5 = self.TOTAL_POWER.to('MW')
        z5 = z1+z2+z3+z4
        t5 = np.array([53.48, 74.42, 83.34, 89.09, 90.31, 90.40, 91.12, 91.11, 90.98, 90.56, 90.15, 89.62, 87.70, 87.88, 87.91, 87.06, 85.46, 84.66, 83.78, 83.35, 83.08, 82.16, 82.29, 81.99, 82.88, 81.00, 82.91, 80.89, 81.24, 81.24])
        x5 = np.array([31.41, 46.50, 54.65, 56.01, 57.55, 56.66, 52.38, 52.35, 54.87, 57.85, 58.11, 55.48, 52.42, 51.94, 49.94, 47.00, 44.75, 42.30, 41.19, 39.74, 38.51, 35.87, 36.09, 34.75, 34.52, 33.66, 33.53, 31.93, 31.64, 31.64])
        
        lns1 = ax.plot(x,y5,'k',label="TOTAL ELECTRICAL DEMAND - COPPE [MW]")
        lns1 = ax.plot(x,z5,'r', label="TOTAL ELECTRICAL DEMAND - QUE$TOR [MW]")
        lns1 = ax.plot(x,t5,'g--', label="ELECTRICAL DEMAND - DATA [MW]")
        lns1 = ax.plot(x,t5+x5,'g', label="POWER DEMAND - DATA [MW]")

        plt.grid()
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=1, fancybox=True, shadow=True)
        ax.set_title('TOTAL ELECTRICAL DEMAND')
        ax.grid(color='k', linestyle='-', linewidth=0.1)
        plt.show()
        
        
        fig, ax = plt.subplots()
        
        x = self.time_array

        y5 = self.TOTAL_EMISSIONS.to('metric_ton/year')*ureg('year')
        t5 = np.array([565217.4, 716724, 785979.1, 858122.5, 879228.4, 880496.5, 897125.4, 895806.4, 889257.7, 885304, 860441.3, 846255.1, 810660.2, 825107.9, 808777, 788893.5, 747608.8, 732481.1, 712048.7, 710474.1, 686615.6, 675485.3, 668465, 676393.5, 674660.3, 670654.9, 668541, 660375, 651743.2, 651743.2])*12*ureg('metric_ton')
        z5 = np.ones(30)*(1/30)*66576394*ureg('metric_ton')
        
        lns1 = ax.plot(x,y5,'k',label="TOTAL EMISSIONS - COPPE [metric_ton]")
        lns1 = ax.plot(x,t5,'g', label="TOTAL EMISSIONS - DATA [metric_ton]")
        lns1 = ax.plot(x,z5,'r', label="TOTAL EMISSIONS - QUE$TOR [metric_ton]")

        plt.grid()
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=1, fancybox=True, shadow=True)
        ax.set_title('TOTAL EMISSIONS')
        ax.grid(color='k', linestyle='-', linewidth=0.1)
        plt.show()
        
        
        
        fig, ax = plt.subplots()
        
        x = self.time_array

        parc1 = np.array([67.62518, 82.1361, 87.81653, 95.01861, 93.52084, 95.43841, 93.69271, 95.36477, 93.73277, 94.27104, 92.0119, 88.56653, 87.18527, 86.14491, 83.50506, 84.37504, 81.6787, 83.79683, 81.02927, 80.70525, 80.98568, 81.11069, 75.23415, 78.3139, 73.55982, 72.90893, 75.28651, 74.03559, 72.02469, 74.24649])
        parc2 = np.array([38.44334, 45.78835, 48.75329, 48.31655, 46.54373, 47.30037, 48.82793, 51.01354, 50.06837, 55.34365, 50.95427, 49.82512, 46.42252, 45.39971, 43.02707, 43.40365, 42.11122, 41.71362, 41.96494, 41.50468, 40.46444, 42.15014, 40.81361, 40.5158, 38.60834, 38.39255, 38.59294, 38.09718, 37.47242, 37.44782])
        
        lns1 = ax.plot(x,parc1,'r--', label="ELECTRICAL DEMAND REINJ PARCIAL [MW]")
        lns1 = ax.plot(x,parc1+parc2,'r', label="POWER DEMAND REINJ PARCIAL [MW]")        


        total1 = np.array([53.48, 74.42, 83.34, 89.09, 90.31, 90.40, 91.12, 91.11, 90.98, 90.56, 90.15, 89.62, 87.70, 87.88, 87.91, 87.06, 85.46, 84.66, 83.78, 83.35, 83.08, 82.16, 82.29, 81.99, 82.88, 81.00, 82.91, 80.89, 81.24, 81.24])
        total2 = np.array([31.41, 46.50, 54.65, 56.01, 57.55, 56.66, 52.38, 52.35, 54.87, 57.85, 58.11, 55.48, 52.42, 51.94, 49.94, 47.00, 44.75, 42.30, 41.19, 39.74, 38.51, 35.87, 36.09, 34.75, 34.52, 33.66, 33.53, 31.93, 31.64, 31.64])
        

        lns1 = ax.plot(x,total1,'g--', label="ELECTRICAL DEMAND REINJ TOTAL [MW]")
        lns1 = ax.plot(x,total1+total2,'g', label="POWER DEMAND REINJ TOTAL [MW]")

        plt.grid()
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=1, fancybox=True, shadow=True)
        ax.set_title('TOTAL ELECTRICAL DEMAND')
        ax.grid(color='k', linestyle='-', linewidth=0.1)
        plt.show()
        
        
        
   
    def compare_compressors(self, comparation_file_path):
        
        questor = pd.read_excel(comparation_file_path)
        
        fig, ax = plt.subplots()
        scale = 'linear'
        
        x = self.time_array
        y1 = self.POWER_DICT['GAS_INJ'].to('MW')
        y2 = self.POWER_DICT['GAS_LIFT'].to('MW')
        #y3 = self.POWER_DICT['Gas Processing'].to('MW')
        y4 = self.POWER_DICT['FLASH_COMP'].to('MW')
        y5 = self.POWER_DICT['Separator'].to('MW')
        y6 = y1+y2+y4+y5
        y7 = np.array(questor['Compressor demand'])*ureg('MW')
        
        lns1 = ax.plot(x,y1,'y',label="GAS INJECTION [MW]")
        lns2 = ax.plot(x,y2,'g',label="GAS LIFT [MW]")
       # lns3 = ax.plot(x,y3,'r',label="GAS PROCESSING [MW]")
        lns4 = ax.plot(x,y4,'b',label="FLASH GAS [MW]")
        lns5 = ax.plot(x,y5,'m',label="SEPARATOR [MW]")
        lns6 = ax.plot(x,y6,'k',label="TOTAL DEMAND [MW]")
        lns7 = ax.plot(x,y7,'k--',label="TOTAL DEMAND - QUE$TOR [MW]")
        
        plt.grid()
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=1, fancybox=True, shadow=True)
        ax.set_title('GAS COMPRESSORS DEMAND')
        #ax.set_xticks(self.field.time_array, self.field.time_array, rotation='vertical')
        ax.set_yscale(scale)
        ax.grid(color='k', linestyle='-', linewidth=0.1)
        plt.show()
        
        
        
        
        
        

        
        
        


    
        
        
        
        
