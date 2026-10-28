# -*- coding: utf-8 -*-

from .registry import ureg

import pandas as pd
import numpy as np

def crude_oil_transport(oil_transp):

    components = ['diluent', 'LNG', 'crude', 'petrocoke']
    modes = ['tanker', 'barge', 'pipeline', 'rail', 'truck']

    components = {
        'diluent':0,
        'LNG':0,
        'crude':0,
        'petrocoke':0,
    }

    fuel_edensity = { 
        'crude': 40.46,
        'diesel': 40.55,
        'LPG': 44.17,
        'residual_oil':37.41
    } #LHV/density

    ghg = {
        'diesel':78901,
        'residual_oil':85783,
        'ng':62048,
    } # [g/mmbtu]

    ghg_list = [78901, 85783, 62048]    # diesel, residual oil, ng
    emission_factors = [78901*ureg('metric_ton/Btu'), 85783*ureg('metric_ton/Btu'), 68194*ureg('metric_ton/Btu')]    # diesel, residual oil, ng

    fraction_transported_each_mode = [1,0,0] # tanker, barge, pipeline
    distances = [90*ureg('mi'), 0*ureg('mi'), 0*ureg('mi')] # tanker, barge, pipeline

    tanker_fuel = 'diesel'
    barge_fuel = 'diesel'

    load_factor_origin_tanker = 0.8
    load_factor_origin_barge = 0.8
    load_factor_dest_tanker = 0.7
    load_factor_dest_barge = 0.6

    barge_origin_econsumption = (14.42/load_factor_origin_barge+350)*0.735*fuel_edensity[barge_fuel]*ureg('Btu/(hp*hr)')
    tanker_origin_econsumption = (14.42/load_factor_origin_tanker+150)*0.735*fuel_edensity[tanker_fuel]*ureg('Btu/(hp*hr)')
    barge_dest_econsumption = (14.42/load_factor_dest_barge+350)*0.735*fuel_edensity[barge_fuel]*ureg('Btu/(hp*hr)')
    tanker_dest_econsumption = (14.42/load_factor_dest_tanker+150)*0.735*fuel_edensity[tanker_fuel]*ureg('Btu/(hp*hr)')

    barge_horsepower = 5600*ureg('hp')
    tanker_horsepower = 34320*ureg('hp')
    tanker_size = 35000*ureg('metric_ton')
    barge_size = 22500*ureg('metric_ton')
    tanker_speed = 19*ureg('mile/hr')
    barge_speed = 5*ureg('mile/hr')

    tanker_eintensity = (tanker_origin_econsumption*load_factor_origin_tanker + tanker_dest_econsumption*load_factor_dest_tanker)*tanker_horsepower/(tanker_speed*tanker_size)
    barge_eintensity = (barge_origin_econsumption*load_factor_origin_barge + barge_dest_econsumption*load_factor_dest_barge)*barge_horsepower/(barge_speed*barge_size)
    pipeline_eintensity = (240*0.55+270*0.36+260*0.09)*ureg('Btu/(mile*metric_ton)')

    tanker_econsumption = tanker_eintensity/18280*(1000/2)*ureg('metric_ton')
    barge_econsumption = barge_eintensity/18280*(1000/2)*ureg('metric_ton')
    pipeline_econsumption = pipeline_eintensity/18280*(1000/2)*ureg('metric_ton')

    fuel_share = {'fuel': ['Tanker', 'Barge', 'Pipeline'],
            'diesel': [0, 0, 0.45],
            'residual_oil': [1.0, 1.0, 0],
            'ng': [0, 0, 0.55]}

    fuel = {'fuel':['diesel','residual_oil','ng'],
            'Tanker': [1,0,0],
            'Barge': [1,0,0],
            'Pipeline': [0.45,0,0.55]}

    fuel_table = pd.DataFrame(fuel)

    tanker_gconsumption = (tanker_econsumption/1000000)*np.inner(list(fuel_table['Tanker']), ghg_list)*ureg('1/oil_bbl')
    barge_gconsumption = (barge_econsumption/1000000)*np.inner(list(fuel_table['Barge']), ghg_list)*ureg('1/oil_bbl')
    pipeline_gconsumption = (pipeline_econsumption/1000000)*np.inner(list(fuel_table['Pipeline']), ghg_list)*ureg('1/oil_bbl')

    econsumptions = [tanker_econsumption, barge_econsumption, pipeline_econsumption]

    ng_energy = sum([val * distances[i] * fraction_transported_each_mode[i] * fuel_share['ng'][i] for i, val in enumerate(econsumptions)]) * oil_transp * 1e-6 * (5.5*ureg('1/oil_bbl'))
    diesel_energy = sum([val * distances[i] * fraction_transported_each_mode[i] * fuel_share['diesel'][i] for i, val in enumerate(econsumptions)]) * oil_transp *1e-6 * (5.5*ureg('1/oil_bbl'))
    residual_oil_energy = sum([val * distances[i] * fraction_transported_each_mode[i] * fuel_share['residual_oil'][i] for i, val in enumerate(econsumptions)]) * oil_transp *1e-6 * (5.5*ureg('1/oil_bbl'))

    energies_consumption = [diesel_energy, residual_oil_energy, ng_energy]

    emissions = sum([val * emission_factors[i] for i, val in enumerate(energies_consumption)])*1e-6

    return emissions
