# -*- coding: utf-8 -*-

globals().clear()

import coppemissions as ce
from coppemissions.registry import ureg
import numpy as np

print('\n############ INICIANDO COPPEMISSIONS (nome provisório)! ############\n\n')

# Field definition
CAMPO = ce.Field('../data/data_buzios_exp_reinj.txt')
Stream_0 = ce.Stream('Initial', CAMPO, CAMPO.Gas, CAMPO.water_prod, CAMPO.oil_prod, CAMPO.gas_prod, CAMPO.wellhead_t, CAMPO.wellhead_p)
print('FIELD DEFINITION: DONE!')

# Separator
SEP = ce.Separator(CAMPO, Stream_0)
Stream_1 = SEP.gas_outstream(); Stream_2 = SEP.oil_outstream(); Stream_3 = SEP.water_outstream()
print('SEPARATOR: DONE!')

# Oil Exportation
EXP_OIL_PUMP = ce.Pump("Oil_exportation", CAMPO, Stream_2.oil_rate, Stream_2.pressure, CAMPO.data['oil_export_pressure'], CAMPO.data['oil_export_pump_efficiency'])
print('OIL EXPORTATION: DONE!')

#Oil Storage
Crude_Storage = ce.Crude_Oil_Storage('Crude_Storage_1', CAMPO, Stream_2)
Flash_Recovered_Stream = Crude_Storage.recovered_stream
Flash_Flared_Stream = Crude_Storage.flared_stream

#FLashing
FLASH_COMP = ce.Compressor('Flashing', CAMPO)
FLASH_COMP.overall_compression_ratio = CAMPO.data['flash_gas_outlet_pressure']/CAMPO.data['flash_gas_initial_pressure']
CAMPO.POWER_DICT['Flashing'] = FLASH_COMP.compressor_energy(CAMPO.data['storage_compressor_efficiency'],FLASH_COMP.overall_compression_ratio,Flash_Recovered_Stream)[0]
print('FLASH GAS: DONE!')

# Seawater capture
Water_Capture_Rate = (CAMPO.water_inj - Stream_3.water_rate)*1 + 20e3*ce.ureg('oil_barrel / day')
WTR_CAP_PUMP = ce.Pump("Water_capture", CAMPO,Water_Capture_Rate,0*ce.ureg('psi'), CAMPO.data['water_capture_pressure'], CAMPO.data['water_capture_pump_efficiency'])
Stream_9 = ce.Stream('9', CAMPO, CAMPO.Gas, CAMPO.water_inj, 0*CAMPO.oil_prod, 0*CAMPO.gas_inj, 300*ce.ureg('K'), CAMPO.data['water_capture_pressure'])
print('WATER CAPTURE: DONE!')

# Water Injection
WI_PUMP = ce.Pump("Water_injection", CAMPO,Stream_9.water_rate,0*ce.ureg('psi'), CAMPO.data['water_injection_pressure'], CAMPO.data['water_injection_efficiency'])
print('WATER INJECTION: DONE!')

# Flare, Vent, Gas Processing
Stream_1.gas_rate += Flash_Recovered_Stream.gas_rate
FLARE = ce.Flare(CAMPO, Stream_1, volume_total=((Stream_1.gas_rate)*0.01+Flash_Flared_Stream.gas_rate))
Stream_4 = FLARE.gas_stream_to_vent()
Vent = ce.Vent("Vent1", CAMPO, Stream_4)
Stream_5 = Vent.outstream


ARQ = ['Gas Gathering', "Dehydration", "Membrane"]
GAS_PROC = ce.GasProcessing(CAMPO, Stream_5)
GAS_PROC.standard_parameters = CAMPO.data
Stream_6 = GAS_PROC.outstream
Stream_7 = GAS_PROC.CO2_outstream
print('GAS PROCESSING: DONE!')

# Gas Injection
Stream_inj = ce.Stream('Stream_inj', CAMPO, CAMPO.Gas,Stream_5.water_rate, Stream_5.oil_rate, CAMPO.gas_inj, Stream_5.temperature, CAMPO.data['gas_injection_inlet_pressure'])
GAS_INJ_COMP = ce.Compressor('Gas_Injection', CAMPO)
GAS_INJ_COMP.overall_compression_ratio = CAMPO.data['gas_injection_outlet_pressure']/CAMPO.data['gas_injection_inlet_pressure']
GAS_INJ_COMP.compressor_energy(CAMPO.data['gas_injection_compressor_efficiency'],GAS_INJ_COMP.overall_compression_ratio,Stream_inj)[0]
print('GAS INJECTION: DONE!')

# Gas Export
Stream_exp = ce.Stream('Stream_exp', CAMPO, CAMPO.Gas,Stream_6.water_rate, Stream_6.oil_rate, Stream_6.gas_rate - CAMPO.gas_inj - CAMPO.gas_lift, Stream_6.temperature, CAMPO.data['gas_export_inlet_pressure'])
EXP_GAS_COMP = ce.Compressor('Gas_Exportation', CAMPO)
EXP_GAS_COMP.overall_compression_ratio = CAMPO.data['gas_export_outlet_pressure']/ CAMPO.data['gas_export_inlet_pressure']
EXP_GAS_COMP.compressor_energy(CAMPO.data['gas_export_compressor_efficiency'],EXP_GAS_COMP.overall_compression_ratio,Stream_exp)[0]
 
#Gas Lift
Stream_lift = ce.Stream('Stream_lift', CAMPO, CAMPO.Gas,Stream_5.water_rate, Stream_5.oil_rate,CAMPO.gas_lift, Stream_5.temperature, CAMPO.data['gas_lift_inlet_pressure'])
GAS_LIFT_COMP = ce.Compressor('Gas_Lift', CAMPO)
GAS_LIFT_COMP.overall_compression_ratio = CAMPO.data['gas_lift_outlet_pressure']/CAMPO.data['gas_lift_inlet_pressure']
GAS_LIFT_COMP.compressor_energy(CAMPO.data['gas_lift_compressor_efficiency'],GAS_LIFT_COMP.overall_compression_ratio, Stream_lift)[0]
print('GAS LIFT: DONE!')

# UTILITIES
Utilities1 = ce.Utilities('Utilities1', CAMPO)
print('UTILITIES: DONE!')

# HEATING
Heating1 = ce.Heating('Heating1', CAMPO)
print('HEATING: DONE!')


print('')
print('TOTAL EMISSIONS = ' + str(sum(CAMPO.TOTAL_EMISSIONS.to('metric_ton/year'))*ureg('year')))
print('TOTAL FLARING EMISSIONS = ' + str(sum(CAMPO.EMISSIONS['Flaring'].to('metric_ton/year'))*ureg('year')))
print('TOTAL COMBUSTION EMISSIONS = ' + str(sum(CAMPO.TOTAL_POWER_EMISSIONS.to('metric_ton/year'))*ureg('year')))
print('TOTAL VENTING & FUGITIVES EMISSIONS = ' + str(sum(CAMPO.TOTAL_VF_EMISSIONS.to('metric_ton/year'))*ureg('year')))

print('')
print('IGEE (annual) = ' + str(CAMPO.IGEE_ANNUAL))
print('')
print('IGEE = ' + str(CAMPO.IGEE))

CAMPO.graph_all()