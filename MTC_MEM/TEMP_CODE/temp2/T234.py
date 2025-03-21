import numpy as np
import pandas as pd
from prop_calculator import VLEThermo
from utility import HX4Water, HeatExchanger, DistillerOpt, heater, heater_duty, flasher

sub = ['CO2', 'H2', 'Methanol', 'H2O', 'CO', 'N2']
ht_stream = pd.Series([0.000184121,
0.010554783,
0.008137011,
1.53147E-05,
0.001638283,
0,
448.175,
84.16805434
], index=sub + ['T', 'P'])

lt_stream = pd.Series([1.69316E-05,
3.13029E-07,
0.008057836,
1.52735E-05,
1.53289E-07,
0,
309.942,
1.2
], index=sub + ['T', 'P'])

vle_cal = VLEThermo(sub)
try:
    T_dew_in = vle_cal.T_dew(lt_stream['P'], lt_stream[sub].values)
except ValueError:
    T_dew_in = 365.15
print(T_dew_in)
hx_ht = HeatExchanger(ht_stream, lt_stream, U=500)
hx_ht_res = hx_ht.fixed_delta_cold(delta_T=20)
print(hx_ht_res)

hx_ht_res = hx_ht.fixed_T_cold(Tout=T_dew_in+10)
print(hx_ht_res)