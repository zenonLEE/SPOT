import numpy as np
import pandas as pd
from prop_calculator import VLEThermo
from utility import HX4Water, HeatExchanger, DistillerOpt, heater, heater_duty, flasher

sub = ['CO2', 'H2', 'Methanol', 'H2O', 'CO', 'N2']
ht_stream = pd.Series([0,
0,
0,
0.5,
0,
0,
530.9308267,
48
], index=sub + ['T', 'P'])

lt_stream = pd.Series([0.008294132,
0.144731216,
5.48115E-05,
5.72703E-07,
0.001610471,
0,
416.645,
84.18
], index=sub + ['T', 'P'])

hx_ht = HeatExchanger(ht_stream, lt_stream, U=500)
hx_ht_res = hx_ht.fixed_T_cold(Tout=522.39)
print(hx_ht_res)

hx_ht_res = hx_ht.fixed_T_cold(Tout=522.39)
print(hx_ht_res)