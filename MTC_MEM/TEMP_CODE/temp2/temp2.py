import pandas as pd
from prop_calculator import VLEThermo
from utility import HeatExchanger,heater

hs = [0.004077228,
0.001200559,
0,
0.005402534,
0,
0,
467.0370477,
1


]

ls = [0,
0,
0.007725089,
0.00745682,
0,
0,
351.002,
1.2
]
sub = ['CO2', 'H2', 'Methanol', 'H2O', 'CO', 'N2','T','P']
hs_pd = pd.Series(hs,index=sub)
ls_pd = pd.Series(ls, index=sub)

a = VLEThermo(sub[:-2],ref='ev')
Td = a.T_dew(1.2, ls[:-2])
Tp = a.T_pub(1.2, ls[:-2])
print(Td, Tp)
# hx_ht = HeatExchanger(hs_pd, ls_pd)
# hx_ht_res = hx_ht.fixed_delta_hc(2)
# print(hx_ht_res)
# hs_heater = heater(hs_pd,T_out=503.15)
#
# ls_heater = heater(ls_pd,T_out=503.15)
# print(hs_heater)
# print(ls_heater)
# a = heater(hs_pd, T_out=355.102000)
# b = heater(hs_pd, T_out=308.15)
#
# print(a,b)