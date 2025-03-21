import numpy as np
import pandas as pd
from CoolProp.CoolProp import PropsSI
from prop_calculator import mixture_property
from prop_calculator import VLEThermo
from utility import mixer, multi_comp_opt, heater, valve
from reactor_HMT import HMT

F1 = np.array([0.012026285,
               0.078445202,
               0.000357675,
               9.33912E-06,
               0.003516678,
               0,
               513.15,
               50
               ])
sub = ["CO2", "H2", "Methanol", "H2O", "carbon monoxide"]
F2 = np.array([0, 0.001281415, 0, 0, 0, 353.15, 30])
F2 = pd.Series(F2, sub + ['T', 'P'])
res = valve(F2, P=1)

a = VLEThermo(comp=sub)
F2_valve = res['Fo']
ex_bf = a.cal_E(F2[-2],F2[-1],F2[:-2])
ex_valve = a.cal_E(F2_valve[-2],F2_valve[-1],F2_valve[:-2])
#kW
print(ex_bf-ex_valve)
# # [0.000233889,0.000769596,0.00012733,6.88411E-05,3.39642E-05]
# # [0.005779793,0.019129694,0.001479506,0.002374663,0.000895157]

# # F1 = pd.Series(F1, index=sub + ['T', 'P'])
# F2 = pd.Series(F2, index=sub + ['T', 'P'])
#
# a = HMT(sub=sub,phi=0.5,ds=5e-3,Dt=0.035,Dm=0.02)
# res = a.htc4(T=F1[-2],P=F1[-1],F=F1[:-3])
# print(res)
# res = mixer(F1, F2)
# # print(res)
# a = VLEThermo(sub)
# rho = a.cal_rho_mol(res['Fo']['T'], res['Fo']['P'], res['Fo'][:-2])
# # print(rho)
# # F2['P'] = 5
# # comp_res = multi_comp_opt(F2, 75, T_cool=308.15, r_max=3.5, r_min=2.5)
# # print(comp_res)
#
#
# F3 = np.array([8.93E-04, 7.46E-03, 0.00E+00, 0.00E+00, 0.00E+00, 120 + 273.15, 75])
# F3 = pd.Series(F2, index=sub + ['T', 'P'])
# res = heater(F3, 473.15)
# duty = res['Q']
# cpm = 4536.3
# qm = 0.01
# q_source = cpm * qm * 0.12 / 1000
# print(res)
# print(q_source)
