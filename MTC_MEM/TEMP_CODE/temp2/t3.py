import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from prop_calculator import VLEThermo
from utility import HX4Water, HeatExchanger, DistillerOpt, heater, heater_duty

Tlt = 369.546000
Plt = 39.984576

sub = ['CO2', 'H2', 'Methanol', 'H2O', 'CO', 'N2']
ht_stream = pd.Series([0,
0,
0,
0.5,
0,
0,
519.258,
48
], index=sub + ['T', 'P'])
ht_heater = heater(ht_stream, T_out=513.15)
qm = ht_heater['Q']
print(qm)
liq_heated_ult = pd.Series([3.81176E-05,
7.78924E-08,
0.007886573,
0.001047998,
1.38649E-08,
0,
362.298,
1.2
],
                           index=sub + ['T', 'P'])


a = VLEThermo(sub,ref='ev')
Td = a.T_dew(1.2, liq_heated_ult[:-2])
Tp = a.T_pub(1.2, liq_heated_ult[:-2])
print(Td, Tp)
Tco = np.arange(0,50, 10) + liq_heated_ult['T']
res = np.zeros((len(Tco),2))
i = 0
for Tc in Tco:
    lt_heater = heater(liq_heated_ult, T_out=Tc)
    fin = lt_heater['Fo']
    apk_path = r"D:\study\00课题\06多联产系统\甲醇单元\反应器比较\膜反应器\DT_with_MR_opt_boiler.bkp"
    block_name = 'B3'
    feed_name = 'S18'
    heavy_name = 'S8'
    light_name = 'S7'
    sp = DistillerOpt(apk_path, block_name=block_name, feed_name=feed_name,
                      light_name=light_name, heavy_name=heavy_name)

    sp_res = sp.run_RF(stream=fin, valid_comp=sub, ht_stream=ht_stream)
    q_duty = sp_res['block']['HD']
    res[i] = [sp_res['ht']['T'], sp_res['block']['HD']]
    i += 1
    # res.append()
print(res)
plt.plot(Tco, res[:,0])
plt.show()
plt.plot(Tco, res[:,1])
plt.show()
# q_diff = q_duty + qm
# lt_pre_heater = heater_duty(liq_heated_ult, q_diff * 2)
# print(lt_pre_heater)
#
# sp = DistillerOpt(apk_path, block_name=block_name, feed_name=feed_name,
#                   light_name=light_name, heavy_name=heavy_name)
# sp_res2 = sp.run_RF(stream=lt_pre_heater['Fo'], valid_comp=sub, ht_stream=ht_stream)
# print(sp_res2)
# a = VLEThermo(sub)
# F = liq_heated_lt.values[:6]
# F1 = [0, 0, 1, 0, 0, 0]
# print(a.T_dew(1.2, F1))
# print(a.T_pub(1.2, F1))
# F2 = [0, 0, 0, 1, 0, 0]
# print(a.T_dew(1.2, F2))
# print(a.T_pub(1.2, F2))
# F3 = [1.667324e-05, 0,  3.877844e-03, 1.355306e-03, 0, 0]
# print(a.T_dew(1.2, F3))
# print(a.T_pub(1.2, F3))
# # h1 = a.cal_H(355)
# ht_hx = HX4Water(ht_stream, liq_heated_lt, U=500)
# ht_hx_res = ht_hx.fixed_T_cold(Tout=373.15)
# print(ht_hx_res)
