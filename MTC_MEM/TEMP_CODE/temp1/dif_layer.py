import matplotlib.pyplot as plt
import numpy as np
from CoolProp.CoolProp import PropsSI
from insulator import Insulation

Din = 0.034
Dd = 0.006
Th = 477
P = 70
Tc = 368
# Dds = np.arange(0.003, 0.01, 0.001)
# Tcs = np.arange(323, 403, 10)
res = []
# for Tc in Tcs:
x = np.array([0.005779845, 0.01912985, 0.001479454, 0.002374611, 0.000895157])
a = Insulation(Din + Dd * 2, Din, 1, 0)
r = a.flux(Th, 70, x, Tc)
m_ch3oh, m_h2o = r['mflux'][2], r['mflux'][3]
h_sen, h_lat = r['hflux'], r['hlg']
q = abs(h_sen / 1000 / (m_ch3oh * 32.04))
res.append([m_ch3oh, m_h2o, h_sen, q])  # kJ/g
print(np.linspace(Din / 2, (Din + Dd * 2) / 2, 200))
# res = np.array(res).round(4)
# # print(res[:, 2])
# # print(res[:,0])
# # q = res[:,2]/1000/(abs(res[:, 0])*32.04)
# plt.plot(Tcs, res[:,-1])
# plt.show()
# def func(x):
#     return x / (np.exp(x) - 1)
#
#
# # c = np.arange(-200, -0.05, 0.05)
# # res = np.zeros(len(c))
# #
# # for i in range(len(c)):
# #     res[i] = func(c[i])
# #
# # plt.plot(c, res)
# # plt.show()
#
# Ts = np.arange(323, 403, 10)
# Ps = np.zeros(len(Ts))
# for i in range(len(Ts)):
#     Ps[i] = PropsSI('Pp', 'T', Ts[i], 'Q', 1, 'Methanol') / 1e5
#
# alpha = Ps / Ts
# plt.plot(Ts, alpha)
# plt.show()
