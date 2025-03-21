import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from CoolProp.CoolProp import PropsSI

from prop_calculator import VLEThermo, mixture_property

# F = np.array([0.015607824, 0.048960983, 0.0068335, 0.007209079, 0.001068755])
# F = np.array([0.013779679, 0.042448227, 0.007806613, 0.007006583, 0.000554595])
# F = np.array([0.015919268,0.050251206,0.006528082,0.007214911,0.001246702])
sub2 = ["CO2", "H2", "Methanol", "H2O", "CO"]
F = np.array([0.014955934, 0.046716216, 0.0074739200, 0.007482738, 0.000924207])
F2 = np.array([0.024463369, 0.073390108, 0, 0, 0])
F3 = np.array([0.02006272, 0.065285, 0.00185223, 0.00440065, 0.00254842])
F4 = np.array([58, 11, 12, 7, 12]) / np.array([44, 2, 32, 18, 32])
print(F4/np.sum(F4))
P = 75
a = VLEThermo(comp=sub2, ref='ev')
print(a.T_dew(P, x=F4))
flasher_gas, flasher_liq, sf = a.flash(60+273, P, F4)
flasher_gas = np.array(flasher_gas)
flasher_liq = np.array(flasher_liq)
print(sf)
F_gas = np.sum(F) * sf * flasher_gas
F_liq = np.sum(F) * (1 - sf) * flasher_liq
print(F_liq / F)


print(a.T_dew(P=25, x=F3))
Hs = []
Hs2 = []
Ts = np.arange(403, 493, 1)
for T in Ts:
    Hs.append(a.cal_H(T, P, F))

Ts2 = np.arange(393, 483, 1)
for T in Ts2:
    Hs2.append(a.cal_H(T, P, F2))
Hs = np.array(Hs) - Hs[0]
Hs2 = np.array(Hs2) - Hs2[0] + 0.27
plt.plot(Hs, Ts, 'r')
plt.plot(Hs2, Ts2, 'b')
plt.show()

T_dew = a.T_dew(P=P, x=F)
phi = a.phi(T=503, P=75, x=F)
xi = F / np.sum(F)
fi = P * (F / np.sum(F)) * phi
print(fi)
print(xi)
print(T_dew)
flasher_gas, flasher_liq, sf = a.flash(393, P, F)
flasher_gas = np.array(flasher_gas)
flasher_liq = np.array(flasher_liq)
print(sf)
F_gas = np.sum(F) * sf * flasher_gas
F_liq = np.sum(F) * (1 - sf) * flasher_liq

print(F_liq)
print(F_liq / F)
print(F_gas / np.sum(F_gas))
