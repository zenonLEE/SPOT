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

P = 75
a = VLEThermo(comp=sub2, ref='ev')
print(a.cal_cp(523, 50, F))
print(a.cal_cp2(523, 50, F))
print(a.cal_cp3(523, 50, F))
# flasher_gas, flasher_liq, sf = a.flash(100 + 273, Pp, F)

# phi1 = a.phi(100 + 273, Pp, flasher_gas)
# print(phi1)
# print(phi1 * (flasher_gas / np.sum(flasher_gas) * Pp))

# kijs = np.array([[0, -0.0789],
#                  [0, -0.0789]])
F_liq = np.array([0, 0, 0.4, 0.6, 0])
x_liq = F_liq / np.sum(F_liq)

# sub_liq = ["Methanol", "H2O"]
# a_liq = VLEThermo(comp=sub_liq, kij=kijs)
print(a.p_bub(100 + 273, F_liq), a.p_dew(100 + 273, F_liq))
phi = a.phi(100 + 273, P, F_liq)
print(phi * (F_liq / np.sum(F_liq) * P))
