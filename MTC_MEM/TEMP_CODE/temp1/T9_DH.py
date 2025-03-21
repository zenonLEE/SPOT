import numpy as np

from prop_calculator import mixture_property, VLEThermo
from thermo import (ChemicalConstantsPackage, SRKMIX, FlashVL, CEOSLiquid, CEOSGas, HeatCapacityGas,
                    FlashVLN)
import pandas as pd

Tr = 483
P = 70

comp_list = ["CO2", "H2", "Methanol", "H2O", "carbon monoxide"]

F2 = np.array([0.008154456, 0.024463369, 0, 0, 0])
F1 = np.array([0.005811954, 0.019464184, 0.001328341, 0.002342503, 0.001014162])
T1, T2 = 547.344686, 483

F3 = np.array([0.006962447, 0.022439809, 0.000415775, 0.00119201, 0.000776235])
T3 = 476.2420782
a = VLEThermo(comp_list)
H1 = a.cal_H(T1, P, F1)
H2 = a.cal_H(T2, P, F2)
print(H2)
print(H1)

H3 = a.cal_H(T3, P, F3)
print(H3)

# ASPEN sim
F4 = np.array(
    [0.00694128792946486, 0.0222581752550055, 0.000496012837229712, 0.00121316807053521, 0.000717155233305423])
T4 = 211.373 + 273.15
H4 = a.cal_H(T4, P, F4)
print((H4-H1)/H1)

F5 = np.array([0.006962934, 0.022441264, 0.000415291, 0.001191523, 0.000776232])
T5 = 476.21548
H5 = a.cal_H(T5, P, F5)
print(H5)

F6 = np.array([0.005820089, 0.019007731, 0.001560635, 0.002334368, 0.000773733])
T6 = 538.8602348923272
H6 = a.cal_H(T6, P, F6)
print(H6)

F7 = np.array([0.006849922, 0.022009185, 0.000574825, 0.001304534, 0.000729709])
T7 = 486.6960719
H7 = a.cal_H(T7, P, F7)
print(H7)

F8 = np.array([0.006850455,0.022010793,0.000574287,0.001304002,0.000729715])
T8 = 486.6691731
H8 = a.cal_H(T8, P, F8)
print((H8-H1)/H1)