import numpy as np
import pandas as pd
from prop_calculator import mixture_property, VLEThermo
from reactor_HMT import HMT

T = 511
P= 69.2#80
# F = [0.006229502,
#      0.02962911,
#      0.000145672,
#      2.86782E-05,
#      0.001728497,
#      0
#      ]
F = [0.0593,
     0.8609,
     0.0260,
     0.0267,
     0.0271,
     0
     ]
F = F/np.sum(F)
ds = 6e-3
Ft = 1e5*25.9e-3/8.314/273.15/4 #0.49 # mol/s

F = Ft*F
# rho_mol = 944.47
# rho_kg = 11.819

sub = ['CO2', 'H2', 'Methanol', 'H2O', 'CO', 'N2']
a = VLEThermo(comp=sub)

b = HMT(sub=sub, phi=0.55, ds=ds, Dt=0.042)
print(b.htc(T=T, P=P, F=F))
# print(b.htc2(T=503, P=40, F=F))
print(b.htc4(T=T, P=P, F=F))

# 0.4432887043471933
# 7.853981633974483e-05
# 1377.4903566217004
# -0.04755227349834889
