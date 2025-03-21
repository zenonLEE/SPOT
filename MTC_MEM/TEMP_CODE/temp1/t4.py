import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import thermo
from fluids.constants import R

comps = ["CO2", "H2", "Methanol", "H2O", "carbon monoxide"]  # ["Methanol"]#

Ts = np.arange(298.15, 599.15, 10)  # (298.15, 599.15, 10)
P = 7E6
sub = 'H2O'

const, correlations = thermo.ChemicalConstantsPackage.from_IDs(IDs=comps)

for i in range(5):
    cp_cal = thermo.HeatCapacityGas(CASRN=const.CASs[i], MW=const.MWs[i], method='POLING_POLY')
    print(cp_cal.correlations)
    print(comps[i])
    print(const.Hfgs[i])
    print(const.Sfgs[i])
    print(cp_cal.T_dependent_property(298.15))


