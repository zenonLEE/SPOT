from thermo import (ChemicalConstantsPackage, SRKMIX, FlashVL, CEOSLiquid, CEOSGas, HeatCapacityGas,
                    GibbsExcessLiquid, MSRKMIX, MSRKMaMIX)
from thermo.nrtl import NRTL
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.constants import R
from CoolProp.CoolProp import PropsSI
from t5 import *

subs = ["CO2", "H2", "Methanol", "H2O", "carbon monoxide"]
constants, correlations = ChemicalConstantsPackage.from_IDs(subs)
n_count = len(subs)

kijs_srk = np.array([[0, -0.3462, 0.0148, 0.0737, 0],
                     [-0.3462, 0, 0, 0, 0.0804],
                     [0.0148, 0, 0, -0.0789, 0],
                     [0.0737, 0, -0.0789, 0, 0],
                     [0, 0.0804, 0, 0, 0]])

kijs_msrk = np.array([[0, 0.1164, 0.1, 0.3, 0.1164],
                      [0.1164, 0, -0.125, -0.745, -0.0007],
                      [0.1, -0.125, 0, -0.075, -0.37],
                      [0.3, -0.745, -0.075, 0, -0.474],
                      [0.1164, -0.0007, -0.37, -0.474, 0]])
#
# kijs = np.array([[0, -0.075],
#                  [-0.075, 0]])  # -0.075

p = [0, 0, 0.2359, 0.1277, 0]

eos_kwargs_msrk = dict(Tcs=np.array(constants.Tcs).tolist(), Pcs=np.array(constants.Pcs).tolist(),
                       omegas=np.array(constants.omegas).tolist(), kijs=kijs_msrk.tolist(), S2s=np.array(p).tolist())
eos_kwargs_srk = dict(Tcs=np.array(constants.Tcs).tolist(), Pcs=np.array(constants.Pcs).tolist(),
                      omegas=np.array(constants.omegas).tolist(), kijs=kijs_srk.tolist())

# thermo.HeatCapacityGas(CASRN=sub_CAS, MW=constants.MWs[0], method='TRCIG')
cp_cal = [HeatCapacityGas(CASRN=constants.CASs[i], MW=constants.MWs[i], method='TRCIG') for i in range(len(subs))]

# Pp = 5e6
Ts = np.arange(298.15, 598.15, 10)
zs = np.array([0.5, 0.5])  # np.array([0, 0, 0.5, 0.5, 0])  # np.array([0.25, 0.75, 0, 0, 0])
zs_l = np.array([0.3, 0.55, 0, 0.1, 0.05])

T = 523
P = 5E6
gas = CEOSGas(MSRKMIX, HeatCapacityGases=cp_cal, eos_kwargs=eos_kwargs_msrk, zs=zs_l.tolist(), T=T, P=P)
liquid = CEOSLiquid(MSRKMIX, HeatCapacityGases=cp_cal, eos_kwargs=eos_kwargs_msrk, zs=zs_l.tolist(), T=T, P=P)
liquid_srk = CEOSLiquid(SRKMIX, HeatCapacityGases=cp_cal, eos_kwargs=eos_kwargs_srk)
gas_srk = CEOSGas(SRKMIX, HeatCapacityGases=cp_cal, eos_kwargs=eos_kwargs_srk)

flasher_srk = FlashVL(constants, correlations, liquid=liquid_srk, gas=gas_srk)
P_dew = flasher_srk.flash(zs=zs_l.tolist(), T=T, VF=1).Pp / 1E5
print(P_dew)
PT = flasher_srk.flash(T=T, P=P, zs=zs_l.tolist())


# print(PT.gas.zs)
# print(PT.liquid0.zs)
# print(PT.VF)
print(gas.H(), liquid.H())
print(PT.H(), PT.G())


g_z = PT.gas.zs
l_z = PT.liquid0.zs
# print(g_z)

gas_f = CEOSGas(SRKMIX, HeatCapacityGases=cp_cal, eos_kwargs=eos_kwargs_srk, zs=g_z, T=T, P=P)
liquid_f = CEOSLiquid(SRKMIX, HeatCapacityGases=cp_cal, eos_kwargs=eos_kwargs_srk, zs=l_z, T=T, P=P)
print(gas_f.H(), gas_f.G())
print(liquid_f.H(), liquid_f.G())
print(gas_f.H() * PT.VF + liquid_f.H() * (1 - PT.VF))
#
# P_dew_msrk = flasher_msrk.flash(zs=zs_l, T=T, VF=1).Pp
#
# flasher_srk = FlashVL(constants, correlations, liquid=liquid_srk, gas=gas_srk)
# P_dew_srk = flasher_srk.flash(zs=zs_l, T=T, VF=1).Pp

# print(P_dew_msrk / 1e6, P_dew_srk / 1e6)
