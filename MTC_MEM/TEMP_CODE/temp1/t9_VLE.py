from thermo import (ChemicalConstantsPackage, SRKMIX, FlashVL, CEOSLiquid, CEOSGas, HeatCapacityGas,
                    GibbsExcessLiquid, MSRKMIX,MSRKMaMIX)
from thermo.nrtl import NRTL
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.constants import R
from CoolProp.CoolProp import PropsSI
from t5 import *

subs = ["Methanol", "H2O"]
constants, correlations = ChemicalConstantsPackage.from_IDs(subs)
n_count = len(subs)

# kijs = np.array([[0, -0.3462, 0.0148, 0.0737, 0],
#                  [-0.3462, 0, 0, 0, 0.0804],
#                  [0.0148, 0, 0, -0.0789, 0],
#                  [0.0737, 0, -0.0789, 0, 0],
#                  [0, 0.0804, 0, 0, 0]])

kijs_srk = np.array([[0, -0.0789],
                     [-0.0789, 0]])

kijs = np.array([[0, -0.075],
                 [-0.075, 0]])  # -0.075

p = [0.2359, 0.1277]

eos_kwargs_msrk = dict(Tcs=np.array(constants.Tcs), Pcs=np.array(constants.Pcs),
                       omegas=np.array(constants.omegas), kijs=kijs, S2s=np.array(p))
eos_kwargs_srk = dict(Tcs=np.array(constants.Tcs), Pcs=np.array(constants.Pcs),
                      omegas=np.array(constants.omegas), kijs=kijs_srk)

# thermo.HeatCapacityGas(CASRN=sub_CAS, MW=constants.MWs[0], method='TRCIG')
cp_cal = [HeatCapacityGas(CASRN=constants.CASs[i], MW=constants.MWs[i], method='TRCIG') for i in range(len(subs))]

P = 5e6
Ts = np.arange(298.15, 598.15, 10)
zs = np.array([0.5, 0.5])  # np.array([0, 0, 0.5, 0.5, 0])  # np.array([0.25, 0.75, 0, 0, 0])
zs_l = np.array([0.5, 0.5])

T = 523
gas = CEOSGas(MSRKMaMIX, HeatCapacityGases=cp_cal, eos_kwargs=eos_kwargs_msrk)
liquid = CEOSLiquid(MSRKMaMIX, HeatCapacityGases=cp_cal, eos_kwargs=eos_kwargs_msrk)
gas_msrk = CEOSGas(MSRKMIX, HeatCapacityGases=cp_cal, eos_kwargs=eos_kwargs_msrk)
liquid_msrk = CEOSLiquid(MSRKMIX, HeatCapacityGases=cp_cal, eos_kwargs=eos_kwargs_msrk)
liquid_srk = CEOSLiquid(SRKMIX, HeatCapacityGases=cp_cal, eos_kwargs=eos_kwargs_srk)
gas_srk = CEOSGas(SRKMIX, HeatCapacityGases=cp_cal, eos_kwargs=eos_kwargs_srk)

flasher = FlashVL(constants, correlations, liquid=liquid, gas=gas)
# flasher.d


ref_path = r"D:\document\00Study\05多联产系统\甲醇单元\Gibbs\热力学性质_ASPEN.xlsx"
ref_data = pd.read_excel(ref_path, sheet_name='VLE_523')  # CO2_H2  CH3OH_H2O

z1, z2, P_d, P_b = flasher.plot_Pxy(T=T, pts=10, show=False, values=True, zmax=1)
# print(P_d)

plt.scatter(ref_data['Zb'], ref_data['Pb'] * 1e6, color='w', edgecolors='orange', marker='o')
plt.scatter(ref_data['Zd'], ref_data['Pd'] * 1e6, color='orange', edgecolors='orange', marker='s')
plt.plot(z1, P_d, 'b',label='msrka_d')
plt.plot(z1, P_b, 'r',label='msrka_b')

flasher_srk = FlashVL(constants, correlations, liquid=liquid_srk, gas=gas_srk)
z1_srk, z2_srk, P_d_srk, P_b_srk = flasher_srk.plot_Pxy(T=T, pts=10, show=False, values=True, zmax=1)
plt.plot(z1_srk, P_d_srk, linestyle='--', color='blue',label='srk_d')
plt.plot(z1_srk, P_b_srk, linestyle='--', color='red',label='srk_b')

flasher_msrk = FlashVL(constants, correlations, liquid=liquid_msrk, gas=gas_msrk)
z1_msrk, z2_msrk, P_d_msrk, P_b_msrk = flasher_msrk.plot_Pxy(T=T, pts=10, show=False, values=True, zmax=1)
plt.plot(z1_msrk, P_d_msrk, linestyle=':', color='blue',label='msrk_d')
plt.plot(z1_msrk, P_b_msrk, linestyle=':', color='red',label='msrk_b')
# print(P_d_srk)
# plt.xlim(0,1)
# plt.ylim(1e6,10e6)
plt.legend()
plt.show()
