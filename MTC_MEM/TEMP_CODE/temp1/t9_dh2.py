from thermo import (ChemicalConstantsPackage, SRKMIX, FlashVL, CEOSLiquid, CEOSGas, HeatCapacityGas,
                    GibbsExcessLiquid, MSRKMIX, SRKMIXTranslated, MSRKMaMIX)
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

kijs_srk = np.array([[0, -0.0789],
                     [-0.0789, 0]])

kijs = np.array([[0, -0.075],
                 [-0.075, 0]])  # -0.075

p = [0.2359, 0.1277]
zra = [0.29056 - 0.08775 * constants.omegas[i] for i in range(n_count)]
# # print(zra)
c_vol = [-0.40768 * 8.314 * constants.Tcs[i] / constants.Pcs[i] * (0.29441 - zra[i]) for i in range(n_count)]
eos_kwargs_msrk = dict(Tcs=np.array(constants.Tcs), Pcs=np.array(constants.Pcs),
                       omegas=np.array(constants.omegas), kijs=kijs, S2s=np.array(p))
eos_kwargs_srk = dict(Tcs=np.array(constants.Tcs), Pcs=np.array(constants.Pcs),
                      omegas=np.array(constants.omegas), kijs=kijs_srk)  # , cs=-np.array(c_vol)

# thermo.HeatCapacityGas(CASRN=sub_CAS, MW=constants.MWs[0], method='TRCIG')
cp_cal = [HeatCapacityGas(CASRN=constants.CASs[i], MW=constants.MWs[i], method='TRCIG') for i in range(len(subs))]

P = 70e5
Ts = [473.15]  # np.arange(298.15, 458.15, 10)
j = 0
zs = np.array([0.8, 0.2])  # np.array([0, 0, 0.5, 0.5, 0])  # np.array([0.25, 0.75, 0, 0, 0])
zs_l = np.array([0.5, 0.5])

z_ch3oh = np.arange(0, 1.1, 0.1)

# H_msrk = np.zeros((len(Ts)))
# H_srk = np.zeros((len(Ts)))
# H_msrk_p = np.zeros((len(Ts), 2))
# H_srk_p = np.zeros((len(Ts), 2))

n_frac = len(z_ch3oh)
n_Ts = len(Ts)
H_msrk = np.zeros(n_frac)
H_srk = np.zeros(n_frac)

V_msrk = np.zeros(n_frac)
V_srk = np.zeros(n_frac)

# H_msrk_p = np.zeros((n_frac, 2))
# H_srk_p = np.zeros((n_frac, 2))
# T = 473.15
for i in range(len(z_ch3oh)):  # len(z_ch3oh)
    # gas = CEOSGas(MSRKMIX, HeatCapacityGases=cp_cal, eos_kwargs=eos_kwargs_msrk, T=Ts[i], Pp=Pp, zs=np.array([]))

    liquid = CEOSLiquid(MSRKMaMIX, HeatCapacityGases=cp_cal, eos_kwargs=eos_kwargs_msrk, T=Ts[j], P=P,
                        zs=np.array([z_ch3oh[i], 1 - z_ch3oh[i]]))
    liquid_srk = CEOSLiquid(SRKMIX, HeatCapacityGases=cp_cal, eos_kwargs=eos_kwargs_srk, T=Ts[j], P=P,
                            zs=np.array([z_ch3oh[i], 1 - z_ch3oh[i]]))
    # gas_srk = CEOSGas(SRKMIX, HeatCapacityGases=cp_cal, eos_kwargs=eos_kwargs_srk, T=Ts[i], Pp=Pp, zs=zs)

    ch3oh_msrk = CEOSLiquid(MSRKMaMIX, HeatCapacityGases=cp_cal, eos_kwargs=eos_kwargs_msrk, T=Ts[j], P=P,
                            zs=np.array([1, 0]))
    ch3oh_srk = CEOSLiquid(SRKMIX, HeatCapacityGases=cp_cal, eos_kwargs=eos_kwargs_srk, T=Ts[j], P=P,
                           zs=np.array([1, 0]))

    h2o_msrk = CEOSLiquid(MSRKMaMIX, HeatCapacityGases=cp_cal, eos_kwargs=eos_kwargs_msrk, T=Ts[j], P=P,
                          zs=np.array([0, 1]))
    h2o_srk = CEOSLiquid(SRKMIX, HeatCapacityGases=cp_cal, eos_kwargs=eos_kwargs_srk, T=Ts[j], P=P, zs=np.array([0, 1]))

    H_msrk[i] = liquid.H() - (ch3oh_msrk.H() * z_ch3oh[i] + h2o_msrk.H() * (1 - z_ch3oh[i]))
    H_srk[i] = liquid_srk.H() - (ch3oh_srk.H() * z_ch3oh[i] + h2o_srk.H() * (1 - z_ch3oh[i]))

    V_msrk[i] = liquid.V() - (ch3oh_msrk.V() * z_ch3oh[i] + h2o_msrk.V() * (1 - z_ch3oh[i]))
    V_srk[i] = liquid_srk.V() - (ch3oh_srk.V() * z_ch3oh[i] + h2o_srk.V() * (1 - z_ch3oh[i]))




# flasher.d


ref_path = r"D:\document\00Study\05多联产系统\甲醇单元\Gibbs\热力学性质_ASPEN.xlsx"
ref_data = pd.read_excel(ref_path, sheet_name='SIMO')  # CO2_H2  CH3OH_H2O
ref_data_p = ref_data[
    (P / 1E6 * 0.95 <= ref_data['Pp']) & (ref_data['Pp'] <= P / 1E6 * 1.05) &
    (Ts[j] * 0.95 <= ref_data['T']) & (ref_data['T'] <= Ts[j] * 1.05)
    ]  #

# plt.plot(z_ch3oh, V_msrk*1E6, 'b', label='msrk')
# plt.plot(z_ch3oh, V_srk*1E6, 'r', label='srk')
# plt.scatter(ref_data_p['X'], -ref_data_p['V']/100, color='w', edgecolors='orange', marker='o')

plt.plot(z_ch3oh, H_msrk, 'b', label='msrk')
plt.plot(z_ch3oh, H_srk, 'r', label='srk')
plt.scatter(ref_data_p['X'], -ref_data_p['H'], color='w', edgecolors='orange', marker='o')
# plt.plot(z_ch3oh, H_nrtl*1e6, 'g')



# plt.xlim(0,1)
# plt.ylim(1e6,10e6)
plt.legend()
plt.show()
