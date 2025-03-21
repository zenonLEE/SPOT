from thermo import (ChemicalConstantsPackage, SRKMIX, FlashVL, CEOSLiquid, CEOSGas, HeatCapacityGas,
                    GibbsExcessLiquid, MSRKMIX)
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

# kijs = np.array([[0, -0.3462, 0.0148, 0.0737, 0],
#                  [-0.3462, 0, 0, 0, 0.0804],
#                  [0.0148, 0, 0, -0.0789, 0],
#                  [0.0737, 0, -0.0789, 0, 0],
#                  [0, 0.0804, 0, 0, 0]])

kijs = np.array([[0, 0.1164, 0.1, 0.3, 0.1164],
                 [0.1164, 0, -0.125, -0.745, -0.0007],
                 [0.1, -0.125, 0, -0.075, -0.37],
                 [0.3, -0.745, -0.075, 0, -0.474],
                 [0.1164, -0.0007, -0.37, -0.474, 0]])

p = [0, 0, 0.2359, 0.1277, 0]

eos_kwargs_msrk = dict(Tcs=np.array(constants.Tcs), Pcs=np.array(constants.Pcs),
                       omegas=np.array(constants.omegas), kijs=kijs, S2s=p)
eos_kwargs_srk = dict(Tcs=np.array(constants.Tcs), Pcs=np.array(constants.Pcs),
                       omegas=np.array(constants.omegas), kijs=kijs)

# thermo.HeatCapacityGas(CASRN=sub_CAS, MW=constants.MWs[0], method='TRCIG')
cp_cal = [HeatCapacityGas(CASRN=constants.CASs[i], MW=constants.MWs[i], method='TRCIG') for i in range(len(subs))]

P = 5e6
Ts = np.arange(298.15, 598.15, 10)
zs = np.array([0, 0, 0.5, 0.5, 0])  # np.array([0, 0, 0.5, 0.5, 0])  # np.array([0.25, 0.75, 0, 0, 0])
zs_l = np.array([0, 0, 0.5, 0.5, 0])

# paras for NRTL
alpha_D = tau_E = tau_F = tau_G = tau_H = np.zeros((n_count, n_count))
tau_A = np.zeros((n_count, n_count))
tau_B = np.zeros((n_count, n_count))
alpha_C = np.ones((n_count, n_count)) * 0.3
tau_A[3, 2], tau_A[2, 3] = 2.7322, -0.693
tau_B[3, 2], tau_B[2, 3] = -617.269, 172.987
nrtl_kwargs = (tau_A, tau_B, tau_E, tau_F, tau_G, tau_H, alpha_C, alpha_D)

# liquid = CEOSLiquid(SRKMIX, HeatCapacityGases=cp_cal, eos_kwargs=eos_kwargs)
H_l, G, S = np.zeros(len(Ts)), np.zeros(len(Ts)), np.zeros(len(Ts))
H_g = np.zeros(len(Ts))
G2 = np.zeros(len(Ts))
H2_l = np.zeros(len(Ts))
H_srk_l = np.zeros(len(Ts))
i = 0
for T in Ts:
    gas = CEOSGas(MSRKMIX, HeatCapacityGases=cp_cal, eos_kwargs=eos_kwargs_msrk, T=T, P=P, zs=zs)
    liquid = CEOSLiquid(MSRKMIX, HeatCapacityGases=cp_cal, eos_kwargs=eos_kwargs_msrk, T=T, P=P, zs=zs)
    liquid_srk = CEOSLiquid(SRKMIX, HeatCapacityGases=cp_cal, eos_kwargs=eos_kwargs_srk, T=T, P=P, zs=zs)
    # liquid = CEOSLiquid(GibbsExcessLiquid, HeatCapacityGases=cp_cal,T=T, Pp=Pp, zs=zs,
    #                     VaporPressures=)
    liquid_nrtl = NRTL(T, xs=zs_l, ABEFGHCD=nrtl_kwargs)
    liquid_nrtl.HE()

    liquid2 = GibbsExcessLiquid(VaporPressures=correlations.VaporPressures,
                                HeatCapacityGases=correlations.HeatCapacityGases,
                                VolumeLiquids=correlations.VolumeLiquids,
                                GibbsExcessModel=liquid_nrtl,
                                T=T, P=P, zs=zs)  # equilibrium_basis='PhiSat', caloric_basis='PhiSat',
    # liquid.H()
    # print(liquid.GE())
    # gas.Tsat
    # print(liquid.HE())
    H_l[i] = liquid.H()
    H2_l[i] = liquid2.H()
    H_srk_l[i] = liquid_srk.H()
    # cal_HMIG(subs, T, Pp, zs_l) * 1000  # + liquid.HE() # liquid.H_from_phi() if T < 512 else
    H_g[i] = gas.H()
    # gas.H_from_phi()
    G[i] = liquid.G()  # + liquid.GE()  # liquid.G() if T < 512 else gas.G()

    S[i] = liquid.S()  # + liquid.SE()  # gas.S()
    i += 1
H_ref = np.dot(np.array(constants.Hfgs), zs)
S_ref = np.dot(np.array(constants.Sfgs), zs)  # thermo has considered the terms x*lnx
H_l = (H_l + H_ref) / 1e3
H_g = (H_g + H_ref) / 1e3
H2_l = (H2_l + H_ref) / 1e3
H_srk_l = (H_srk_l+H_ref)/1E3
S = (S + S_ref) / 1e3
G = (G + H_ref - Ts * S_ref) / 1e3
# G1 = H - Ts * S

ref_path = r"D:\document\00Study\05多联产系统\甲醇单元\Gibbs\热力学性质_ASPEN.xlsx"
ref_data = pd.read_excel(ref_path, index_col='T', sheet_name='CH3OH_H2O')  # CO2_H2  CH3OH_H2O

plt.plot(Ts, H_l, 'r+')
# plt.plot(Ts, H_g, 'r*')
# plt.plot(Ts, H2_l, 'g')
plt.plot(Ts, H_srk_l, 'y')
# plt.plot(Ts, G1, 'r')
plt.plot(ref_data.index, ref_data['H'], 'b')
# plt.plot(ref_data.index, ref_data['GV'], 'b')

plt.show()
