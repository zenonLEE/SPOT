import numpy as np

import thermo
from fluids.constants import R
from t5 import *
import matplotlib.pyplot as plt
import pandas as pd

comps = ["CO2", "H2", "Methanol", "H2O", "carbon monoxide"]
constants, properties = thermo.ChemicalConstantsPackage.from_IDs(comps)
# kijs = thermo.interaction_parameters.IPDB.get_ip_asymmetric_matrix('ChemSep PR', constants.CASs, 'kij')
kijs = np.array([[0, -0.3462, 0.0148, 0.0737, 0],
                 [-0.3462, 0, 0, 0, 0.0804],
                 [0.0148, 0, 0, -0.0789, 0],
                 [0.0737, 0, -0.0789, 0, 0],
                 [0, 0.0804, 0, 0, 0]])
# thermo.interaction_parameters.InteractionParameterDB.validate_table()
eos_kwargs = {'Pcs': constants.Pcs, 'Tcs': constants.Tcs, 'omegas': constants.omegas, 'kijs': kijs}

P = 5E6
Ts = np.arange(298.15, 599.15, 10)
gas_frac = [0.25, 0.75, 0, 0, 0]
liq_frac = [0, 0, 0.5, 0.5, 0]
frac = gas_frac
H, Hig = np.zeros(len(Ts)), np.zeros(len(Ts))
S, Sig = np.zeros(len(Ts)), np.zeros(len(Ts))
G, GIG = np.zeros(len(Ts)), np.zeros(len(Ts))
i = 0
for T in Ts:
    eos = thermo.eos_mix.SRKMIX(T=T, P=P, Tcs=constants.Tcs, Pcs=constants.Pcs, omegas=constants.omegas,
                                zs=frac, kijs=kijs)
    try:
        fi = eos.fugacities_g
    except AttributeError:
        fi = eos.fugacities_l

    gamma_i = np.array(fi) / cal_FiIG(T, P, comps, frac)
    gamma_i = np.where(gamma_i == 0, 1, gamma_i)
    # print(gamma_i)
    # print(gamma_i)
    Hig[i] = cal_HMIG(comps, T, P, frac)
    Sig[i] = cal_SMIG(comps, T, P, frac)
    H[i] = Hig[i]
    S[i] = Sig[i]

    GIG[i] = Hig[i] - T * Sig[i]
    G[i] = GIG[i] + np.sum(R * T * np.log(gamma_i) * frac)/1000

    i += 1

# G = H - Ts*S
ref_path = r"D:\document\00Study\05多联产系统\甲醇单元\Gibbs\热力学性质_ASPEN.xlsx"
ref_data = pd.read_excel(ref_path, index_col='T', sheet_name='CO2_H2')

plt.plot(Ts, H, 'r+')
# plt.plot(Ts, GIG, 'r')
plt.plot(ref_data.index, ref_data['H'], 'b+')
plt.show()
