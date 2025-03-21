import matplotlib.pyplot as plt
from thermo import (ChemicalConstantsPackage, SRKMIX, FlashVL, CEOSLiquid, CEOSGas, HeatCapacityGas,
                    GibbsExcessLiquid, MSRKMIX, MSRKMaMIX, FlashPureVLS)
import numpy as np

subs = ["H2O"]  # ["CO2", "H2", "Methanol", "H2O", "carbon monoxide"]
n_count = len(subs)

T, P = 473, 7E6
Ts = np.arange(298.15, 599.15, 10)
G = np.zeros(len(Ts))
for i in range(n_count):
    for j in range(len(Ts)):
        constants, correlations = ChemicalConstantsPackage.from_IDs([subs[i]])
        cp_cal = HeatCapacityGas(CASRN=constants.CASs[0], MW=constants.MWs[0], method='TRCIG')
        eos_kwargs_srk = dict(Tcs=np.array(constants.Tcs), Pcs=np.array(constants.Pcs),
                              omegas=np.array(constants.omegas))
        liquid_srk = CEOSLiquid(SRKMIX, HeatCapacityGases=[cp_cal], eos_kwargs=eos_kwargs_srk, T=Ts[j], P=P)
        gas_srk = CEOSGas(SRKMIX, HeatCapacityGases=[cp_cal], eos_kwargs=eos_kwargs_srk, T=Ts[j], P=P)
        flasher_srk = FlashPureVLS(constants, correlations, liquids=[liquid_srk], gas=gas_srk, solids=[])
        PT = flasher_srk.flash(T=Ts[j], P=P)
        G[j] = PT.G() / 1E3
        # print(subs[i])
        # print(PT.G()/1E3)
plt.plot(Ts, G, 'r+')
plt.show()
