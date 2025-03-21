from thermo import (ChemicalConstantsPackage, SRKMIX, FlashVL, CEOSLiquid, CEOSGas, HeatCapacityGas,
                    GibbsExcessLiquid, MSRKMIX, MSRKMaMIX)
import numpy as np
from insulator import Insulation

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
p = [0, 0, 0.2359, 0.1277, 0]
cp_cal = [HeatCapacityGas(CASRN=constants.CASs[i], MW=constants.MWs[i], method='TRCIG') for i in range(len(subs))]
comps = [0, 0, 0.5, 0.5, 0]  # [0.25,0.75,0,0,0]#[0.2, 0.55, 0.1, 0.1, 0.05]
te = np.array([0, 0, 1, 1, 0])  # np.array([0.5461247, 1.65821958, 0.44395256, 0.4538753, 0.00992274])
print(np.sum(te))
eos_kwargs_msrk = dict(Tcs=np.array(constants.Tcs), Pcs=np.array(constants.Pcs),
                       omegas=np.array(constants.omegas), kijs=kijs_msrk.tolist(), S2s=np.array(p),
                       HeatCapacityGas=cp_cal)
eos_kwargs_srk = dict(Tcs=np.array(constants.Tcs), Pcs=np.array(constants.Pcs),
                      omegas=np.array(constants.omegas), kijs=kijs_srk, HeatCapacityGas=cp_cal)

T, P = 473, 7E6

liquid_srk = CEOSLiquid(SRKMIX, eos_kwargs=eos_kwargs_srk, zs=np.array(te), T=T, P=P)
gas_srk = CEOSGas(SRKMIX, eos_kwargs=eos_kwargs_srk, zs=np.array(te), T=T, P=P)
flasher_srk = FlashVL(constants, correlations, liquid=liquid_srk, gas=gas_srk)
PT = flasher_srk.flash(zs=te / np.sum(te), T=T, P=P)
print(PT.G() * np.sum(te))

x = [0.43907971, 1.36908887, 0.27656504, 0.20172921, 0.02592486]
