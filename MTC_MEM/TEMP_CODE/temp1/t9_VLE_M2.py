from thermo import (ChemicalConstantsPackage, SRKMIX, FlashVL, CEOSLiquid, CEOSGas, HeatCapacityGas,
                    GibbsExcessLiquid, MSRKMIX, MSRKMaMIX)
import numpy as np

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

eos_kwargs_msrk = dict(Tcs=np.array(constants.Tcs), Pcs=np.array(constants.Pcs),
                       omegas=np.array(constants.omegas), kijs=kijs_msrk.tolist(), S2s=np.array(p))
eos_kwargs_srk = dict(Tcs=np.array(constants.Tcs), Pcs=np.array(constants.Pcs),
                      omegas=np.array(constants.omegas), kijs=kijs_srk)

# thermo.HeatCapacityGas(CASRN=sub_CAS, MW=constants.MWs[0], method='TRCIG')
cp_cal = [HeatCapacityGas(CASRN=constants.CASs[i], MW=constants.MWs[i], method='TRCIG') for i in range(len(subs))]

# Pp = 5e6
Ts = np.arange(298.15, 598.15, 10)
zs = np.array([0.5, 0.5])  # np.array([0, 0, 0.5, 0.5, 0])  # np.array([0.25, 0.75, 0, 0, 0])
zs_l = np.array([0.2, 0.55, 0.1, 0.1, 0.05])

T = 453
P = 7E6

liquid_srk = CEOSLiquid(SRKMIX, HeatCapacityGases=cp_cal, eos_kwargs=eos_kwargs_srk)
gas_srk = CEOSGas(SRKMIX, HeatCapacityGases=cp_cal, eos_kwargs=eos_kwargs_srk)

flasher_srk = FlashVL(constants, correlations, liquid=liquid_srk, gas=gas_srk)
flasher_srk.flash(zs=zs_l)
PT = flasher_srk.flash(T=T, P=P, zs=zs_l)

print(PT.VF)
print(PT.H(), PT.G())

g_z = np.array(PT.gas.zs)
l_z = np.array(PT.liquid0.zs)

gas_f = CEOSGas(SRKMIX, HeatCapacityGases=cp_cal, eos_kwargs=eos_kwargs_srk, zs=g_z, T=T, P=P)
liquid_f = CEOSLiquid(SRKMIX, HeatCapacityGases=cp_cal, eos_kwargs=eos_kwargs_srk, zs=l_z, T=T, P=P)
print(l_z)
print(gas_f.H(), gas_f.G())
print(liquid_f.H(), liquid_f.G())
print(gas_f.H() * PT.VF + liquid_f.H() * (1 - PT.VF))

print("*"*10)

liquid_msrk = CEOSLiquid(MSRKMIX, HeatCapacityGases=cp_cal, eos_kwargs=eos_kwargs_msrk)
gas_msrk = CEOSGas(MSRKMIX, HeatCapacityGases=cp_cal, eos_kwargs=eos_kwargs_msrk)
flasher_msrk = FlashVL(constants, correlations, liquid=liquid_msrk, gas=gas_msrk)
PT = flasher_msrk.flash(T=T, P=P, zs=zs_l)

print(PT.VF)
print(PT.H(), PT.G())

g_z = np.array(PT.gas.zs)
l_z = np.array(PT.liquid0.zs)
print(l_z)

gas_f = CEOSGas(MSRKMIX, HeatCapacityGases=cp_cal, eos_kwargs=eos_kwargs_msrk, zs=g_z, T=T, P=P)
liquid_f = CEOSLiquid(MSRKMIX, HeatCapacityGases=cp_cal, eos_kwargs=eos_kwargs_msrk, zs=l_z, T=T, P=P)
print(gas_f.H(), gas_f.G())
print(liquid_f.H(), liquid_f.G())
print(gas_f.H() * PT.VF + liquid_f.H() * (1 - PT.VF))
print(gas_f.G() * PT.VF + liquid_f.G() * (1 - PT.VF))