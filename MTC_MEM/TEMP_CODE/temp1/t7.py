import thermo.heat_capacity
from scipy.constants import R
from thermo import (ChemicalConstantsPackage, HeatCapacityGas, CEOSGas,
                    SRKMIX)
import chemicals
import numpy as np
import matplotlib.pyplot as plt

comps = ["CO2"]
constants, properties = ChemicalConstantsPackage.from_IDs(comps)
# kijs = thermo.interaction_parameters.IPDB.get_ip_asymmetric_matrix('ChemSep PR', constants.CASs, 'kij')
kijs = np.array([[0, -0.3462, 0.0148, 0.0737, 0],
                 [-0.3462, 0, 0, 0, 0.0804],
                 [0.0148, 0, 0, -0.0789, 0],
                 [0.0737, 0, -0.0789, 0, 0],
                 [0, 0.0804, 0, 0, 0]])
# thermo.interaction_parameters.InteractionParameterDB.validate_table()
eos_kwargs = {'Pcs': constants.Pcs, 'Tcs': constants.Tcs, 'omegas': constants.omegas}

cp_coe = np.array([[1.9291E+01, 7.7562E-02, -6.8177E-05, 3.1746E-08, -6.1234E-12],
                   [2.6782E+01, 1.2238E-02, -2.2989E-05, 1.9404E-08, -5.2609E-12]]) / R
cp_cal = HeatCapacityGas(POLING_POLY=(50, 1000, cp_coe[0]))

sub_CAS = constants.CASs[0]
cp_cal2 = HeatCapacityGas(CASRN=sub_CAS, MW=constants.MWs[0], method='TRCIG')

Ts = np.arange(298.15, 598.15, 10)
H1, H2 = np.zeros(len(Ts)), np.zeros(len(Ts))
H3 = np.zeros((len(Ts)))
# chemicals.heat_capacity.Poling()
i = 0
for T in Ts:
    eos = CEOSGas(SRKMIX, eos_kwargs=eos_kwargs, HeatCapacityGases=[cp_cal2], zs=[1], T=T, P=1E6)
    # print(eos.H())
    H1[i] = cp_cal.T_dependent_property_integral(298.15, T)
    # H1[i] = chemicals.heat_capacity.Poling_integral(T, cp_coe[0],cp_coe[1],cp_coe[2],cp_coe[3],cp_coe[4])-\
    #         chemicals.heat_capacity.Poling_integral(298.15, cp_coe[0],cp_coe[1],cp_coe[2],cp_coe[3],cp_coe[4])
    H2[i] = cp_cal2.T_dependent_property_integral(298.15, T)
    H3[i] = eos.H()
    i += 1

plt.plot(Ts, H1, 'r')
plt.plot(Ts, H2, 'g')
plt.plot(Ts, H3, 'b')
plt.show()


