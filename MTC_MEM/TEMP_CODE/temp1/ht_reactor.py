import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
from prop_calculator import VLEThermo
from reactor import Reaction

# comp = ["CO2", "H2", "Methanol", "H2O", "CO", 'Ar']
# a = VLEThermo(comp)
# x = np.array([3,82,0,0,4,11])/100
# z= a.z(493.2, 50,x)
# print(z)
# def __init__(self, L, D, Dc, n, phi, rho, chem_para, T0, P0, F0, eos, qmh=0):
L = 1
Dc = 0.02
n = 1
phi = 0.67
rho = 1950
T0 = 493
P0 = 70
Dt = 0.04
F0 = np.array([0.008154456, 0.024463369, 0, 0, 0, 0])
chem_info = {'stoichiometry': {'1': [-1, -3, 1, 1, 0, 0],
                               '2': [-1, -1, 0, 1, 1, 0]},
             'kad': {'H2': [0.499, 17197], 'H2O': [6.62e-11, 124119], 'H2O/H2': [3453.38, 0]},
             'kr': {'1': [1.07, 36696], '2': [12200000000.0, -94765]},
             'keq': {'1': [3066, -10.592], '2': [-2073, 2.029]},
             'heat_reaction': {'1': [35890.0, 40047000.0], '2': [9177.0, -44325000.0]},
             'kn_model': 'BU'}

a = Reaction(L, D=Dt, Dc=Dc, n=n, phi=phi, rho=rho, chem_para=chem_info, T0=T0, P0=P0, F0=F0, eos=1, qmh=1)
T = 513
P = 70

F_ns = np.arange(1, 10, 2).tolist() + np.around(np.arange(0.1, 0.9, 0.3), 1).tolist()
# 0.012292527328704943
v = 0.012292527328704943

Us = []  # np.array(F_ns)*v
for F_n in F_ns:
    Us.append(1 / a.htr(T=T, P=P, F_dict=F0 * F_n))
plt.scatter(F_ns, Us)

plt.grid()
plt.show()
# Ut =
