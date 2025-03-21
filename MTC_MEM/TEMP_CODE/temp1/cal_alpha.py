import numpy as np
from CoolProp.CoolProp import PropsSI
import matplotlib.pyplot as plt

D = 1E-6
ke = 0.2
Tr = 493
xb = 0.15
pt = 7e6
kes = np.arange(0.05, 0.1, 0.01).round(2)
Tcs = np.arange(313, 413, 1)
alpha = np.zeros((len(Tcs), len(kes)))
beta = np.zeros((len(Tcs), len(kes)))
Hr = 49e3
# j=0
# for Tc in Tcs:
#     T_ave = (Tr + Tc) / 2
#     ps1 = PropsSI('Pp', 'T', Tc, 'Q', 1, 'Methanol')
#     ps2 = PropsSI('Pp', 'T', Tc, 'Q', 1, 'H2O')
#     Hlg_1 = PropsSI('HMOLAR', 'T', Tc, 'Q', 1, 'Methanol') - PropsSI('HMOLAR', 'T', Tc, 'Q', 0, 'Methanol')
#     Hlg_2 = PropsSI('HMOLAR', 'T', Tc, 'Q', 1, 'H2O') - PropsSI('HMOLAR', 'T', Tc, 'Q', 0, 'H2O')
#     Hlg = (Hlg_1 + Hlg_2) / 2
#     beta[j] = Hr / Hlg
#     j += 1
#
# plt.plot(Tcs,beta)
# plt.show()
drive = np.zeros((len(Tcs), len(kes)))
i = 0
for Tc in Tcs:
    j = 0
    for ke in kes:
        T_ave = (Tr + Tc) / 2
        ps1 = PropsSI('Pp', 'T', Tc, 'Q', 1, 'Methanol')
        ps2 = PropsSI('Pp', 'T', Tc, 'Q', 1, 'H2O')
        Hlg_1 = PropsSI('HMOLAR', 'T', Tc, 'Q', 1, 'Methanol') - PropsSI('HMOLAR', 'T', Tc, 'Q', 0, 'Methanol')
        Hlg_2 = PropsSI('HMOLAR', 'T', Tc, 'Q', 1, 'H2O') - PropsSI('HMOLAR', 'T', Tc, 'Q', 0, 'H2O')
        Hlg = 50e3 #(Hlg_1+Hlg_2)/2
        ps = (ps1 + ps2) / 2
        xs = ps / pt
        x_ave = (xb - xs) / np.log((1 - xs) / (1 - xb))
        T_diff = Tr - Tc
        p_diff = pt * xb - ps
        alpha[i, j] = 1/((ke / D) * (T_diff / p_diff) * (8.314 * T_ave / Hlg) * x_ave)
        beta[i,j] = Hr/Hlg
        j += 1
    i += 1

for i, ke in enumerate(kes):
    plt.plot(Tcs, alpha[:,i], label=f'ke={ke}')
# plt.xlabel('Tc')
# plt.ylabel('alpha')
# plt.title('Alpha vs Tc for Different ke')
plt.legend(loc='upper left')
# plt.grid(True)
plt.show()
# plt.plot(Tcs, drive)
# plt.show()
