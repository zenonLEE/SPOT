import numpy as np


def bi_diff_fu(T, P):
    mass = np.array([44, 2, 32, 18, 28])
    v = [26.9, 6.12, 29.9, 13.1, 18]  # [26.9, 7.07, 29.901, 12.7, 18.9]
    mass_mix, D_bi = np.zeros((5, 5)), np.zeros((5, 5))
    for i in range(5):
        for j in range(5):
            mass_mix[i, j] = 2 / (1 / mass[i] + 1 / mass[j])
            D_bi[i, j] = 1e-4 * 0.00143 * T ** 1.75 / (
                    P * mass_mix[i, j] ** 0.5 * (v[i] ** (1 / 3) + v[j] ** (1 / 3)) ** 2)
    return D_bi


T = (523 + 403) / 2
D = bi_diff_fu(T, 30)
print(D)
D_cm_id = 1 / (0.25 / D[2, 0] + 0.75 / D[2, 1])
print(D_cm_id)
