from prop_calculator import VLEThermo
import numpy as np


def mixer_real(F1, T1, F2, T2, P, species):
    """
    real mixer
    ref: Modelling, Estimation and Optimization of the Methanol Synthesis with Catalyst Deactivation
    :param F1: component of input gas 1, mol/s; ndarray
    :param T1: temperature of input gas 1, K
    :param F2: component of input gas 2, mol/s; ndarray
    :param T2: temperature of input gas 2, K
    :param P: pressure of input gas, bar
    :param species: component of input gas, list
    :return: molar flux of components, temperature
    """
    cal = VLEThermo(species)
    P = P * 1E5
    H_in = cal.cal_H(T1, P, F1) + cal.cal_H(T2, P, F2)
    F_out = F1 + F2
    H_diff = 100000
    if abs(T1 - T2) < 0.2:
        T_out = (T1 + T2) / 2
    else:
        for T in np.arange(min(T1 - 10, T2 - 10), max(T1 + 10, T2 + 10), 0.1):
            H_o = cal.cal_H(T, P, F_out)
            cal_diff = abs(H_o - H_in)
            if cal_diff < H_diff:
                H_diff = cal_diff
                T_out = T
            if cal_diff / H_in < 0.001 / 1e3:
                T_out = T
                break
    return F_out, T_out