import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import thermo
from fluids.constants import R


def cal_HIG(sub, T):
    constants, correlations = thermo.ChemicalConstantsPackage.from_IDs(IDs=[sub])
    sub_CAS = constants.CASs[0]
    Hfs = constants.Hfgs[0] / 1000
    cp_cal = thermo.HeatCapacityGas(CASRN=sub_CAS, MW=constants.MWs[0], method='TRCIG')
    Hig = cp_cal.T_dependent_property_integral(298.15, T) / 1000 + Hfs
    return Hig



def cal_H(sub, T, P):
    constants, correlations = thermo.ChemicalConstantsPackage.from_IDs(IDs=[sub])
    sub_CAS = constants.CASs[0]
    Hfs = constants.Hfgs[0] / 1000
    Gfs = constants.Gfgs[0] / 1000
    Sfs = constants.Sfgs[0] / 1000

    # calculate the properties for ideal gas
    cp_cal = thermo.HeatCapacityGas(CASRN=sub_CAS, MW=constants.MWs[0], method='TRCIG')
    Hig = cp_cal.T_dependent_property_integral(298.15, T) / 1000 + Hfs
    Sig = cp_cal.T_dependent_property_integral_over_T(298.15, T) / 1000 + Sfs
    Gig = Hig - T * Sig
    zra = 0.29056 - 0.08775 * constants.omegas[0]
    c_vol = -0.40768 * 8.314 * (constants.Tcs[0]) / (constants.Pcs[0]) * (0.29441 - zra)

    # calculate the departure enthalpy
    eos = thermo.SRK(Tc=constants.Tcs[0], Pc=constants.Pcs[0], omega=constants.omegas[0],
                     T=T, P=P)

    Tsat = eos.Tsat(P, True)

    try:
        h_dep = eos.H_dep_l if T < Tsat else eos.H_dep_g  # eos.H_dep_g - eos.Hvap(T)
        s_dep = eos.S_dep_l if T < Tsat else eos.S_dep_g
    except AttributeError:
        h_dep = eos.H_dep_l
        s_dep = eos.S_dep_l

    H = h_dep / 1000 + Hig
    S = (s_dep - R * np.log(P / 1e5)) / 1000 + Sig
    G = H - T * S
