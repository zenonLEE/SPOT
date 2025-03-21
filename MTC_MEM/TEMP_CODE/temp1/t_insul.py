# def __init__(self, Do, Din, n, location):

import matplotlib.pyplot as plt
from thermo import ChemicalConstantsPackage

from reactor import vof, ks
import numpy as np
import pandas as pd
import scipy
from CoolProp.CoolProp import PropsSI
from prop_calculator import VLE, mixture_property, VLEThermo

R = 8.314


class Insulation:
    def __init__(self, Do, Din, n, location):

        # insulator parameters
        self.nit = n  # tube number of the insulator
        self.location = location  # 0 for hot in cold out; 1 for cold in hot out
        self.Do, self.Din = Do, Din
        self.thick = (self.Do - self.Din) / 2
        self.comp_list = ["CO2", "H2", "Methanol", "H2O", "CO"]

    @staticmethod
    def bi_diff(T, P):
        k = 1.38e-23
        mass = np.array([44, 2, 32, 18, 28])
        sigma = np.array([3.941, 2.827, 3.626, 2.641, 3.690])
        epsilon = np.array([195.2, 59.7, 481.8, 809.1, 91.7])
        sigma_mix, epsilon_mix = np.zeros((5, 5)), np.zeros((5, 5))
        mass_mix, D_bi = np.zeros((5, 5)), np.zeros((5, 5))
        for i in range(5):
            for j in range(5):
                sigma_mix[i, j] = (sigma[i] + sigma[j]) / 2
                epsilon_mix[i, j] = (epsilon[i] * epsilon[j]) ** 0.5
                mass_mix[i, j] = 2 / (1 / mass[i] + 1 / mass[j])
                T_star = k * T / (epsilon_mix[i, j] * k)
                omega = 1.06036 / T_star ** 0.1561 + 0.193 / np.exp(0.47635 * T_star) \
                        + 1.03587 / np.exp(1.52996 * T_star) + 1.76474 / np.exp(3.89411 * T_star)
                D_bi[i, j] = 1e-4 * 0.00266 * T ** 1.5 / (P * mass_mix[i, j] ** 0.5 * sigma_mix[i, j] ** 2 * omega)
        return D_bi

    @staticmethod
    def bi_md(T, P):
        Dij_523 = np.array([[4.36E-07, 5.81E-06],  # CH3OH-CO2 CH3OH-H2 1.5
                            [5.80E-07, 5.74E-06]])  # H2O-CO2 H2O-H2 1.5
        # Dij_523 = np.array([[4.25E-07, 5.81E-06],
        #                     [5.62E-7, 5.74E-06]])
        return Dij_523 * (T / 523) ** 1.5 * (70 / P)

    def ode_multi(self, x_in, inner_cond, outer_cond, P, properties):
        """
        ode for the concentration distribution along the channel, two condensates
        :param x_in: molar fraction at input; pd
        :param inner_cond: temperature, molar fraction, and radius at inside;list
        :param outer_cond: temperature, molar fraction, and radius at outside; list
        :param P: pressure of mixture, bar
        :param properties: heat capacity, diffusion coefficient, thermal conductivity of mixture; list
        :return: concentration and its slop
        """

        P = P * 1e5  # convert bar to pa
        x_main = x_in[["CO2", "H2"]] / np.sum(x_in[["CO2", "H2"]])
        # x_main = x_in[["CO", "H2"]] / np.sum(x_in[["CO", "H2"]])
        [cp_c, cp_d, k] = properties
        [T1, x_c1, x_d1, r1] = inner_cond
        [T2, x_c2, x_d2, r2] = outer_cond

        D_bi = self.bi_md((T1 + T2) / 2, P / 1e5)

        D_31, D_32 = D_bi[0, 0], D_bi[0, 1]
        D_41, D_42 = D_bi[1, 0], D_bi[1, 1]

        D_cm_id = 1 / (0.25 / D_31 + 0.75 / D_32)
        D_dm_id = 1 / (0.25 / D_41 + 0.75 / D_42)
        Nc_i = -(2 * P / R / (T2 + T1)) * D_cm_id * (x_c1 - x_c2) / (r1 - r2)
        Nd_i = -(2 * P / R / (T2 + T1)) * D_dm_id * (x_d1 - x_d2) / (r1 - r2)
        r = Nc_i / Nd_i

        diff_prop = [[x_c1, x_c2, cp_c, D_31, D_32, Nc_i],  # CH3OH
                     [x_d1, x_d2, cp_d, D_41, D_42, Nd_i]]  # H2O
        diff_prop = np.array(diff_prop)

        if x_c1 >= 3e-3 and x_d1 >= 1e-3:
            def model(z, y):
                [xc, xd, Nd, T, dTdz] = y
                Nc = r * Nd
                D_cm = -(xc * (Nc + Nd) - Nc) / (x_main["CO2"] * Nc / D_31 + x_main["H2"] * Nc / D_32)
                D_dm = -(xd * (Nc + Nd) - Nd) / (x_main["CO2"] * Nd / D_41 + x_main["H2"] * Nd / D_42)
                # D_cm = -(xc * (Nc + Nd) - Nc) / (x_main["CO"] * Nc / D_31 + x_main["H2"] * Nc / D_32)
                # D_dm = -(xd * (Nc + Nd) - Nd) / (x_main["CO"] * Nd / D_41 + x_main["H2"] * Nd / D_42)
                dxc_dz = -(Nc - xc * (Nc + Nd)) / (D_cm * (P / R / T))
                dxd_dz = -(Nd - xd * (Nc + Nd)) / (D_dm * (P / R / T))

                # dNc_dz = -Nc / z
                dNd_dz = -Nd / z
                d2T_dz2 = -dTdz / z + dTdz * cp_c * Nc / k + dTdz * cp_d * Nd / k
                return np.vstack((dxc_dz, dxd_dz, dNd_dz, dTdz, d2T_dz2))

            def bound(ya, yb):
                return np.array([ya[0] - x_c1, ya[1] - x_d1, ya[3] - T1,
                                 yb[0] - x_c2, yb[3] - T2])  # yb[1] - x_d2,

            xa, xb = r1, r2
            xini = np.linspace(xa, xb, 200)
            yini = np.zeros((5, xini.size))
            yini[0] = np.linspace(x_c1, x_c2, xini.size)
            yini[1] = np.linspace(x_d1, x_d2, xini.size)
            yini[2] = -(2 * P / R / (T2 + T1)) * D_dm_id * (x_d1 - x_d2) / (r1 - r2)
            yini[3] = np.linspace(T1, T2, xini.size)
            yini[4] = (T1 - T2) / (r1 - r2)
            res = scipy.integrate.solve_bvp(model, bound, xini, yini, tol=1e-8, max_nodes=5000)
            xsol = np.linspace(xa, xb, 200)
            ysol = res.sol(xsol)
            ysol = np.insert(ysol, 2, ysol[2] * r, axis=0)

        else:
            [diff_spec, stat_spec] = [0, 1] if x_c1 >= 3e-3 else [1, 0]
            [x_1, x_2, cp, D_1, D_2, N_id] = diff_prop[diff_spec]

            def model(z, y):
                [x, N, T, dTdz] = y
                # D_dm = (1 - x) / (x_main["CO2"] / D_1 + x_main["H2"] / D_2)
                D_dm = (1 - x) / (x_main["CO2"] / D_1 + x_main["H2"] / D_2)
                dxd_dz = -(N - x * N) / (D_dm * (P / R / T))

                dNd_dz = -N / z
                d2T_dz2 = -dTdz / z + dTdz * cp * N / k

                return np.vstack((dxd_dz, dNd_dz, dTdz, d2T_dz2))

            def bound(ya, yb):
                return np.array([ya[0] - x_1, ya[2] - T1,
                                 yb[0] - x_2, yb[2] - T2])

            xa, xb = r1, r2
            xini = np.linspace(xa, xb, 200)
            yini = np.zeros((4, xini.size))
            yini[0] = np.linspace(x_c1, x_c2, xini.size)
            yini[1] = N_id
            yini[2] = np.linspace(T1, T2, xini.size)
            yini[3] = (T1 - T2) / (r1 - r2)
            res = scipy.integrate.solve_bvp(model, bound, xini, yini, tol=1e-8, max_nodes=5000)
            xsol = np.linspace(xa, xb, 200)
            ysol = res.sol(xsol)
            # [xc, Nc, xd, Nd, T, dTdz]
            ysol = np.insert(ysol, stat_spec, diff_prop[stat_spec, 0], axis=0)
            ysol = np.insert(ysol, stat_spec + 2, 0, axis=0)

        return ysol

    def flux(self, Th, P, F_dict, Tc):
        """
        calculate the diffusional flux
        :param Tc: temperature of water in the cold side of reactor, K
        :param Th: temperature of gas in the reactor, K
        :param P: pressure of gas in the reactor, bar
        :param F_dict: gas component in the reactor, mol/s; ndarray
        :return: molar flux, heat flux, temperature variation and deviation of sim per length; dict
        """
        # calculate the correction to volumetric flow rate (m3/s)
        # calculate the partial pressure

        Ft = np.sum(F_dict)
        # insulator parameter
        radium = [self.Din / 2, self.Do / 2]
        # calculate the molar fraction of the mix
        xi_h = pd.Series(F_dict / Ft, index=self.comp_list)
        Pi_h = xi_h * P
        P_sat_CH3OH = PropsSI('Pp', 'T', Tc, 'Q', 1, "Methanol") * 1e-5
        P_sat_H2O = PropsSI('Pp', 'T', Tc, 'Q', 1, "H2O") * 1e-5
        P_sat = [P_sat_CH3OH, P_sat_H2O]  # [0,0]#
        condensates = ["Methanol", 'H2O']

        if xi_h["Methanol"] < 3e-3 and xi_h['H2O'] < 1e-3:
            # no reacted gas, end the calculation
            res = {
                "mflux": np.zeros(len(self.comp_list)), "hlg": 0,
                "hflux": 0, "Tvar": 0, "dev": 0
            }
            return res

        elif (xi_h["Methanol"] > 3e-3 and xi_h['H2O'] < 1e-3) or (xi_h["Methanol"] < 3e-3 and xi_h['H2O'] > 1e-3):
            # only one condensable spec is formed
            [diff_spec, stat_spec] = [0, 1] if xi_h["Methanol"] >= 3e-3 else [1, 0]

            if Pi_h[diff_spec + 2] < P_sat[diff_spec]:
                # condensation does not occur
                Tw = Th - 0.1
                mix_pro_ave = mixture_property((Tc + Tw) / 2, xi_h, P)  # rho is not used, hence z=1 is used
                k_e = mix_pro_ave["k"] * vof + ks * (1 - vof)  # effective heat conductivity of the insulator
                # heat conduction along the insulator
                qcv_cond = -2 * np.pi * k_e * (Tc - Tw) / np.log(radium[1 - self.location] / radium[self.location])
                property_h = mixture_property(Th, xi_h, P)  # rho is not used, hence z=1 is used
                dT = qcv_cond / Ft / property_h["cp_m"]
                # dT = -dT if self.location == 0 else dT
                na_CH3OH, na_H20, dev = 0, 0, 0
                qcv_lat = 0
            else:
                Pi_h_other_sum = P - Pi_h[diff_spec]
                Pi_c_other_sum = P - P_sat[diff_spec]
                Pi_c = Pi_h * (Pi_c_other_sum / Pi_h_other_sum)
                Pi_c[diff_spec + 2] = P_sat[diff_spec]
                xi_c = Pi_c / np.sum(Pi_c)

                Tw = Th - 0.1
                # gas properties inside the insulator
                property_c = mixture_property(Tc, xi_c, P)
                property_w = mixture_property(Tw, xi_h, P)
                mix_pro_ave = (property_w + property_c) / 2
                k_e = mix_pro_ave["k"] * vof + ks * (1 - vof)  # effective heat conductivity of the insulator
                # calculate the diffusional flux inside the insulator
                cold_cond = [Tc, xi_c["Methanol"], xi_c["H2O"], radium[1 - self.location]]
                hot_cond = [Th, xi_h["Methanol"], xi_h["H2O"], radium[self.location]]

                cond_list = [hot_cond, cold_cond]
                cal_property = [mix_pro_ave["cp_Methanol"], mix_pro_ave["cp_H2O"], k_e]
                # [xc, xd, Nc,Nd, T, dTdz]
                ode_res = self.ode_multi(xi_h, cond_list[self.location], cond_list[self.location - 1], P, cal_property)

                # mass flux inside the insulator, mol/(s m)
                na_H20 = ode_res[3][-self.location] * radium[-self.location] * 2 * np.pi * vof
                na_CH3OH = ode_res[2][-self.location] * radium[-self.location] * 2 * np.pi * vof

                # heat conduction inside the insulator, W/m
                qcv_cond = -k_e * ode_res[5][-self.location] * radium[-self.location] * 2 * np.pi
                property_h = mixture_property(Th, xi_h, P)  # rho is not used, hence z=1 is used
                dT = qcv_cond / Ft / property_h["cp_m"]
                dev = 0

                # calculate the latent heat
                qcv_lat = (PropsSI('H', 'T', Tc, 'Q', 1, 'water') -
                           PropsSI('H', 'T', Tc, 'Q', 0, 'water')) / 1000 * 18 * na_H20 + \
                          (PropsSI('H', 'T', Tc, 'Q', 1, 'Methanol') -
                           PropsSI('H', 'T', Tc, 'Q', 0, 'Methanol')) / 1000 * 32 * na_CH3OH

        else:
            # vle calculation, determine the dew pressure
            # vle_c = VLE(Tc, comp=self.comp_list)
            # P_dew = (1 / (xi_h["Methanol"] / P_sat_CH3OH + xi_h['H2O'] / P_sat_H2O))

            vle_c = VLEThermo(comp=self.comp_list)
            P_dew = vle_c.p_dew(Tc, xi_h.values)

            if P < P_dew:
                qcv_delta = 1e5
                Tw = Th - 0.1
                while qcv_delta > 20:
                    # condensation does not occur
                    # gas properties inside the insulator
                    mix_pro_ave = mixture_property((Tc + Tw) / 2, xi_h, P)  # rho is not used, hence z=1 is used
                    k_e = mix_pro_ave["k"] * vof + ks * (1 - vof)  # effective heat conductivity of the insulator

                    # heat conduction along the insulator
                    qcv_cond = -2 * np.pi * k_e * (Tc - Tw) / np.log(radium[1 - self.location] / radium[self.location])
                    # heat convection inside the reactor
                    # qcv_conv = -self.convection(Th, Pp, F_dict) * radium[self.location] * 2 * np.pi * (Tw - Th)
                    qcv_delta = 5  # abs(qcv_cond - qcv_conv)
                    Tw -= 1
                # temperature variation inside the reactor
                property_h = mixture_property(Th, xi_h, P)  # rho is not used, hence z=1 is used
                dT = qcv_cond / Ft / property_h["cp_m"]
                qcv_lant = 0
                # dT = -dT if self.location == 0 else dT
                na_CH3OH, na_H20, dev = 0, 0, 0
            else:
                # condensation occurs
                # calculate the composition of vapor and liquid phase

                # flash_comp, _ = vle_c.flash(Pp=Pp, mix=xi_h)
                # xi_c = flash_comp.loc['V']
                flash_g_comp, flash_l_comp = vle_c.flash(T=Tc, P=P, x=xi_h.values)
                r_ch3oh_h2o_fl = flash_l_comp[2] / flash_l_comp[3]
                xi_c = pd.Series(flash_g_comp, index=self.comp_list)
                re_ch3oh_h2o = 1
                while re_ch3oh_h2o > 0.05:
                    property_c = mixture_property(Tc, xi_c, P)
                    property_w = mixture_property(Th, xi_h, P)
                    mix_pro_ave = (property_w + property_c) / 2
                    k_e = mix_pro_ave["k"] * vof + ks * (1 - vof)  # effective heat conductivity of the insulator
                    # calculate the diffusional flux inside the insulator
                    cold_cond = [Tc, xi_c["Methanol"], xi_c["H2O"], radium[1 - self.location]]
                    hot_cond = [Th, xi_h["Methanol"], xi_h["H2O"], radium[self.location]]
                    cond_list = [hot_cond, cold_cond]
                    cal_property = [mix_pro_ave["cp_Methanol"], mix_pro_ave["cp_H2O"], k_e]
                    # [xc, xd, Nc,Nd, T, dTdz]
                    ode_res = self.ode_multi(xi_h, cond_list[self.location], cond_list[self.location - 1], P,
                                             cal_property)

                    # mass flux inside the insulator, mol/(s m)
                    na_H20 = ode_res[3][-self.location] * radium[-self.location] * 2 * np.pi * vof
                    na_CH3OH = ode_res[2][-self.location] * radium[-self.location] * 2 * np.pi * vof

                    # heat conduction inside the insulator, W/m
                    qcv_cond = -k_e * ode_res[5][-self.location] * radium[-self.location] * 2 * np.pi

                    gap_min = ode_res[1][-1] - cond_list[self.location - 1][2]
                    property_h = mixture_property(Th, xi_h, P)
                    dT = qcv_cond / Ft / property_h["cp_m"]  # k/m
                    dev = gap_min / xi_h["H2O"]

                    # update the comps of liquid
                    r_ch3oh_h2o = na_CH3OH / na_H20
                    re_ch3oh_h2o = abs(r_ch3oh_h2o_fl - r_ch3oh_h2o) / r_ch3oh_h2o_fl
                    re_ch3oh_h2o = 0.01
                    # print(re_ch3oh_h2o)

            # [xc, Nc, xd, Nd, T, dTdz]
            # xsol = np.linspace(0.02, 0.03, 200)
            # plt.plot(xsol, ode_res[0], label='H2O')
            # plt.plot(xsol, ode_res[2], label='CH3OH')
            # plt.legend()
            # plt.show()

            if dev > 0.1: print([na_CH3OH, na_H20], dev, 'dev too big')

        dF = np.zeros_like(F_dict)
        dF[2:4] = [na_CH3OH, na_H20]
        # print(k_e)
        if self.location == 0:
            dF = -1 * dF
            dT = -1 * dT
            if na_H20 < 0 or na_CH3OH < 0:
                print('err')
                print(cond_list, cal_property)
                print("*" * 10)
        else:
            if na_H20 > 0 or na_CH3OH > 0:
                print('err')
                print(cond_list, cal_property)
                print("*" * 10)

        # calculate the latent heat
        qcv_lat = (PropsSI('H', 'T', Tc, 'Q', 1, 'water') -
                   PropsSI('H', 'T', Tc, 'Q', 0, 'water')) / 1000 * 18 * na_H20 + \
                  (PropsSI('H', 'T', Tc, 'Q', 1, 'Methanol') -
                   PropsSI('H', 'T', Tc, 'Q', 0, 'Methanol')) / 1000 * 32 * na_CH3OH
        res = {
            "mflux": dF, "hflux": qcv_cond, "hlg": qcv_lat,
            "Tvar": dT, "dev": dev
        }
        return res  # dF, dT * h_phi, dev


Din, Dd = 0.05, 0.01
Do = Din + 2 * Dd
a = Insulation(Do, Din, 1, 0)

# flux(self, Th, Pp, F_dict, Tc)
Th = 503
P = 70
F = np.array([0.004211526, 0.013809958, 0.00153705, 0.001018258, 0.00058769])
Tc = 353

b = VLEThermo(["CO2", "H2", "Methanol", "H2O", "carbon monoxide"])
f = b.flash(Tc, P, x=F)
# PT = f.flash(T=Tc, Pp=Pp, zs=zs_l.tolist())
r_ch3oh_h2o = 0.5383897284155227 / 0.4409601246632461
print(r_ch3oh_h2o)
res = a.flux(Th, P, F, Tc)
mflux = res["mflux"]
print(mflux[2] / mflux[3])

subs = ["CO2", "H2", "Methanol", "H2O", "carbon monoxide"]
constants, correlations = ChemicalConstantsPackage.from_IDs(subs)
n_count = len(subs)

kijs_srk = np.array([[0, -0.3462, 0.0148, 0.0737, 0],
                     [-0.3462, 0, 0, 0, 0.0804],
                     [0.0148, 0, 0, -0.0789, 0],
                     [0.0737, 0, -0.0789, 0, 0],
                     [0, 0.0804, 0, 0, 0]])

#

eos_kwargs_srk = dict(Tcs=np.array(constants.Tcs).tolist(), Pcs=np.array(constants.Pcs).tolist(),
                      omegas=np.array(constants.omegas).tolist(), kijs=kijs_srk.tolist())

cp_cal = [HeatCapacityGas(CASRN=constants.CASs[i], MW=constants.MWs[i], method='TRCIG') for i in range(len(subs))]

liquid_srk = CEOSLiquid(SRKMIX, HeatCapacityGases=cp_cal, eos_kwargs=eos_kwargs_srk)
gas_srk = CEOSGas(SRKMIX, HeatCapacityGases=cp_cal, eos_kwargs=eos_kwargs_srk, T=T, P=P, zs=zs_l.tolist())

flasher_srk = FlashVL(constants, correlations, liquid=liquid_srk, gas=gas_srk)