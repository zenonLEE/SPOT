import time
import warnings
import numpy as np
import pandas as pd
import win32com.client
from CoolProp.CoolProp import PropsSI
import scipy.optimize as opt

warnings.filterwarnings('ignore')
from thermo import (ChemicalConstantsPackage, SRKMIX, FlashVL, CEOSLiquid, CEOSGas, HeatCapacityGas,
                    FlashVLN, heat_capacity)

# from thermo.unifac import DOUFSG, DOUFIP2016, UNIFAC

# T_gas = 200  # 310.93  # K
# P_gas = 30  # 6.3  # bar
R = 8.314


def bi_diff(T, P):
    k = 1.38e-23
    mass = np.array([44.01, 2.016, 32.042, 18.015, 28.01])
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
            D_bi[i, j] = 1e-4 * 0.00266 * T ** 1.5 / (
                    P * mass_mix[i, j] ** 0.5 * sigma_mix[i, j] ** 2 * omega)
    return D_bi


def D_fu(T, P):
    """
    Fuller method to calculate diffusion coe
    :param T: K
    :param P: bar
    :return: m2/s
    """
    sub = ['H2', 'CO', 'CO2', 'H2O', 'MEOH']
    vf = np.array([6.12, 18, 26.7, 13.1, 29.9])
    # np.array([26.9, 7.07, 29.9, 12.7, 18.9])  # [26.7, 6.12, 29.9,13.1, 18]
    # [26.9, 6.12, 31.25, 13.1 ,18] in comsol
    mass = np.array([2.016, 28.01, 44.01, 18.015, 32.042])  # np.array([44.01, 2.016, 32.042, 18.015, 28.01])
    num = len(vf)
    res = np.zeros((num, num))
    for i in range(num):
        for j in range(num):
            mass_mix = 2 / (1 / mass[i] + 1 / mass[j])
            res[i, j] = 1e-4 * 0.00143 * T ** 1.75 / (P * mass_mix ** 0.5 * (vf[i] ** (1 / 3) + vf[j] ** (1 / 3)) ** 2)
    return pd.DataFrame(res, index=sub, columns=sub)


def Dm_fu(stream_source):
    """
    calculate dif coe in mix using Fuller method
    :param stream_source:
    :return:
    """
    D_bi = D_fu(stream_source['T'], stream_source['P'])
    sub = ['H2', 'CO', 'CO2', 'H2O', 'MEOH']
    xi = stream_source[sub] / stream_source[sub].sum()
    Dm = pd.Series(index=sub)
    for i in sub:
        temp = 0
        for j in sub:
            if j != i:
                temp += xi[j] / D_bi.loc[i, j]
        Dm[i] = 1 / temp
    return Dm


def trans_prop(T, stream):
    """
    calculate the properties of gas mixture
    :param T: gas temperature, K
    :param xi_gas: molar fraction; pd.Serize
    :param Pt: total pressure, bar
    :param z: compression factor
    :param rho_only: compute rho only
    :param phis: fugacity coes
    :return: thermal conductivity W/(m K), viscosity Pa s, heat capacity J/mol/K; pd.series
    """
    # prepare data for calculation

    new_index = ['Methanol' if i == 'MEOH' else i for i in stream.index.tolist()]
    stream.index = new_index
    n = len(stream.index) - 2  # number of gas species
    xi_gas = stream.drop['T', 'P'] / stream.drop['T', 'P'].sum()
    [cp, k, vis, M, rho] = np.ones((5, n)) * 1e-5

    Ti_sat = pd.Series(np.ones(n) * 100, index=xi_gas.index)

    if 'Methanol' in xi_gas.index:
        try:
            Ti_sat['Methanol'] = PropsSI('T', 'P', fi_gas['Methanol'], 'Q', 1, 'Methanol')
        except ValueError:
            Ti_sat['Methanol'] = 300
    if "H2O" in xi_gas.index:
        try:
            Ti_sat['H2O'] = PropsSI('T', 'P', fi_gas['H2O'], 'Q', 1, 'H2O')
        except ValueError:
            Ti_sat['H2O'] = 300
    i = 0
    for comp in xi_gas.index:
        M[i] = PropsSI('MOLARMASS', 'T', T, 'P', 1e5, comp)  # molar weight, kg/mol
        i += 1

    i = 0
    # calculate the properties of pure gases
    for comp in xi_gas.index:
        gas = "N2" if comp == "CO" else comp  # "CO" is not available in CoolProp
        if xi_gas[comp] > 0:
            if T > Ti_sat[comp] * 1.01:
                # thermal conductivity, W/(m K)
                k[i] = PropsSI('L', 'T', T, 'P', stream['P'], gas)
                # viscosity, Pa S
                vis[i] = PropsSI('V', 'T', T, 'P', stream['P'], gas)
                # heat capacity, J/(mol K)
                cp[i] = PropsSI('CPMOLAR', 'T', T, 'P', stream['P'], gas)
                # density, kg/m3
            else:
                cp[i] = PropsSI('CPMOLAR', 'T', T, 'Q', 1, gas)
                k[i] = PropsSI('L', 'T', T, 'Q', 1, gas)
                vis[i] = PropsSI('V', 'T', T, 'Q', 1, gas)
        else:
            # thermal conductivity, W/(m K)
            k[i] = 0
            # viscosity, Pa S
            vis[i] = 1e-10
            # heat capacity, J/(mol K)
            cp[i] = 0
            # density, kg/m3
            rho[i] = 0
        i += 1
    # calculate the properties of mixture
    phi, denominator = np.ones((n, n)), np.ones((n, n))  # Wilke coefficient
    vis_m, k_m = 0, 0
    for i in range(n):
        for j in np.arange(n):
            phi[i, j] = (1 + (vis[i] / vis[j]) ** 0.5 * (M[j] / M[i]) ** 0.25) ** 2 / (8 * (1 + M[i] / M[j])) ** 0.5
            denominator[i, j] = xi_gas[j] * phi[i, j]  # if i != j else 0
        vis_m += xi_gas[i] * vis[i] / np.sum(denominator[i])

    return pd.Series([k_m, vis_m], index=["k", "vis"])


def mixture_property(T, xi_gas, Pt, z=1, rho_only=False, phis=None):
    """
    calculate the properties of gas mixture
    :param T: gas temperature, K
    :param xi_gas: molar fraction; pd.Serize
    :param Pt: total pressure, bar
    :param z: compression factor
    :param rho_only: compute rho only
    :param phis: fugacity coes
    :return: thermal conductivity W/(m K), viscosity Pa s, heat capacity J/mol/K; pd.series
    """
    # prepare data for calculation
    if phis is not None:
        phi_cal = VLEThermo(xi_gas.index.tolist())
        phi_gas = phi_cal.phi(T, Pt, xi_gas.values)
    else:
        phi_gas = np.ones(len(xi_gas))
    new_index = ['CO' if i == 'carbon monoxide' else i for i in xi_gas.index.tolist()]
    xi_gas.index = new_index
    n = len(xi_gas.index)  # number of gas species
    xi_gas = xi_gas / xi_gas.sum()
    [cp, k, vis, M, rho] = np.ones((5, n)) * 1e-5

    fi_gas = xi_gas * Pt * 1e5 * phi_gas  # convert bar to pa
    Ti_sat = pd.Series(np.ones(n) * 100, index=xi_gas.index)

    if 'Methanol' in xi_gas.index:
        try:
            Ti_sat['Methanol'] = PropsSI('T', 'P', fi_gas['Methanol'], 'Q', 1, 'Methanol')
        except ValueError:
            Ti_sat['Methanol'] = 300
    if "H2O" in xi_gas.index:
        try:
            Ti_sat['H2O'] = PropsSI('T', 'P', fi_gas['H2O'], 'Q', 1, 'H2O')
        except ValueError:
            Ti_sat['H2O'] = 300
    i = 0
    for comp in xi_gas.index:
        M[i] = PropsSI('MOLARMASS', 'T', T, 'P', 1e5, comp)  # molar weight, kg/mol
        i += 1
    # print(Ti_sat)
    # print(fi_gas)
    M_m = np.sum(M * xi_gas)  # molar weight of mixture, kg/mol
    # print(M_m)
    rho_m = Pt * 1E5 * M_m / (z * R * T)  # kg/m3 np.sum(rho)
    # print(rho_m/M_m)
    if rho_only:
        return pd.Series([0, 0, rho_m, cp[2], cp[3], 0],
                         index=["k", "vis", 'rho', 'cp_' + xi_gas.index[2], 'cp_' + xi_gas.index[3], "cp_m"])

    i = 0
    # calculate the properties of pure gases
    for comp in xi_gas.index:
        gas = "N2" if comp == "CO" else comp  # "CO" is not available in CoolProp
        if fi_gas[comp] > 1000:
            if T > Ti_sat[comp] * 1.01:
                # thermal conductivity, W/(m K)
                k[i] = PropsSI('L', 'T', T, 'P', Pt, gas)
                # viscosity, Pa S
                vis[i] = PropsSI('V', 'T', T, 'P', Pt, gas)
                # heat capacity, J/(mol K)
                cp[i] = PropsSI('CPMOLAR', 'T', T, 'P', Pt, gas)
                # density, kg/m3
                # rho[i] = PropsSI('D', 'T', T, 'Pp', xi_gas[comp], gas)
            else:
                cp[i] = PropsSI('CPMOLAR', 'T', T, 'Q', 1, gas)
                k[i] = PropsSI('L', 'T', T, 'Q', 1, gas)
                vis[i] = PropsSI('V', 'T', T, 'Q', 1, gas)
                # print(comp, T, cp[i])
                # rho[i] = PropsSI('D', 'T', T, 'Q', 1, gas)
        else:
            # thermal conductivity, W/(m K)
            k[i] = 0
            # viscosity, Pa S
            vis[i] = 1e-10
            # heat capacity, J/(mol K)
            cp[i] = 0
            # density, kg/m3
            rho[i] = 0
        i += 1
    # print(cp)
    # calculate the properties of mixture
    cp_m = np.sum(cp * xi_gas)
    phi, denominator = np.ones((n, n)), np.ones((n, n))  # Wilke coefficient
    vis_m, k_m = 0, 0
    for i in range(n):
        for j in np.arange(n):
            phi[i, j] = (1 + (vis[i] / vis[j]) ** 0.5 * (M[j] / M[i]) ** 0.25) ** 2 / (8 * (1 + M[i] / M[j])) ** 0.5
            denominator[i, j] = xi_gas[j] * phi[i, j]  # if i != j else 0
        vis_m += xi_gas[i] * vis[i] / np.sum(denominator[i])
        k_m += xi_gas[i] * k[i] / np.sum(denominator[i])
    prop = np.array([k_m, vis_m, rho_m, cp_m, M_m])
    prop = np.hstack((prop, cp))
    return pd.Series(prop,
                     index=["k", "vis", 'rho', "cp_m", "M"] + ['cp_' + xi_gas.index[i] for i in range(n)])


class VLE:
    """
    calculate the VLE properties including:
    fugacity, fugacity coe, dew point of mixture, flash calculation
    """

    def __init__(self, T, comp):
        """
        Initialize the VLE object.
        :param T: Temperature (K)
        :param comp: component of mix, name in list
        """
        self.T = T
        self.comp = comp
        self.num = len(self.comp)
        self.index = pd.Series(np.arange(self.num), index=self.comp)
        self.Tc, self.Pc, self.Psat, self.Omega = self._calculate_properties()
        self.a, self.b = self._calculate_parameters()
        self.k_a = self._mix_rule()

    def _calculate_properties(self):
        """
        Calculate critical properties and Psat for each component.
        """
        Tc = np.zeros(self.num)
        Pc = np.zeros(self.num)
        Psat = np.zeros(self.num)
        Omega = np.zeros(self.num)

        for i in range(self.num):
            Tc[i] = PropsSI('Tcrit', self.comp[i])
            Pc[i] = PropsSI('Pcrit', self.comp[i]) * 1e-5
            Psat[i] = PropsSI('Pp', 'T', self.T, 'Q', 1, self.comp[i]) * 1e-5 if self.T < Tc[i] else 1e5
            Omega[i] = PropsSI('acentric', self.comp[i])

        return Tc, Pc, Psat, Omega

    def _calculate_parameters(self):
        """
        Calculate a and b parameters for each component.
        """
        alpha = (1 + (0.48 + 1.574 * self.Omega - 0.176 * self.Omega ** 2)
                 * (1 - (self.T / self.Tc) ** 0.5)) ** 2
        a = 0.42748 * (R * 10) ** 2 * self.Tc ** 2 * alpha / self.Pc
        b = 0.08664 * (R * 10) * self.Tc / self.Pc
        return a, b

    def _mix_rule(self):
        """
        define the mixing rule parameter
        """
        if self.comp == ["CO2", "H2", "Methanol", "H2O", "CO"]:
            # mixing rule parameter CO2 H2 CH3OH H2O CO

            k = np.array([[0, -0.3462, 0.0148, 0.0737, 0],
                          [-0.3462, 0, 0, 0, 0.0804],
                          [0.0148, 0, 0, -0.0789, 0],
                          [0.0737, 0, -0.0789, 0, 0],
                          [0, 0.0804, 0, 0, 0]])
            # k = np.array([[0, -0.0007, 0.1, 0.3, 0.1164],
            #               [0.1164, 0, -0.125, -0.745, -0.0007],
            #               [0.1, -0.125, 0, -0.075, -0.37],
            #               [0.3, -0.745, -0.075, 0, -0.474],
            #               [0.1164, -0.0007, -0.37, -0.474, 0]])
        else:
            k = np.zeros((self.num, self.num))
        return k

    @staticmethod
    def para_mix(comp, a, b, kij):
        """
        Calculate the mixing parameters for EoS.
        :param comp: Molar fraction of fluid (numpy array)
        :param a: a parameters (numpy array)
        :param b: b parameters (numpy array)
        :param kij: Binary interaction coefficients (numpy array)
        :return: a_mix, b_mix (floats)
        """
        num = len(comp)
        b_mix = np.sum(b * comp)
        a_mix = 0
        for m in range(num):
            for n in range(num):
                a_mix += comp[m] * comp[n] * (a[m] * a[n]) ** 0.5 * (1 - kij[m, n])
        return a_mix, b_mix

    @staticmethod
    def func_z(beta, q, status):
        # phase 1 refers to liquid phase, 0 refers to vapor phase
        if status == 1:
            def func(x):
                return beta + x * (x + beta) * (1 + beta - x) / q / beta - x

            return func
        elif status == 0:
            def func(x):
                return 1 + beta - q * beta * (x - beta) / x / (x + beta) - x

            return func

    @staticmethod
    def func_v(z, K):
        def func(x):
            return np.sum(z * (K - 1) / (1 + x * (K - 1)))

        return func

    @staticmethod
    def func_l(z, K):
        def func(x):
            return np.sum(z * (1 - K) / (x + (1 - x) * K))

        return func

    def phi(self, comp, P, phase=0):
        """
        Calculate fugacity coefficients (phi) for the components in a mixture.
        :param comp: Molar fraction of fluid (pandas Series)
        :param P: Total pressure, bar
        :param phase: Phase (0 for vapor, 1 for liquid)
        :return: Fugacity coefficients (numpy array), compression factor of mixture
        """
        # extract the fluid parameters
        index_list = comp.index.tolist()
        index = self.index.loc[index_list]
        num = len(index)
        a, b = self.a[index], self.b[index]
        k_mix = self.k_a[index][:, index]
        [q_ba, a_ba] = np.zeros((2, num))

        # the mix para in EoS for vapor phase
        a_mix, b_mix = self.para_mix(comp, a, b, k_mix)
        beta_mix = b_mix * P / (R * 10) / self.T

        q_mix = a_mix / b_mix / (R * 10) / self.T
        Z_guess = 1e-3 if phase == 1 else 0.8
        Z_mix = opt.fsolve(self.func_z(beta_mix, q_mix, phase), [Z_guess])[0]
        Z_mix = 1e-3 if Z_mix < 0 else Z_mix
        I_mix = np.log((Z_mix + beta_mix) / Z_mix)
        ln_phi = np.empty(num)
        for j in range(num):
            # cycle for each component
            a_ba[j] = 0
            for m in range(num):
                a_ba[j] += (a[j] * a[m]) ** 0.5 * (1 - k_mix[j, m]) * comp[m] * 2
            a_ba[j] -= a_mix
            q_ba[j] = q_mix * (1 + a_ba[j] / a_mix - b[j] / b_mix)
            ln_phi[j] = b[j] * (Z_mix - 1) / b_mix - np.log(Z_mix - beta_mix) - q_ba[j] * I_mix
        phi = np.exp(ln_phi)
        return phi, Z_mix

    def dew_p(self, y, x_guess=None):
        """
        Calculate dew point pressure for a given component or set of components and initial liquid phase composition.

        :param y: Molar fraction of equilibrium vapor (pandas Series)
        :param x_guess: Initial liquid phase composition
        :return: Dew point pressure and composition (dict)
        """
        num = len(y)
        comp = pd.DataFrame(index=['V', 'L1'], columns=y.index)
        comp.iloc[0] = y.values

        # find the equilibrium pressure and mol fraction of liquid phase
        comp.iloc[1] = x_guess if x_guess is not None else self.Psat / np.sum(self.Psat)  # y.values
        P_min = np.min(self.Psat) if len(y) == 2 else 48  # when only CH3OH\H2O is considered, P_min decided by Psat
        step = 0.01 if (len(y) == 2 and self.T < 373) else 1
        for P in np.arange(P_min, 200, step):
            delta_K_sum = 1e5
            K_sum_K_pre = 10
            while delta_K_sum > 0.05:
                phi = np.empty((2, num))
                for i in range(2):
                    # cycle for each phase
                    phi[i], _ = self.phi(comp.iloc[i], P, phase=i)
                K = phi[1] / phi[0]
                if np.isnan(np.sum(K)): return None
                K_sum_cal = np.sum(comp.iloc[0].values / K)
                comp.iloc[1] = comp.iloc[0] / K / K_sum_cal
                delta_K_sum = abs(K_sum_cal - K_sum_K_pre)
                K_sum_K_pre = K_sum_cal

            if abs(K_sum_cal - 1) < 0.01:
                res = {'Pp': P, "K": K, "phi": phi, "comp": comp}
                return res
            elif K_sum_cal > 1:
                res = {'Pp': P, "K": K, "phi": phi, "comp": comp}
                return res

    def flash(self, P, mix):
        """
        Perform flash calculations to find the phase equilibrium composition at a specified pressure.

        :param P: Pressure
        :param mix: Molar fraction of mix fluid (pandas Series)
        :return: Phase equilibrium composition (DataFrame)
        """
        comp = pd.DataFrame(index=['V', 'L1'], columns=mix.index)
        fi = np.zeros((2, self.num))
        fi[1] = 1
        K = np.exp(5.37 * (1 + self.Omega) * (1 - 1 / (self.T / self.Tc))) / (P / self.Pc)
        m, dev = 0, 1
        tol, max_iter, vol_guess = 1e-3, 50, [1e-3, 5e-3, 1e-2, 0.05, 0.1]
        i = 0
        while dev > tol:
            vol = opt.fsolve(self.func_l(mix.values, K), vol_guess[i])[0]
            comp.iloc[1] = mix.values / (vol + (1 - vol) * K)
            comp.iloc[0] = K * comp.iloc[1]
            phi, _ = self.phi(comp.iloc[0], P, 0)
            fi[0] = phi * P * comp.iloc[0]
            phi, _ = self.phi(comp.iloc[1], P, 1)
            fi[1] = phi * P * comp.iloc[1]
            K = fi[1] * K / fi[0]
            dev = np.max(abs(fi[0] - fi[1]) / fi[0])
            if m % max_iter == 0 and m != 0:
                print(abs(fi[0] - fi[1]) / fi[0])
                print("reach max iteration")
                i += 1
                # return comp, phi
            m += 1
        return comp, phi


class Thermo:
    """
    calculate thermo-physical properties including
    formation of enthalpy, formation of Gibbs energy
    at given comps, temperature
    """

    def __init__(self, hf_para=None, gf_para=None):

        # fit para for calculate enthalpy of formation, and Gibbs energy of formation
        self.gf_paras = self.__Gf_para() if gf_para is None else gf_para
        self.hf_paras = self.__Hf_para() if hf_para is None else hf_para
        self.gf_fit, self.hf_fit = self.fit_para

    @staticmethod
    def __Gf_para():
        # # Parameters for calculating gibbs energy of formation
        # # CO2 H2 CH3OH H2O CO
        Ts = np.arange(300, 1100, 100)
        Gfs = np.array([[-394.379, -394.656, -394.914, -395.152, -395.367, -395.558, -395.724, -395.865],
                        [0, 0, 0, 0, 0, 0, 0, 0],
                        [-162.057, -148.509, -134.109, -119.125, -103.737, -88.063, -72.188, -56.170],
                        [-228.5, -223.9, -219.05, -214.008, -208.814, -203.501, -198.091, -192.603],
                        [-137.333, -146.341, -155.412, -164.480, -173.513, -182.494, -191.417, -200.281]])

        Gfs = pd.DataFrame(Gfs, index=["CO2", "H2", "Methanol", "H2O", "CO"], columns=Ts)
        return Gfs

    @staticmethod
    def __Hf_para():
        # # Parameters for calculating enthalpy of formation
        # # CO2 H2 CH3OH H2O CO
        Ts = np.arange(300, 1100, 100)

        Hfs = np.array([[-393.511, -393.586, -393.672, -393.791, -393.946, -394.133, -394.343, -394.568],
                        [0, 0, 0, 0, 0, 0, 0, 0],
                        [-201.068, -204.622, -207.750, -210.387, -212.570, -214.350, -215.782, -216.916],
                        [-241.844, -242.845, -243.822, -244.751, -245.620, -246.424, -247.158, -247.820],
                        [-110.519, -110.121, -110.027, -110.157, -110.453, -110.870, -111.378, -111.952]])
        Hfs = pd.DataFrame(Hfs, index=["CO2", "H2", "Methanol", "H2O", "CO"], columns=Ts)
        return Hfs

    @property
    def fit_para(self):
        """
        fit para for calculating the energy (Gibbs energy or enthalpy) of formation for a species
        """
        gf_fit_paras = {}
        hf_fit_paras = {}
        for comp in self.gf_paras.index:
            gf_fit_paras[comp] = np.polyfit(self.gf_paras.loc[comp].index, self.gf_paras.loc[comp].values, 2)
            hf_fit_paras[comp] = np.polyfit(self.hf_paras.loc[comp].index, self.hf_paras.loc[comp].values, 2)
        return gf_fit_paras, hf_fit_paras

    def sat(self, pi):
        """
        :param pi: partial pressure (pd.Series)
        :return:
        """
        pc = pd.Series(index=pi.index)
        for comp in pi.index:
            pc = PropsSI('Pcrit', comp) * 1e-5

    def H(self, T, comps):
        """
        calculate formation of enthalpy at T for given comp
        :param T: K
        :param comps: list
        :return: molar enthalpy (pd.Series)
        """
        H = pd.Series(index=comps)
        for comp in comps:
            H.loc[comp] = np.polyval(self.hf_fit[comp], T)
        return H

    def G(self, T, pi):
        """
        calculate formation of gibbs energy at T for given partial pressure
        :param T: K
        :param pi: partial pressure (pd.Series)
        :return: molar gibbs energy (pd.Series)
        """
        G = pd.Series(index=pi.index)
        for comp in pi.index:
            if pi[comp] < 1e-15:
                G.loc[comp] = 0
            else:
                G.loc[comp] = np.polyval(self.gf_fit[comp], T) + 8.314 * T * np.log(pi[comp] / 1) / 1000
        return G


class VLEThermo:
    def __init__(self, comp, ref='ch', kij=None):
        self.comp = comp.copy()
        self.num = len(comp)
        for i in range(self.num):
            self.comp[i] = 'carbon monoxide' if self.comp[i] == 'CO' else self.comp[i]
        self.kij = kij
        self.const, self.cor = ChemicalConstantsPackage.from_IDs(self.comp)
        self.cp, self.eos_kw = self.eos_paras()
        # determine the properties at ref state
        self.ref = ref
        self.Sref, self.Href, self.Gref = self.ref_data()

    def ref_data(self):
        if self.ref == 'ch':
            return self.__cfs()
        elif self.ref == 'ev':
            return self.__tds()
        else:
            raise ValueError('the reference state is not available')

    def eos_paras(self, eos='SRK'):
        # n_count = len(self.comp)
        cp_cal = []
        for i in range(self.num):
            try:
                cp_cal.append(HeatCapacityGas(CASRN=self.const.CASs[i], MW=self.const.MWs[i], method='TRCIG'))
            except ValueError:
                cp_cal.append(HeatCapacityGas(CASRN=self.const.CASs[i], MW=self.const.MWs[i], method='POLING_POLY'))

        if eos == 'SRK':
            if self.num == 1:
                kijs = np.array([0])
            else:
                kijs = np.zeros((self.num, self.num))
                kijs[:5, :5] = np.array([[0, -0.3462, 0.0148, 0.0737, 0],
                                         [-0.3462, 0, 0, 0, 0.0804],
                                         [0.0148, 0, 0, -0.0789, 0],
                                         [0.0737, 0, -0.0789, 0, 0],
                                         [0, 0.0804, 0, 0, 0]]) if self.kij is None else self.kij

            eos_kwargs = dict(Tcs=np.array(self.const.Tcs), Pcs=np.array(self.const.Pcs),
                              omegas=np.array(self.const.omegas), kijs=kijs)



        elif eos == 'MSRK':
            p = [0, 0, 0.2359, 0.1277, 0, 0]
            kijs = np.zeros((self.num, self.num))
            kijs[:5, :5] = np.array([[0, 0.1164, 0.1, 0.3, 0.1164],
                                     [0.1164, 0, -0.125, -0.745, -0.0007],
                                     [0.1, -0.125, 0, -0.075, -0.37],
                                     [0.3, -0.745, -0.075, 0, -0.474],
                                     [0.1164, -0.0007, -0.37, -0.474, 0]])
            eos_kwargs = dict(Tcs=np.array(self.const.Tcs), Pcs=np.array(self.const.Pcs),
                              omegas=np.array(self.const.omegas), kijs=kijs, S2s=p)

        return cp_cal, eos_kwargs

    def __local_db(self):
        pass

    def __tds(self):
        """
        the environmental reference state
        10.1016/j.apenergy.2016.10.103
        """

        s = [0.21389569, 0.130699312, 0.126593996, 0.070001677, 0.197809827] + [0] * (self.num - 5)
        h = [83.772, 274.158, 753.974, 20.871, 334.259] + [0] * (self.num - 5)

        s = np.array(s) * 1e3
        h = np.array(h) * 1e3
        g = h - 298.15 * s
        s_dep, h_dep = np.zeros(len(self.comp)), np.zeros(len(self.comp))

        h_dep[2:4] = np.array([-201 - (-239.2), -241.8 - (-285.8)]) * 1000
        s_dep[2:4] = np.array([239.9 - 126.8, 188.8 - 70])
        s += s_dep
        h += h_dep
        return np.array(s), np.array(h), np.array(g) * 1e3

    def __cfs(self):
        """
        the formation energy at std state ideal gas state
        """
        g = np.array(self.const.Hfgs) - 298.15 * np.array(self.const.Sfgs)
        return np.array(self.const.Sfgs), np.array(self.const.Hfgs), g

    def _init_apk(self):
        try:
            # Ensure the Aspen Plus application is accessible
            self.aspen = win32com.client.gencache.EnsureDispatch("Apwn.Document")
            # Open an existing Aspen Plus file
            self.aspen.InitFromArchive2(self.ASPEN_apk)
            # Make Aspen Plus visible
            self.aspen.Visible = False  # Set to True if you want to see the Aspen Plus GUI
            # self.aspen.Quit()
        except Exception as e:
            print(f"An error occurred: {e}")
            self.aspen.Quit()
            exit(1)

    def cal_cp_ig(self, T):
        cps = []
        for i in self.cp:
            cps.append(i.T_dependent_property(T))
        return np.array(cps)

    def cal_H_ig(self, T, F):
        Higs = []
        for i in self.cp:
            Higs.append(i.T_dependent_property_integral(298.15, T) / 1000)
        Higs = np.array(Higs) + self.Href / 1000
        return np.dot(Higs, F)

    def cal_G(self, T, P, x):
        frac = x / np.sum(x)
        gas = CEOSGas(SRKMIX, HeatCapacityGases=self.cp, eos_kwargs=self.eos_kw)
        liquid = CEOSLiquid(SRKMIX, HeatCapacityGases=self.cp, eos_kwargs=self.eos_kw)
        flasher = FlashVL(self.const, self.cor, liquid=liquid, gas=gas)
        PT = flasher.flash(zs=frac, T=T, P=P * 1E5)
        G_dep = PT.G()
        H_ref = np.dot(self.Href, frac)
        S_ref = np.dot(self.Sref, frac)
        return (G_dep + (H_ref - T * S_ref)) * np.sum(x) / 1000

    def cal_S(self, T, P, x):
        """
        calculate the entropy flux of stream
        :param T:
        :param P:
        :param x: molar flux of each component
        :return: kW/K
        """
        frac = x / np.sum(x)
        gas = CEOSGas(SRKMIX, HeatCapacityGases=self.cp, eos_kwargs=self.eos_kw)
        liquid = CEOSLiquid(SRKMIX, HeatCapacityGases=self.cp, eos_kwargs=self.eos_kw)
        flasher = FlashVL(self.const, self.cor, liquid=liquid, gas=gas)
        PT = flasher.flash(zs=frac, T=T, P=P * 1E5)
        S_dep = PT.S()
        S_ref = np.dot(self.Sref, frac)
        return (S_dep + S_ref) * np.sum(x) / 1000

    def cal_H(self, T, P, x):
        """
        calculate the enthalpy flux of stream
        :param T:
        :param P: bar
        :param x: molar flux of each component
        :return: kW
        """
        "The reference state for most subclasses is an ideal-gas enthalpy of zero at 298.15 K and 101325 Pa"
        frac = x / np.sum(x)
        gas = CEOSGas(SRKMIX, HeatCapacityGases=self.cp, eos_kwargs=self.eos_kw)
        liquid = CEOSLiquid(SRKMIX, HeatCapacityGases=self.cp, eos_kwargs=self.eos_kw)
        flasher = FlashVL(self.const, self.cor, liquid=liquid, gas=gas)
        PT = flasher.flash(zs=frac, T=T, P=P * 1E5)
        H_dep = PT.H()
        H_ref = np.dot(self.Href, frac)
        return (H_dep + H_ref) * np.sum(x) / 1000

    def cal_E(self, T, P, x):
        """
        :return: the exergy of material, kW
        """
        return self.cal_H(T, P, x) - self.cal_S(T, P, x) * 298.15

    def cal_alpha(self, T, P, x):
        """
        :return: the energy quality of material
        """
        return self.cal_E(T, P, x) / self.cal_H(T, P, x)

    def p_dew(self, T, x):
        frac = x / np.sum(x)
        gas = CEOSGas(SRKMIX, HeatCapacityGases=self.cp, eos_kwargs=self.eos_kw)
        liquid = CEOSLiquid(SRKMIX, HeatCapacityGases=self.cp, eos_kwargs=self.eos_kw)
        flasher = FlashVL(self.const, self.cor, liquid=liquid, gas=gas)
        try:
            return flasher.flash(zs=frac, T=T, VF=1).Pp / 1E5
        except UnboundLocalError:
            return 15e5

    def T_dew(self, P, x):
        frac = x / np.sum(x)
        gas = CEOSGas(SRKMIX, HeatCapacityGases=self.cp, eos_kwargs=self.eos_kw)
        liquid = CEOSLiquid(SRKMIX, HeatCapacityGases=self.cp, eos_kwargs=self.eos_kw)
        flasher = FlashVL(self.const, self.cor, liquid=liquid, gas=gas)
        try:
            return flasher.flash(zs=frac, P=P * 1e5, VF=1).T
        except UnboundLocalError:
            return 300

    def T_pub(self, P, x):
        frac = x / np.sum(x)
        gas = CEOSGas(SRKMIX, HeatCapacityGases=self.cp, eos_kwargs=self.eos_kw)
        liquid = CEOSLiquid(SRKMIX, HeatCapacityGases=self.cp, eos_kwargs=self.eos_kw)
        flasher = FlashVL(self.const, self.cor, liquid=liquid, gas=gas)
        try:
            return flasher.flash(zs=frac, P=P * 1e5, VF=0).T
        except UnboundLocalError:
            return 300

    def p_bub(self, T, x):
        frac = x / np.sum(x)
        gas = CEOSGas(SRKMIX, HeatCapacityGases=self.cp, eos_kwargs=self.eos_kw)
        liquid = CEOSLiquid(SRKMIX, HeatCapacityGases=self.cp, eos_kwargs=self.eos_kw)
        flasher = FlashVL(self.const, self.cor, liquid=liquid, gas=gas)
        try:
            try:
                return flasher.flash(zs=frac, T=T, VF=0).Pp / 1E5
            except AttributeError:
                return flasher.flash(zs=frac, T=T, VF=0).P / 1e5
        except UnboundLocalError:
            # print('no bubble point')
            return 15e5

    def phi(self, T, P, x):
        frac = x / np.sum(x)
        gas = CEOSGas(SRKMIX, HeatCapacityGases=self.cp, eos_kwargs=self.eos_kw, T=T, P=P * 1e5, zs=frac)
        liquid = CEOSLiquid(SRKMIX, HeatCapacityGases=self.cp, eos_kwargs=self.eos_kw, T=T, P=P * 1e5, zs=frac)
        flasher = FlashVL(self.const, self.cor, liquid=liquid, gas=gas)
        TP = flasher.flash(zs=frac, T=T, P=P * 1E5)
        phi_res = TP.phis()
        return phi_res

    def z(self, T, P, x):
        frac = x / np.sum(x)
        gas = CEOSGas(SRKMIX, HeatCapacityGases=self.cp, eos_kwargs=self.eos_kw, T=T, P=P * 1e5, zs=frac)
        liquid = CEOSLiquid(SRKMIX, HeatCapacityGases=self.cp, eos_kwargs=self.eos_kw)
        flasher = FlashVL(self.const, self.cor, liquid=liquid, gas=gas)
        TP = flasher.flash(zs=frac, T=T, P=P * 1E5)
        # print(frac)
        return TP.Z()

    def flash(self, T, P, x):
        """
        flash model realize volatile component separation
        :param T: flash temperature
        :param P: flash pressure
        :param x: input molar flux for each componet
        :return: molar fraction of output gas, molar fraction of output liquid, separation ratio of gas
        """
        frac = x / np.sum(x)
        gas = CEOSGas(SRKMIX, HeatCapacityGases=self.cp, eos_kwargs=self.eos_kw, T=T, P=P * 1E5, zs=frac)
        liquid = CEOSLiquid(SRKMIX, HeatCapacityGases=self.cp, eos_kwargs=self.eos_kw)
        flasher = FlashVLN(self.const, self.cor, liquids=[liquid, liquid], gas=gas)
        TP = flasher.flash(zs=frac, T=T, P=P * 1E5)
        P_b = self.p_bub(T, x)
        sf = TP.VF
        try:
            flasher_gas = TP.gas.zs
        except AttributeError as e:
            # print(e)
            flasher_gas = None if P > P_b else TP.liquid1.zs
            sf = 0 if P > P_b else TP.liquids_betas[1]
            if flasher_gas is None:
                print('no gas in flash')
                exit(0)
        try:
            flasher_liq = TP.liquid0.zs
        except AttributeError as e:
            # print(e)
            flasher_liq = None
        return flasher_gas, flasher_liq, sf

    def cal_Hlg(self, T, x):
        x = np.array(x)
        frac = x / np.sum(x)
        gas = CEOSGas(SRKMIX, HeatCapacityGases=self.cp, eos_kwargs=self.eos_kw)
        liquid = CEOSLiquid(SRKMIX, HeatCapacityGases=self.cp, eos_kwargs=self.eos_kw)
        flasher = FlashVL(self.const, self.cor, liquid=liquid, gas=gas)
        TV = flasher.flash(zs=frac, T=T, VF=1)
        TL = flasher.flash(zs=frac, T=T, VF=0)
        dH = TV.H() - TL.H()
        # print(frac)
        return dH

    def flash_q(self, T, P, x):
        frac = x / np.sum(x)
        gas = CEOSGas(SRKMIX, HeatCapacityGases=self.cp, eos_kwargs=self.eos_kw, T=T, P=P * 1E5, zs=frac)
        liquid = CEOSLiquid(SRKMIX, HeatCapacityGases=self.cp, eos_kwargs=self.eos_kw)
        flasher = FlashVL(self.const, self.cor, liquid=liquid, gas=gas)
        TP = flasher.flash(zs=frac, T=T, P=P * 1E5)
        dH = gas.H() - TP.H()
        # print(frac)
        return dH

    def cal_cp(self, T, P, x):
        frac = x / np.sum(x)
        gas = CEOSGas(SRKMIX, HeatCapacityGases=self.cp, eos_kwargs=self.eos_kw, T=T, P=P * 1E5, zs=frac)
        return gas.Cp()  # J/mol K

    def cal_T_from_H(self, H, P, x):
        """
        calculate temperature for constant pressure and enthalpy
        :param H: enthalpy flux of fluid, kw
        :param P: fluid pressure, bar
        :param x: molar flux of component, mol/s
        :return:
        """
        frac = x / np.sum(x)
        H_ref = np.dot(self.Href, frac)
        gas = CEOSGas(SRKMIX, HeatCapacityGases=self.cp, eos_kwargs=self.eos_kw)
        liquid = CEOSLiquid(SRKMIX, HeatCapacityGases=self.cp, eos_kwargs=self.eos_kw)
        flasher = FlashVL(self.const, self.cor, liquid=liquid, gas=gas)
        PT = flasher.flash(zs=frac, H=(H * 1000) / np.sum(x) - H_ref, P=P * 1E5)
        return np.round(PT.T, 3)

    def cal_T_from_S(self, S, P, x):
        """
        calculate temperature for constant pressure and enthalpy
        :param S: entropy flux of fluid, kw
        :param P: fluid pressure, bar
        :param x: molar flux of component, mol/s
        :return:
        """
        frac = x / np.sum(x)
        S_ref = np.dot(self.Sref, frac)
        gas = CEOSGas(SRKMIX, HeatCapacityGases=self.cp, eos_kwargs=self.eos_kw)
        liquid = CEOSLiquid(SRKMIX, HeatCapacityGases=self.cp, eos_kwargs=self.eos_kw)
        flasher = FlashVL(self.const, self.cor, liquid=liquid, gas=gas)
        PT = flasher.flash(zs=frac, S=(S * 1000) / np.sum(x) - S_ref, P=P * 1E5)
        return np.round(PT.T, 3)

    def cal_rho(self, T, P, x):
        """
        calculate the mass density of mixture, kg/m3
        """
        frac = x / np.sum(x)
        gas = CEOSGas(SRKMIX, HeatCapacityGases=self.cp, eos_kwargs=self.eos_kw, zs=frac)
        liquid = CEOSLiquid(SRKMIX, HeatCapacityGases=self.cp, eos_kwargs=self.eos_kw)
        flasher = FlashVL(self.const, self.cor, liquid=liquid, gas=gas)
        PT = flasher.flash(zs=frac, T=T, P=P * 1E5)
        return PT.MW() / PT.V() / 1000

    def cal_rho_mol(self, T, P, x):
        """
        calculate the molar density of mixture, mol/m3
        """
        frac = x / np.sum(x)
        gas = CEOSGas(SRKMIX, HeatCapacityGases=self.cp, eos_kwargs=self.eos_kw, zs=frac)
        liquid = CEOSLiquid(SRKMIX, HeatCapacityGases=self.cp, eos_kwargs=self.eos_kw)
        flasher = FlashVL(self.const, self.cor, liquid=liquid, gas=gas)
        PT = flasher.flash(zs=frac, T=T, P=P * 1E5)
        return 1 / PT.V()

    def vis(self, T, P, x):
        """
        calculate the molar density of mixture, mol/m3
        """
        frac = x / np.sum(x)
        gas = CEOSGas(SRKMIX, HeatCapacityGases=self.cp, eos_kwargs=self.eos_kw, zs=frac)
        liquid = CEOSLiquid(SRKMIX, HeatCapacityGases=self.cp, eos_kwargs=self.eos_kw)
        flasher = FlashVL(self.const, self.cor, liquid=liquid, gas=gas)
        PT = flasher.flash(zs=frac, T=T, P=P * 1E5)
        print(PT.MW())
        return PT.mu()

    def mw(self, T, P, x):
        """
        calculate the molar density of mixture, mol/m3
        """
        frac = x / np.sum(x)
        gas = CEOSGas(SRKMIX, HeatCapacityGases=self.cp, eos_kwargs=self.eos_kw, zs=frac)
        liquid = CEOSLiquid(SRKMIX, HeatCapacityGases=self.cp, eos_kwargs=self.eos_kw)
        flasher = FlashVL(self.const, self.cor, liquid=liquid, gas=gas)
        PT = flasher.flash(zs=frac, T=T, P=P * 1E5)
        return PT.MW()


class ThermoASPEN:
    def __init__(self, apk_name):

        self._init_apk(apk_name)

        self.apk_comp = ['H2', 'CO', 'CO2', 'H2O', 'MEOH']
        self.Sref, self.Href, self.Gref = self.__tds()
        self.Sref_ig, self.Href_ig, self.Gref_ig = self.__igs()
        self.Ms = pd.Series([2.016, 28.01, 44.01, 18.015, 32.042], index=self.apk_comp)

    def __igs(self):
        h = [3.7253E-15, -110.53, -393.51, -241.818, -200.94]
        s = [0, 0.0892839, 0.00288445, -0.0444273, -0.1295321]
        s = np.array(s)
        h = np.array(h)
        g = h - 298.15 * s
        return pd.Series(s, index=self.apk_comp), pd.Series(h, index=self.apk_comp), pd.Series(g, index=self.apk_comp)

    def __tds(self):
        """
        the environmental reference state
        10.1016/j.apenergy.2016.10.103
        """
        # the thermodynamic value at 298.15 K for ideal state
        s = [0.130700773, 0.197832827, 0.21398957, 0.196497077, 0.250353096]  # kJ/mol/K
        h = [274.1571269, 334.2663, 83.8124, 67.9146, 795.7183]  # kJ/mol

        s = np.array(s)
        h = np.array(h)
        g = h - 298.15 * s
        return pd.Series(s, index=self.apk_comp), pd.Series(h, index=self.apk_comp), pd.Series(g, index=self.apk_comp)

    @staticmethod
    def convert_index(stream_source):
        stream = stream_source.copy()
        new_index = []
        for i in stream.index.tolist():
            if i == 'carbon monoxide':
                sub = 'CO'
            elif i == 'Methanol':
                sub = 'MEOH'
            else:
                sub = i
            new_index.append(sub)
        stream.index = new_index
        return stream

    def _init_apk(self, apk):
        try:
            # Ensure the Aspen Plus application is accessible
            self.aspen = win32com.client.gencache.EnsureDispatch("Apwn.Document")
            # Open an existing Aspen Plus file
            self.aspen.InitFromArchive2(apk)
            # Make Aspen Plus visible
            self.aspen.Visible = False  # Set to True if you want to see the Aspen Plus GUI
            # self.aspen.Quit()
        except Exception as e:
            print(f"An error occurred: {e}")
            self.aspen.Quit()
            exit(1)

    def set_M(self, stream_source, prop):
        """
        Set the feed parameters for the specified stream in Aspen Plus.
        :param stream_source: Stream info, pd.Seize.
        :param prop: property to be calculated, DHM or DSM
        """
        prop_name = 'HMX' if prop == 'H' else 'SMX'
        try:
            # Access the feed stream node
            stream = self.convert_index(stream_source)
            stream[self.apk_comp] = stream[self.apk_comp] / stream[self.apk_comp].sum()
            temp_node = self.aspen.Tree.FindNode(f"\\Data\\Properties\\Analysis\\{prop_name}\\Input\\LIST\\#1\\#0")
            pres_node = self.aspen.Tree.FindNode(f"\\Data\\Properties\\Analysis\\{prop_name}\\Input\\LIST\\#0\\#0")
            temp_node.Value = stream['T']
            pres_node.Value = stream['P']
            for comp in self.apk_comp:
                # Application.Tree.FindNode("\Data\Properties\Analysis\DHM\Input\FLOW\CO")
                comp_node = self.aspen.Tree.FindNode(f"\\Data\\Properties\\Analysis\\{prop_name}\\Input\\FLOW\\{comp}")
                try:
                    comp_node.Value = stream[comp]
                except KeyError:
                    comp_node.Value = 0

        except Exception as e:
            print(f"Error setting feed parameters: {e}")
            self.aspen.Quit()
            exit(1)

    def set_DM(self, stream_source, prop):
        """
        Set the feed parameters for the specified stream in Aspen Plus.
        :param stream_source: Stream info, pd.Seize.
        :param prop: property to be calculated, DHM or DSM
        """
        prop_name = 'DHM' if prop == 'H' else 'DSM'
        try:
            # Access the feed stream node
            stream = self.convert_index(stream_source)
            stream[self.apk_comp] = stream[self.apk_comp] / stream[self.apk_comp].sum()
            temp_node = self.aspen.Tree.FindNode(f"\\Data\\Properties\\Analysis\\{prop_name}\\Input\\LIST\\#1\\#0")
            pres_node = self.aspen.Tree.FindNode(f"\\Data\\Properties\\Analysis\\{prop_name}\\Input\\LIST\\#0\\#0")
            temp_node.Value = stream['T']
            pres_node.Value = stream['P']
            for comp in self.apk_comp:
                # Application.Tree.FindNode("\Data\Properties\Analysis\DHM\Input\FLOW\CO")
                comp_node = self.aspen.Tree.FindNode(f"\\Data\\Properties\\Analysis\\{prop_name}\\Input\\FLOW\\{comp}")
                try:
                    comp_node.Value = stream[comp]
                except KeyError:
                    comp_node.Value = 0

        except Exception as e:
            print(f"Error setting feed parameters: {e}")
            self.aspen.Quit()
            exit(1)

    def set_MIG(self, stream_source, prop):
        """
        Set the feed parameters for the specified stream in Aspen Plus.
        :param stream_source: Stream info, pd.Seize.
        :param prop: property to be calculated, DHM or DSM
        """
        prop_name = 'HMIG' if prop == 'H' else 'SMIG'
        try:
            # Access the feed stream node
            stream = self.convert_index(stream_source)
            stream[self.apk_comp] = stream[self.apk_comp] / stream[self.apk_comp].sum()
            temp_node = self.aspen.Tree.FindNode(f"\\Data\\Properties\\Analysis\\{prop_name}\\Input\\LIST\\#1\\#0")
            temp_node.Value = stream['T']
            for comp in self.apk_comp:

                comp_node = self.aspen.Tree.FindNode(f"\\Data\\Properties\\Analysis\\{prop_name}\\Input\\FLOW\\{comp}")
                try:
                    comp_node.Value = stream[comp]
                except KeyError:
                    comp_node.Value = 0

        except Exception as e:
            print(f"Error setting feed parameters: {e}")
            self.aspen.Quit()
            exit(1)

    def set_IG(self, stream_source, prop):
        """
        Set the feed parameters for the specified stream in Aspen Plus.
        :param stream_source: Stream info, pd.Seize.
        :param prop: property to be calculated, DHM or DSM
        """
        prop_name = 'HIG' if prop == 'H' else 'SIG'
        try:
            # Access the feed stream node
            stream = stream_source.copy()
            temp_node = self.aspen.Tree.FindNode(f"\\Data\\Properties\\Analysis\\{prop_name}\\Input\\TLIST\\#0")
            pres_node = self.aspen.Tree.FindNode(f"\\Data\\Properties\\Analysis\\{prop_name}\\Input\\PLIST\\#0")
            temp_node.Value = stream['T']
            pres_node.Value = stream['P']

        except Exception as e:
            print(f"Error setting feed parameters: {e}")
            self.aspen.Quit()
            exit(1)

    def access_to_DM(self, prop_name):
        prop_key = 'DHM' if prop_name == 'H' else 'DSM'
        prop_value = 'DHMX' if prop_name == 'H' else 'DSMX'
        try:
            node_path = f"\\Data\\Properties\\Analysis\\{prop_key}\\Output\\PROP_DATA\\{prop_value}\\TOTAL\\1"
            res_node = self.aspen.Tree.FindNode(node_path)
            res = res_node.Value
            return res
        except Exception as e:
            print(f"Error accessing stream parameters: {e}")
            self.aspen.Quit()
            exit(1)

    def access_to_M(self, prop_name):
        prop_key = 'HMX' if prop_name == 'H' else 'SMX'
        prop_value = 'HMX' if prop_name == 'H' else 'SMX'
        try:
            # Application.Tree.FindNode("\Data\Properties\Analysis\SMX\Output\PROP_DATA\SMX\TOTAL\1")
            # \Data\Properties\Analysis\SMX\Output\PROP_DATA\SMX\TOTAL\1
            node_path = f"\\Data\\Properties\\Analysis\\{prop_key}\\Output\\PROP_DATA\\{prop_value}\\TOTAL\\1"

            res_node = self.aspen.Tree.FindNode(node_path)
            res = res_node.Value
            return res
        except Exception as e:
            print(f"Error accessing stream parameters: {e}")
            self.aspen.Quit()
            exit(1)

    def access_to_MIG(self, prop_name):
        prop_key = 'HMIG' if prop_name == 'H' else 'SMIG'
        prop_value = 'HIGMX' if prop_name == 'H' else 'SIGMX'
        try:
            node_path = f"\\Data\\Properties\\Analysis\\{prop_key}\\Output\\PROP_DATA\\{prop_value}\\VAPOR\\1"
            res_node = self.aspen.Tree.FindNode(node_path)
            res = res_node.Value
            return res
        except Exception as e:
            print(f"Error accessing stream parameters: {e}")
            self.aspen.Quit()
            exit(1)

    def access_to_IG(self, stream_source, prop_name):
        stream = self.convert_index(stream_source)
        prop_key = 'HIG' if prop_name == 'H' else 'SIG'

        try:
            node_path = f"\\Data\\Properties\\Analysis\\{prop_key}\\Output\\Prop\nData\\PROPTAB"
            res = pd.Series(index=self.apk_comp)
            for i in self.apk_comp:
                res_node = node_path + f"\\VAPOR\nHIG\n{i}\n1"
                res[i] = self.aspen.Tree.FindNode(res_node).Value
            res_value = (res * stream[:-2]).sum()
            return res_value
        except Exception as e:
            print(f"Error accessing stream parameters: {e}")
            self.aspen.Quit()
            exit(1)

    def cal_H(self, stream):
        f_total = stream[:-2].sum()
        stream_flow = self.convert_index(stream[:-2])
        self.set_DM(stream, prop='H')
        self.set_MIG(stream, prop='H')
        self.aspen.Engine.Run2()
        H_dep = self.access_to_DM(prop_name='H') * f_total
        HIG = self.access_to_MIG('H') * f_total
        HIG_ref = (self.Href_ig * stream_flow).sum()
        HM = H_dep + HIG - HIG_ref + (stream_flow * self.Href).sum()
        return HM

    def cal_S(self, stream):
        f_total = stream[:-2].sum()
        stream_flow = self.convert_index(stream[:-2])
        self.set_DM(stream, prop='S')
        self.set_MIG(stream, prop='S')
        self.aspen.Engine.Run2()
        S_dep = self.access_to_DM(prop_name='S') * f_total
        SIG = self.access_to_MIG('S') * f_total
        SIG_ref = (self.Sref_ig * stream_flow).sum()
        SM = S_dep + SIG - SIG_ref + (stream_flow * self.Sref).sum()
        return SM

    def cal_E(self, stream):
        f_total = stream[:-2].sum()
        stream_flow = self.convert_index(stream[:-2])
        self.set_DM(stream, prop='S')
        self.set_MIG(stream, prop='S')
        self.set_DM(stream, prop='H')
        self.set_MIG(stream, prop='H')
        self.aspen.Engine.Run2()
        S_dep = self.access_to_DM(prop_name='S') * f_total
        SIG = self.access_to_MIG('S') * f_total
        SIG_ref = (self.Sref_ig * stream_flow).sum()
        SM = S_dep + SIG - SIG_ref + (stream_flow * self.Sref).sum()

        H_dep = self.access_to_DM(prop_name='H') * f_total
        HIG = self.access_to_MIG('H') * f_total
        HIG_ref = (self.Href_ig * stream_flow).sum()
        HM = H_dep + HIG - HIG_ref + (stream_flow * self.Href).sum()
        return HM - SM * 298.15

    def cal_prop2(self, stream, quit_aspen=True):
        # print(stream)
        stream_flow = self.convert_index(stream)
        stream_flow = stream_flow[self.apk_comp]
        f_total = stream_flow[self.apk_comp].sum()

        self.set_M(stream, prop='S')
        self.set_MIG(stream, prop='S')
        # self.aspen.Engine.Run2()

        self.set_M(stream, prop='H')
        self.set_MIG(stream, prop='H')
        self.aspen.Engine.Run2()

        S = self.access_to_M(prop_name='S') * f_total
        SIG_ref = (self.Sref_ig * stream_flow).sum()
        SM = S - SIG_ref + (stream_flow * self.Sref).sum()

        H = self.access_to_M(prop_name='H') * f_total
        HIG_ref = (self.Href_ig * stream_flow).sum()
        HM = H - HIG_ref + (stream_flow * self.Href).sum()
        EM = HM - SM * 298.15
        GM = HM - SM * stream['T']

        # print(HM,SM,EM)
        if quit_aspen:
            self.aspen.Quit()
        return pd.Series([HM, SM, EM, GM], index=['H', 'S', 'E', 'G'])

    def eva_streams(self, streams):
        """
        evaluate the properties of multiply stream
        :return:
        """
        for i in streams.index.tolist():
            properties = self.cal_prop2(streams.loc[i], quit_aspen=False)
            streams.loc[i, ['H', 'S', 'E', 'G']] = properties
        self.aspen.Quit()
        return streams

    def convert_to_kg(self, stream_source):

        stream = self.convert_index(stream_source[:-2])
        f_total = stream[self.apk_comp].sum()
        M_ave = (stream[self.apk_comp] * self.Ms[self.apk_comp]).sum() / f_total  # g/mol
        res = self.cal_prop2(stream_source)
        return res / (M_ave / 1000)  # kj/kg

    def perform(self):
        self.aspen.Engine.Run2()

    # def run_RF(self, stream, conds, quit_aspen=True):
    #     stream_num = len(conds)
    #     prop_res = pd.DataFrame(columns=['T', 'P', 'H', 'S', 'G'], index=range(stream_num + 1))
    #     sts_info = self.generate_sts(stream, conds=conds)
    #     for i in range(stream_num + 1):
    #         self.set_feed_parameters(stream_source=sts_info[i])
    #         self.aspen.Engine.Run2()
    #         prop_res.iloc[i] = self.access_to_st()
    #
    #     if quit_aspen:
    #         self.aspen.Quit()
    #     return prop_res
