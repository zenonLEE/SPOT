import numpy as np
from CoolProp.CoolProp import PropsSI
import pandas as pd

from prop_calculator import mixture_property, VLEThermo

R = 8.314  # J/mol/K
ks, vof = 0.05, 0.8  # 0.2, 0.8  # 0.2 1.5 for 0.42 1 for 0.3 0.2 for 0.15 # 1, 0.4 for CO exp 0.35 0.47


class Reaction:
    """
    basic simulation of CO2 to CH3OH
    energy and mass balance are calculated
    """

    def __init__(self, L, D, Dc, n, phi, rho, chem_para, T0, P0, F0, comp, eos, qmh=0):

        # 0 for ideal 1 for SRK
        self.eos = eos

        #
        self.qmh = qmh

        # reactor parameters
        self.L, self.Dt, self.n = L, D, n
        self.Dc = Dc
        self.phi, self.rho = phi, rho
        self.ds = 6e-3  # catalyst effective particle diameter, cylinder

        # prescribed chem data of reaction
        self.comp_list = comp
        self.chem_data = chem_para
        self.react_num = len(self.chem_data["kr"])
        self.comps_num = len(self.comp_list)
        self.react_sto = np.empty((self.react_num, self.comps_num))
        self.kn_model = self.chem_data['kn_model']
        for i in range(self.react_num):
            key = str(i + 1)
            if self.comps_num == 5:
                self.react_sto[i] = self.chem_data["stoichiometry"][key]
            else:
                self.react_sto[i] = self.chem_data["stoichiometry"][key] + [0]

        # feed gas parameter
        self.vle_cal = VLEThermo(self.comp_list)
        self.R = 8.314
        self.P0, self.T0 = P0, T0  # P0 bar, T0 K
        self.F0, self.Ft0 = F0, np.sum(F0)
        z0 = self.vle_cal.z(T=T0, P=P0, x=F0)
        self.v0 = z0 * self.Ft0 * self.R * self.T0 / (self.P0 * 1e5)

    @staticmethod
    def react_H(T, in_dict):
        dH = np.zeros(len(in_dict["heat_reaction"].keys()))
        i = 0
        for key, value in in_dict["heat_reaction"].items():
            dH[i] = -(value[0] * T + value[1]) * 1e-6
            i += 1
        # dH[0] = -49.93
        # dH[1] = 41.12
        return dH

    @staticmethod
    def kad(T, in_dict):
        """
        calculate the equilibrium constant of adsorption
        :param T: operating temperature
        :param in_dict: prescribed chemical parameter
        :return: equilibrium constant of adsorption, 1/bar
        """
        adsorption_eq_constant = dict()
        for key, value in in_dict["kad"].items():
            adsorption_eq_constant[key] = value[0] * np.exp(value[1] / T / R)
        return adsorption_eq_constant

    @staticmethod
    def keq(T, in_dict):
        """
        calculate the equilibrium constant
        :param T: operating temperature
        :param in_dict: prescribed chemical parameter
        :return: equilibrium constant
        """
        react_eq_constant = dict()
        for key, value in in_dict["keq"].items():
            react_eq_constant[key] = 10 ** (value[0] / T + value[1])
        return react_eq_constant

    @staticmethod
    def kr(T, in_dict):
        """
        calculate the reaction rate constant
        :param T: operating temperature, K
        :param in_dict: prescribed chemical parameter
        :return: the reaction rate constant, mol kg−1 s−1 bar-1/2
        """
        react_rate_constant = dict()
        for key, value in in_dict["kr"].items():
            react_rate_constant[key] = value[0] * np.exp(value[1] / T / R)
        return react_rate_constant

    def ergun(self, T, P, F_dict):
        """
        energy and material balance in the reactor
        :param T: operating temperature, K
        :param P: operating pressure, bar
        :param F_dict: molar flow rate of each component, mol/s; ndarray
        :return: pressure drop per length, pa/m
        """
        Ft = np.sum(F_dict)
        xi = F_dict / Ft * P
        z = self.vle_cal.z(T=T, P=P, x=F_dict)
        gas_prop = mixture_property(T, pd.Series(xi, index=self.comp_list), P, z, rho_only=False)
        v = Ft / (gas_prop['rho'] / gas_prop['M'])  # m3/s self.v0 * (self.P0 / Pp) * (T / self.T0) * (Ft / self.Ft0)
        Ae = self.Dt ** 2 / 4 if self.qmh == 0 else (self.Dt ** 2 - self.Dc ** 2) / 4
        u = v / (np.pi * Ae)
        if self.eos == 1:
            # properties = VLE(T, comp=self.comp_list)
            # _, z = properties.phi(pd.Series(F_dict / Ft, index=self.comp_list), Pp)
            z = self.vle_cal.z(T=T, P=P, x=F_dict)
        elif self.eos == 0:
            z = 1
        rho = self.vle_cal.cal_rho(T=T, P=P, x=F_dict)
        gas_property = mixture_property(T, pd.Series(F_dict / Ft, index=self.comp_list), P, z, rho_only=False)
        Re = self.ds * u * rho / gas_property['vis'] / self.phi
        drop_per_length = - (150 / Re + 1.75) * self.phi / (1 - self.phi) ** 3 * (
                gas_property['rho'] * u ** 2 / self.ds)  # Pa/m
        # print(drop_per_length)
        return drop_per_length

    def convection(self, T, P, F_dict):
        """
        :param T: temperature of reactor gas, K
        :param P: pressure of reactor, bar
        :param F_dict: molar flow rate of each component, mol/s; ndarray
        :return: convection heat transfer coefficient, W/m2 K
        """
        Ft = np.sum(F_dict)
        xi = F_dict / Ft * P
        mix_property = mixture_property(T, pd.Series(xi, self.comp_list), Pt=P)
        M = 0.25 * 44 + 0.75 * 2
        Pr = mix_property['vis'] * (mix_property['cp_m'] / (M / 1000)) / mix_property['k']
        v = self.v0 * (self.P0 / P) * (T / self.T0) * (Ft / self.Ft0)  # m3/s
        u = v / (np.pi * self.Dt ** 2 / 4)
        Re = u * self.Dt * mix_property['rho'] / mix_property['vis']
        if Re > 1e4:
            Nu = 0.0265 * Re ** 0.8 * Pr ** 0.3
        elif 2300 < Re < 1e4:
            f = (0.79 * np.log(Re) - 1.64) ** -2
            Nu = f / 8 * (Re - 1000) * Pr / (1 + 12.7 * (f / 8) ** 0.5 * (Pr ** (2 / 3) - 1))
        elif Re < 2300:
            Nu = 3.66
        h = Nu * mix_property['k'] / self.Dt  # W/m K
        return h

    def htr2(self, T, P, F_dict):
        """
        calculate the internal heat transfer coe in within the catalyst layer

        :param T: temperature of reactor gas, K
        :param P: pressure of reactor, bar
        :param F_dict: molar flow rate of each component, mol/s; ndarray
        :return: convection heat transfer coefficient, W/m2 K
        """
        Ft = np.sum(F_dict)
        xi = F_dict / Ft * P
        z = self.vle_cal.z(T=T, P=P, x=F_dict)
        gas_prop = mixture_property(T, pd.Series(xi, index=self.comp_list), P, z, rho_only=False)
        kp = 0.38  # 0.21 + 0.00015 * T  # thermal conductivity of particle

        De = self.Dt if self.qmh == 0 else (self.Dt - self.Dc)  # effective diameter of annual tube
        Ae = self.Dt ** 2 / 4 if self.qmh == 0 else (self.Dt ** 2 - self.Dc ** 2) / 4

        v = Ft / (gas_prop['rho'] / gas_prop['M'])  # m3/s self.v0 * (self.P0 / Pp) * (T / self.T0) * (Ft / self.Ft0)
        u = v / (np.pi * Ae)
        r_kg_ks = gas_prop['k'] / kp

        ke0 = gas_prop['k'] * ((1 - self.phi) / 1.5 +
                               self.phi / (0.13 * (1 - self.phi) ** 1.44 + 2 / 3 * r_kg_ks))
        Pe = 8.65 * (1 + 19.4 * (self.ds / De) ** 2)
        keg = u * gas_prop['rho'] * gas_prop['cp_m'] * self.ds / Pe
        ke = ke0 + keg

        hw0 = (2 * (1 - self.phi) +
               self.phi / (2 / 3 * r_kg_ks + 0.0024 * (De / self.ds) ** 1.58)) * gas_prop['k'] / self.ds
        Re = u * gas_prop['rho'] * self.ds / gas_prop['vis']
        if Re < 1200:
            hwg = 0.0835 * Re ** 0.91 * gas_prop['k'] / self.ds
            if Re < 10:
                print('TOO SMALL Re')
        else:
            hwg = 1.23 * Re ** 0.53 * gas_prop['k'] / self.ds

        hw = hwg + hw0
        Rin = (np.log(self.Dt / self.Dc) / (2 * np.pi * ke) + 1 / (np.pi * self.Dc * hw))  # 1/(W/m K)
        Ro = (np.log(self.Dt / self.Dc) / (2 * np.pi * ke) + 1 / (np.pi * self.Dt * hw))
        return Rin, Ro

    def htr(self, T, P, F_dict):
        """
        calculate the internal heat transfer coe in within the catalyst layer
        Cui, CEJ, 2020
        :param T: temperature of reactor gas, K
        :param P: pressure of reactor, bar
        :param F_dict: molar flow rate of each component, mol/s; ndarray
        :return: convection heat transfer coefficient, W/m2 K
        """
        Ft = np.sum(F_dict)
        xi = F_dict / Ft * P
        z = self.vle_cal.z(T=T, P=P, x=F_dict)
        gas_prop = mixture_property(T, pd.Series(xi, index=self.comp_list), P, z, rho_only=False)
        kp = 0.38  # 0.21 + 0.00015 * T  # thermal conductivity of particle
        De = self.Dt if self.qmh == 0 else (self.Dt - self.Dc)  # effective diameter of annual tube
        Ae = self.Dt ** 2 / 4 if self.qmh == 0 else (self.Dt ** 2 - self.Dc ** 2) / 4

        v = Ft / (gas_prop['rho'] / gas_prop['M'])  # m3/s self.v0 * (self.P0 / Pp) * (T / self.T0) * (Ft / self.Ft0)
        u = v / (np.pi * Ae)
        Re = u * gas_prop['rho'] * self.ds / gas_prop['vis']
        Pr = gas_prop['vis'] * (gas_prop['cp_m'] / gas_prop['M']) / gas_prop['k']
        Pe = Re * Pr
        N = De / self.ds
        PeL = 8 * (2 - (1 - 2 / N) ** 2)
        kappa = kp / gas_prop['k']
        B = 1.25 * 2 * (self.phi / (1 - self.phi)) ** 1.11
        right_term2 = (B * (1 - 1 / kappa) / (1 - B / kappa) ** 2 * np.log(kappa / B) -
                       (B - 1) / (1 - B / kappa) + (B + 1) / 2) * 2 * self.phi ** 0.5 / (1 - B / kappa)
        kr = (1 - self.phi ** 0.5) / (1 / gas_prop['k'] - right_term2)
        ke_r = (kr / gas_prop['k'] + Pe / PeL) * gas_prop['k']
        Nu_w0 = (1.3 + 5 / N) * (kr / gas_prop['k'])
        Nu_ws = 0.3 * Pr ** (1 / 3) * Re ** 0.75
        Nu_m = 0.054 * Pr * Re
        Nu_w = Nu_w0 + 1 / (1 / Nu_ws + 1 / Nu_m)
        hw = Nu_w * gas_prop['k'] / self.ds
        Bi = hw * De / 2 / ke_r
        hk = 1 / (De / 2 / 3 / ke_r * (Bi + 3) / (Bi + 4))
        Ut = 1 / (1 / hw + 1 / hk)
        return 1 / Ut  # 1/(W/m2 K)

    def rate_vi(self, T, Pi):
        """
        calculate the reaction rate
        :param T: operating temperature, K
        :param Pi: partial pressure of each component, bar
        :return: reaction rate of each component for each and all reaction; mol/s/kg_cat
        """

        # convert the partial pressure from ndarray to pd.Series
        Pi = pd.Series(Pi, index=self.comp_list)
        K_H2O = 96808 * np.exp(-51979 / 8.314 / T)
        k_r = 11101.2 * np.exp(-117432 / 8.314 / T)
        Ke = 1 / np.exp(-12.11 + 5319 / T + 1.012 * np.log(T) + 1.144 * 10 ** (-4 * T))
        react_rate = k_r * (Pi["CO2"] - Pi["CO"] * Pi["H2O"] / Ke / Pi["H2"]) / (
                1 + K_H2O * Pi["H2O"] / Pi["H2"]) * 1000

        react_comp_rate = np.zeros((3, 5))
        react_comp_rate[1] = react_rate * self.react_sto[1]
        react_comp_rate[2] = react_rate * self.react_sto[1]
        # print(react_comp_rate)
        return react_comp_rate

    def rate_bu(self, T, Pi):
        """
        calculate the reaction rate
        :param T: operating temperature, K
        :param Pi: partial pressure of each component, bar
        :return: reaction rate of each component for each and all reaction; mol/s/kg_cat
        """
        # convert the partial pressure from ndarray to pd.Series
        Pi = pd.Series(Pi, index=self.comp_list)

        # calculate the reaction constant
        rate_const = self.kr(T, self.chem_data)
        ad_const = self.kad(T, self.chem_data)
        eq_const = self.keq(T, self.chem_data)

        # calculate the rate of each reaction
        react_rate = np.zeros(self.react_num)
        driving = rate_const['1'] * Pi['CO2'] * Pi['H2'] * (
                1 - Pi['H2O'] * Pi["Methanol"] / Pi["H2"] ** 3 / Pi['CO2'] / eq_const['1'])
        inhibiting = (1 + ad_const["H2O/H2"] * Pi['H2O'] / Pi['H2'] +
                      ad_const["H2"] * Pi["H2"] ** 0.5 + ad_const["H2O"] * Pi["H2O"])
        react_rate[0] = driving / inhibiting ** 3

        driving = rate_const['2'] * Pi['CO2'] * (1 - Pi['H2O'] * Pi["CO"] / Pi["H2"] / Pi['CO2'] / eq_const['2'])
        react_rate[1] = driving / inhibiting

        # compute the reaction rate for each component in every reaction
        react_comp_rate = self.react_sto * np.repeat(react_rate, self.comps_num).reshape(self.react_num, self.comps_num)
        react_comp_rate = np.vstack((react_comp_rate, np.sum(react_comp_rate, axis=0).T))

        return react_comp_rate

    def rate_sl(self, T, Pi):
        """
        calculate the reaction rate
        :param T: operating temperature, K
        :param Pi: partial pressure of each component, bar
        :return: reaction rate of each component for each and all reaction; mol/s/kg_cat
        """
        # convert the partial pressure from ndarray to pd.Series
        Pi = pd.Series(Pi, index=self.comp_list)

        # calculate the reaction constant
        rate_const = self.kr(T, self.chem_data)
        # print(rate_const)
        ad_const = self.kad(T, self.chem_data)
        eq_const = self.keq(T, self.chem_data)

        # calculate the rate of each reaction
        react_rate = np.zeros(self.react_num)
        driving = rate_const['1'] * Pi['CO2'] * Pi['H2'] ** 2 * (
                1 - Pi['H2O'] * Pi["Methanol"] / Pi["H2"] ** 3 / Pi['CO2'] / eq_const['1'])
        inhibiting = (ad_const["H2"] * Pi['H2'] ** 0.5 +
                      ad_const["H2O"] * Pi["H2O"] + Pi["Methanol"])
        react_rate[0] = driving / inhibiting ** 2

        driving = rate_const['2'] * Pi['CO2'] * (1 - Pi['H2O'] * Pi["CO"] / Pi["H2"] / Pi['CO2'] / eq_const['2'])
        react_rate[1] = driving / inhibiting

        # compute the reaction rate for each component in every reaction
        react_comp_rate = self.react_sto * np.repeat(react_rate, self.comps_num).reshape(self.react_num, self.comps_num)
        react_comp_rate = np.vstack((react_comp_rate, np.sum(react_comp_rate, axis=0).T))
        # react_comp_rate = np.hstack((react_comp_rate, np.array([0, 0, 0]).reshape(3, 1)))

        return react_comp_rate

    def rate_gr(self, T, Pi):
        """
        calculate the reaction rate
        :param T: operating temperature, K
        :param Pi: partial pressure of each component, bar
        :return: reaction rate of each component for each and all reaction; mol/s/kg_cat
        """

        # convert the partial pressure from ndarray to pd.Series
        Pi = pd.Series(Pi, index=self.comp_list)

        # calculate the reaction constant
        rate_const = self.kr(T, self.chem_data)
        ad_const = self.kad(T, self.chem_data)
        eq_const = self.keq(T, self.chem_data)

        # calculate the rate of each reaction
        react_rate = np.zeros(self.react_num)
        driving = rate_const['1'] * ad_const['CO2'] * (
                Pi['CO2'] * Pi['H2'] ** 1.5 - Pi['H2O'] * Pi["Methanol"] / Pi["H2"] ** 1.5 / eq_const['1'])
        inhibiting = (1 + ad_const["CO"] * Pi['CO'] + ad_const["CO2"] * Pi['CO2']) * \
                     (Pi["H2"] ** 0.5 + ad_const["H2O/H2"] * Pi["H2O"])
        react_rate[0] = driving / inhibiting

        driving = rate_const['2'] * ad_const['CO2'] * (Pi['CO2'] * Pi["H2"] - Pi['H2O'] * Pi["CO"] / eq_const['2'])
        react_rate[1] = driving / inhibiting

        driving = rate_const['3'] * ad_const['CO'] * (
                Pi['CO'] * Pi["H2"] ** 1.5 - Pi['Methanol'] / Pi["H2"] ** 0.5 / eq_const['3'])
        react_rate[2] = driving / inhibiting

        # compute the reaction rate for each component in every reaction
        react_comp_rate = self.react_sto * np.repeat(react_rate, self.comps_num).reshape(self.react_num, self.comps_num)
        react_comp_rate = np.vstack((react_comp_rate, np.sum(react_comp_rate, axis=0).T))

        return react_comp_rate

    def balance(self, T, P, F_dict):
        """
        energy and material balance in the reactor
        :param T: operating temperature, K
        :param P: operating pressure, bar
        :param F_dict: molar flow rate of each component, mol/s; ndarray
        :return: temperature and molar flux variation of gas
        """
        Ft = np.sum(F_dict)  # total molar flow rate
        xi = F_dict / Ft  # molar fraction of mix

        # calculate the partial pressure/fugacity
        # calculate the correction to volumetric flow rate (m3/s)
        if self.eos == 1:
            # fugacity coe, compression factor
            # vle_cal = VLE(T, self.comp_list)
            # phi, _ = vle_cal.phi(comp=pd.Series(xi, index=self.comp_list), Pp=Pp, phase=0)
            # vle_cal = VLEThermo(self.comp_list)
            phi = np.array(self.vle_cal.phi(T=T, P=P, x=xi))
        else:
            phi = 1
        v = self.v0 * (self.P0 / P) * (T / self.T0) * (Ft / self.Ft0)
        Pi = xi * P  # F_dict * R * T / v * 1e-5  # bar
        fi = Pi * phi
        # calculate the change of the molar flow rate due to reactions, mol/s/kg_cat
        if self.kn_model == 'GR':
            dF_react = self.rate_gr(T, fi)
        elif self.kn_model == 'BU':
            dF_react = self.rate_bu(T, fi)  # self.rate_vi(T, fi) #
        elif self.kn_model == 'SL':
            dF_react = self.rate_sl(T, fi)
        elif self.kn_model == "VI":
            dF_react = self.rate_vi(T, fi)

        # calculate the change of enthalpy due to reaction, kJ/(kg_cat s)
        dH_react = self.react_H(T, self.chem_data)
        if self.react_num == 3: dH_react[2] = dH_react[0] - dH_react[1]
        dH = np.matmul(dF_react[:-1, 0], dH_react.T)

        res = {
            'mflux': dF_react[-1],
            'hflux': dH * 1e3
        }
        return res
