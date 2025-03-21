import numpy as np
import pandas as pd

from prop_calculator import VLEThermo, mixture_property


class HMT:
    def __init__(self, sub, phi, ds, Dt, Dm=None, Dc=None):

        self.comp_list = sub
        self.phi = 0.6
        self.vle_cal = VLEThermo(comp=sub)
        self.Dt = Dt
        self.Dm = Dm
        self.Dc = Dc
        self.ds = ds
        self._cal_Ae()
        self._cal_De()

    def _cal_Ae(self):
        if self.Dm is None:
            self.Ae = np.pi * self.Dt ** 2 / 4
        else:
            self.Ae = np.pi * self.Dt ** 2 / 4 - np.pi * self.Dm ** 2 / 4

    def _cal_De(self):
        if self.Dm is None:
            self.De = self.Dt
        else:
            self.De = self.Dt - self.Dm

    def ergun(self, T, P, F):
        """
        energy and material balance in the reactor
        :param T: operating temperature, K
        :param P: operating pressure, bar
        :param F: molar flow rate of each component, mol/s; ndarray
        :return: pressure drop per length, pa/m
        """
        y = F / np.sum(F)
        v = np.sum(F) / self.vle_cal.cal_rho_mol(T, P, y)
        u = v / self.Ae
        z = self.vle_cal.z(T, P, y)
        rho_mass = self.vle_cal.cal_rho(T, P, y)
        gas_property = mixture_property(T, pd.Series(y, index=self.comp_list), P, z, rho_only=False)
        Re = self.ds * u * rho_mass / gas_property['vis'] / self.phi
        drop_per_length = - (150 / Re + 1.75) * self.phi / (1 - self.phi) ** 3 * (
                rho_mass * u ** 2 / self.ds)  # Pa/m

        return drop_per_length

    def htc(self, T, P, F):
        """
        calculate the internal heat transfer coe in within the catalyst layer
        Cui, CEJ, 2020
        :param T: temperature of reactor gas, K
        :param P: pressure of reactor, bar
        :param F: molar flow rate of each component, mol/s; ndarray
        :return: convection heat transfer coefficient, W/m2 K
        """
        Ft = np.sum(F)
        y = F / Ft * P
        z = self.vle_cal.z(T=T, P=P, x=F)
        gas_prop = mixture_property(T, pd.Series(y, index=self.comp_list), P, z, rho_only=False)
        kp = 0.21 + 0.00015 * T  # thermal conductivity of particle 0.38

        y = F / np.sum(F)
        v = np.sum(F) / self.vle_cal.cal_rho_mol(T, P, y)
        u = v / self.Ae
        rho_mass = self.vle_cal.cal_rho(T, P, y)
        cp_mol = self.vle_cal.cal_cp(T, P, y)  # J/mol/K
        MW = self.vle_cal.mw(T, P, y)  # g/mol
        Re = self.ds * u * rho_mass / gas_prop['vis']  # / self.phi

        Pr = gas_prop['vis'] * (cp_mol / (MW / 1000)) / gas_prop['k']
        Pe = Re * Pr
        N = self.De / self.ds
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
        Bi = hw * self.De / 2 / ke_r
        hk = 1 / (self.De / 2 / 3 / ke_r * (Bi + 3) / (Bi + 4))
        Ut = 1 / (1 / hw + 1 / hk)

        print(f'hw: {hw}')
        print(f'hk: {hk}')
        return Ut

    def htc2(self, T, P, F):
        """
        calculate the internal heat transfer coe in within the catalyst layer
        Cui, CEJ, 2020
        :param T: temperature of reactor gas, K
        :param P: pressure of reactor, bar
        :param F: molar flow rate of each component, mol/s; ndarray
        :return: convection heat transfer coefficient, W/m2 K
        """
        Ft = np.sum(F)
        y = F / Ft * P
        z = self.vle_cal.z(T=T, P=P, x=F)
        gas_prop = mixture_property(T, pd.Series(y, index=self.comp_list), P, z, rho_only=False)
        y = F / np.sum(F)
        v = np.sum(F) / self.vle_cal.cal_rho_mol(T, P, y)
        u = v / self.Ae
        rho_mass = self.vle_cal.cal_rho(T, P, y)
        cp_mol = self.vle_cal.cal_cp(T, P, y)  # J/mol/K
        MW = self.vle_cal.mw(T, P, y)  # g/mol
        Re = self.ds * u * rho_mass / gas_prop['vis']
        Pr = gas_prop['vis'] * (cp_mol / (MW / 1000)) / gas_prop['k']

        U = 2.03 * Re ** 0.8 / np.exp(6 * self.ds / self.De) * gas_prop['k'] / self.De
        print(U)
        print(f'Re:{Re}')
        ke_r = (5 + 0.061 * Re) * gas_prop['k']
        print(f'(5 + 0.061 * Re) * k:{ke_r}')
        N = self.ds / self.De
        Nu_w = 2.
        # ke_r = (2.894 + 0.068 * Re) * gas_prop['k']
        # print(f'(2.894 + 0.068 * Re) * k:{ke_r}')
        Nu_w = 0.17 * Re ** 0.79
        print(f'Nu_w: {Nu_w}')
        hw = Nu_w * gas_prop['k'] / self.ds
        print(f'hw:{hw}')

        Ut = 1 / (1 / hw + self.De / 6.13 / ke_r)
        print(f'Ut:{Ut}')
        return Ut

    def htc3(self, T, P, F):
        """
        calculate the internal heat transfer coe in within the catalyst layer
        Cui, CEJ, 2020
        :param T: temperature of reactor gas, K
        :param P: pressure of reactor, bar
        :param F: molar flow rate of each component, mol/s; ndarray
        :return: convection heat transfer coefficient, W/m2 K
        """
        Ft = np.sum(F)
        y = F / Ft * P
        z = self.vle_cal.z(T=T, P=P, x=F)
        gas_prop = mixture_property(T, pd.Series(y, index=self.comp_list), P, z, rho_only=False)
        y = F / np.sum(F)
        v = np.sum(F) / self.vle_cal.cal_rho_mol(T, P, y)
        u = v / self.Ae
        rho_mass = self.vle_cal.cal_rho(T, P, y)
        rho_mol = self.vle_cal.cal_rho_mol(T, P, y)
        cp_mol = self.vle_cal.cal_cp(T, P, y)  # J/mol/K

        kf = gas_prop['k']
        vis = gas_prop['vis']
        N = self.De / self.ds
        Re = self.ds * u * rho_mass / vis
        print(Re)
        phi_ke0 = 0.13 * (1 - self.phi) ** 1.44
        kp = 0.21 + 0.00015 * T  # thermal conductivity of particle 0.38
        ke0 = ((1 - self.phi) / 1.5 + self.phi / (phi_ke0 + kf / kp * (2 / 3))) * kf
        B = 1 + 19.4 * (1 / N) ** 2
        Pe = 8.65 * B
        keG = rho_mol * cp_mol * u * self.ds / Pe
        ke = ke0 + keG

        phi_w = 0.0024 * N ** 1.58
        hw_0 = (2 * (1 - self.phi) + self.phi / (kf / kp * (2 / 3) + phi_w)) * kf / self.ds
        hw_g = (0.0835 * Re ** 0.91) * kf / self.ds
        hw = hw_0 + hw_g
        print(hw_0, hw_g)
        Ut = 1 / (1 / hw + self.De / 6.13 / ke)

        return Ut

    def htc4(self, T, P, F):
        """
        calculate the internal heat transfer coe in within the catalyst layer
        Cui, CEJ, 2020
        :param T: temperature of reactor gas, K
        :param P: pressure of reactor, bar
        :param F: molar flow rate of each component, mol/s; ndarray
        :return: convection heat transfer coefficient, W/m2 K
        """
        Ft = np.sum(F)
        y = F / Ft * P
        z = self.vle_cal.z(T=T, P=P, x=F)
        gas_prop = mixture_property(T, pd.Series(y, index=self.comp_list), P, z, rho_only=False)
        y = F / np.sum(F)
        v = np.sum(F) / self.vle_cal.cal_rho_mol(T, P, y)
        u = v / self.Ae
        rho_mass = self.vle_cal.cal_rho(T, P, y)
        rho_mol = self.vle_cal.cal_rho_mol(T, P, y)
        cp_mol = self.vle_cal.cal_cp(T, P, y)  # J/mol/K
        MW = self.vle_cal.mw(T, P, y)  # g/m

        kf = gas_prop['k']
        vis = gas_prop['vis']
        N = self.De / self.ds
        Re = self.ds * u * rho_mass / vis
        Pr = vis * (cp_mol / (MW / 1000)) / kf

        phi_ke0 = 0.13 * (1 - self.phi) ** 1.44
        kp = 0.21 + 0.00015 * T  # thermal conductivity of particle 0.38
        ke0 = ((1 - self.phi) / 1.5 + self.phi / (phi_ke0 + kf / kp * (2 / 3))) * kf
        B = 1 + 19.4 * (1 / N) ** 2
        Pe = 8.65 * B
        keG = rho_mol * cp_mol * u * self.ds / Pe
        ke = ke0 + keG

        Nu_w0 = (1.3 + 5 / N) * (ke0 / kf)
        Nu_ws = 0.3 * Pr ** (1 / 3) * Re ** 0.75
        Nu_m = 0.054 * Pr * Re
        Nu_w = Nu_w0 + 1 / (1 / Nu_ws + 1 / Nu_m)
        hw = Nu_w*kf/self.ds

        Ut = 1 / (1 / hw + self.De / 6.13 / ke)
        Bi = hw*self.De/ke
        # print(Re, ke, hw, Bi)
        # print(u*2/1e-6)
        return Ut
