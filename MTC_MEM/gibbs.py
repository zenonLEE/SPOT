from datetime import datetime

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from CoolProp.CoolProp import PropsSI
from scipy.optimize import minimize
from scipy.integrate import solve_ivp
from prop_calculator import VLEThermo

from insulator import Insulation
from thermo import (ChemicalConstantsPackage, SRKMIX, FlashVL, CEOSLiquid, CEOSGas, HeatCapacityGas,
                    GibbsExcessLiquid, MSRKMIX, SRKMIXTranslated, SRK)


class Gibbs(VLEThermo):
    def __init__(self, comp, ref='ch', elements=None):
        super().__init__(comp, ref)
        self.comp = comp
        # elements counts for each comps
        self.elements = self.__elements() if elements is None else elements
        try:
            self.const, self.cor = ChemicalConstantsPackage.from_IDs(self.comp.tolist())
        except AttributeError:
            self.const, self.cor = ChemicalConstantsPackage.from_IDs(self.comp)
        self.cp, self.eos_kw = self.eos_paras()

    @staticmethod
    def __elements():
        return np.array([[1, 0, 1, 0, 1],
                         [0, 2, 4, 2, 0],
                         [2, 0, 1, 1, 1]])

    def cal_Gr(self, T, P, sto):
        sto = np.array(sto)
        feed = np.where(sto <= 0, sto, 0)
        product = np.where(sto > 0, sto, 0)
        return self.cal_G(T, P, product) + self.cal_G(T, P, feed)

    def cal_Hr(self, T, P, sto):
        sto = np.array(sto)
        feed = np.where(sto <= 0, sto, 0)
        product = np.where(sto > 0, sto, 0)
        return self.cal_H(T, P, product) + self.cal_H(T, P, feed)

    def cal_Sr(self, T, P, sto):
        sto = np.array(sto)
        feed = np.where(sto <= 0, sto, 0)
        product = np.where(sto > 0, sto, 0)
        return self.cal_S(T, P, product) + self.cal_S(T, P, feed)

    def cal_dH(self, T, P, x):
        frac = x / np.sum(x)
        gas = CEOSGas(SRKMIX, HeatCapacityGases=self.cp, eos_kwargs=self.eos_kw, T=T, P=P, zs=frac)
        H_dep_in = gas.H_dep()
        liquid = CEOSLiquid(SRKMIX, HeatCapacityGases=self.cp, eos_kwargs=self.eos_kw)
        flasher = FlashVL(self.const, self.cor, liquid=liquid, gas=gas)
        PT = flasher.flash(zs=frac, T=T, P=P * 1E5)
        H_dep = PT.H()
        return (H_dep_in - H_dep) * np.sum(x) / 1000  # kW

    def cal_dn(self, T, P, x):
        frac = x / np.sum(x)
        gas = CEOSGas(SRKMIX, HeatCapacityGases=self.cp, eos_kwargs=self.eos_kw, T=T, P=P, zs=frac)
        liquid = CEOSLiquid(SRKMIX, HeatCapacityGases=self.cp, eos_kwargs=self.eos_kw)
        flasher = FlashVL(self.const, self.cor, liquid=liquid, gas=gas)
        PT = flasher.flash(zs=frac, T=T, P=P * 1E5)
        # gas_sep = PT.gas.zs
        liq_sep = PT.liquid0.zs

        return (1 - PT.VF) * np.sum(x) * np.array(liq_sep)  # mol/s

    @staticmethod
    def eq_cons(x, feed, element_count):
        """
        conservation of matter
        """
        return [np.sum((x - feed) * element_count[i]) for i in range(3)]

    def mix_gibbs(self, T, P, x):
        """
        calculate the gibbs energy of mixture at give T Pp
        :param T: K
        :param P: bar
        :param x: molar flow of comps in mixture, mol/s (np.array)
        :return: total gibbs energy of mixture, kW
        """
        x = np.array(x)
        if (x < 0).any():
            return 0
        else:
            return self.cal_G(T, P, x)

    def min_eq(self, T, P, feed, x0=None):
        """
        calculate the product by minimizing the Gibbs energy
        :param T: K
        :param P: bar
        :param feed: molar flow rate (pd.Series), mol/s
        :param x0: initial guess of product
        :return:
        """
        # Define the initial guess for the molar fractions of the product
        if x0 is None:
            x0 = feed.copy()
            r_guess, s_guess = 0.5, 0.5
            x0.loc[:] = np.array([(1 - r_guess) * x0[0],
                                  x0[1] - x0[0] * (2 * r_guess * s_guess + r_guess),
                                  r_guess * s_guess * x0[0],
                                  r_guess * x0[0],
                                  r_guess * (1 - s_guess) * x0[0]])
        else:
            x0 = x0
        # Combine equality and inequality constraints
        constraint_rule = [{'type': 'eq', 'fun': lambda x: self.eq_cons(x, feed, self.elements)},
                           {'type': 'ineq', 'fun': lambda x: x}]  # ensure >= 0

        # Solve the optimization problem using SciPy's minimize function
        res = minimize(lambda x: self.mix_gibbs(T, P, x), x0, constraints=constraint_rule,
                       method="SLSQP")
        return res.x, res.fun

    def solve_over_range(self, T, P, feed_comp, save=False):
        """
        calculate the equilibrium conversion over ranges of temperature and pressure
        :param T: K
        :param P: bar
        :param feed_comp: molar flow rate (pd.Series), mol/s
        :param save: whether save data to file
        """
        CO2_R = pd.DataFrame(index=T, columns=P)
        CO_R = pd.DataFrame(index=T, columns=P)
        C_R = pd.DataFrame(index=T, columns=P)
        select = pd.DataFrame(index=T, columns=P)
        for P_eq in P:
            for T_eq in T:
                print(T_eq, P_eq)
                # Calculate equilibrium
                product, min_gibbs_energy = self.min_eq(T_eq, P_eq, feed_comp)
                CO2_R.loc[T_eq, P_eq] = (feed_comp[0] - product[0]) / feed_comp[0]
                CO_R.loc[T_eq, P_eq] = (feed_comp[-1] - product[-1]) / feed_comp[-1]
                C_R.loc[T_eq, P_eq] = (feed_comp[-1] - product[-1]+ feed_comp[0] - product[0]) / (feed_comp[-1]+feed_comp[0])
                print(CO2_R.loc[T_eq, P_eq],CO_R.loc[T_eq, P_eq],C_R.loc[T_eq, P_eq])
                select.loc[T_eq, P_eq] = (product[2] - feed_comp[2]) / (feed_comp[0] - product[0])
        if save:
            res_path = 'eq_SRK_%s_%s_%s_%s.xlsx' % (min(T), max(T), min(P), max(P))
            with pd.ExcelWriter(res_path, engine='openpyxl') as writer:
                CO2_R.to_excel(writer, index=True, header=True, sheet_name='conversion_CO2')
                CO_R.to_excel(writer, index=True, header=True, sheet_name='conversion_CO')
                C_R.to_excel(writer, index=True, header=True, sheet_name='conversion_C')
                select.to_excel(writer, index=True, header=True, sheet_name='select')
        return CO2_R, select

    def series_reactor(self, feed_cond, r_target, sp_paras, sp_type):
        """
        the combination of Gibbs reactor and separator
        :param feed_cond: dict(comp, Tin, Pin) molar flow rate (pd.Series), mol/s
        :param r_target: the target conversion of CO2 for whole process
        :param sp_paras: separator paras
        :param sp_type: type of separator, flash, diff
        :return: product in each stage
        """

        r_t = 0
        products, r_each = [], []
        reactor_feed = feed_cond['comp']
        Tin, Pin = feed_cond['Tin'], feed_cond["Pin"]
        sep_work, dHrs, Ws = [], [], []
        Q, E = [], []
        sep_liq, sep_gas = [], []  # separated product, separated unreacted gas

        while r_t < r_target:
            # print(reactor_feed.values)
            product, _ = self.min_eq(Tin, Pin, reactor_feed)
            # metric of single reactor
            dHr = self.cal_H(Tin, Pin, product) - self.cal_H(Tin, Pin, reactor_feed.values)  # dH during reaction, kW
            dHrs.append(dHr)
            r = (reactor_feed[0] - product[0]) / reactor_feed[0]  # CO2 conversion
            r_each.append([r])

            # generate the feed for the next stage through separator
            sp_feed_paras = [Tin, Pin, product]
            sep_res = self.flash_sep(sp_feed_paras, sp_paras) if sp_type == 'flash' \
                else self.cond_sep(sp_feed_paras, sp_paras)
            reactor_feed = pd.Series(sep_res['gas'], index=self.comp)
            # reactor_feed, Q_diff, prod_CH4O = self.cond_sep(sp_feed_paras, sp_paras)

            # metric of separator
            Q.append([sep_res['H_c'], sep_res['H_h']])
            E.append([sep_res['E_c'], sep_res['E_h']])
            sep_liq.append(sep_res['liq'])
            sep_gas.append(sep_res['gas'])

            # metric of the whole process
            r_t = (feed_cond['comp'][0] - product[0]) / feed_cond['comp'][0]  # total conversion of CO2
            # selectivity of CH4O
            s_t = (feed_cond['comp'][0] - product[0] - product[-1]) / (feed_cond['comp'][0] - product[0])

            products.append(product.tolist())

        process_metric = pd.Series([s_t, r_t], index=['s', 'r'])
        sep_liq = np.array(sep_liq)  # .reshape(len(prods_CH4O), 1)
        sep_gas = np.array(sep_gas)
        dHrs = np.array(dHrs).reshape(len(dHrs), 1)
        Q, E = np.array(Q), np.array(E)
        res = np.hstack((dHrs, Q, E, sep_liq, sep_gas))

        # metric of each reactor
        sim_metric = pd.DataFrame(res, columns=['Hrs', 'H_c', 'H_h', 'E_c', 'E_h'] +
                                               [i + '_l' for i in self.comp.tolist()] +
                                               [i + '_g' for i in self.comp.tolist()])
        # metric of the whole process
        process_metric = self.series_metric(feed_cond, sp_paras, sim_metric, process_metric)
        return process_metric

    @staticmethod
    def save_res(sim_metric, feed_cond, sp_cond, sp_type):
        res_path = f'res_Gibbs/{sp_type}_{datetime.now().date()}.xlsx'

        res_save = pd.concat((pd.Series(feed_cond)[['Tin', 'Pin']], pd.Series(sp_cond), sim_metric))
        res_save = pd.DataFrame(res_save.values.reshape(1, len(res_save)), columns=res_save.index.tolist())

        try:
            with pd.ExcelWriter(res_path, engine='openpyxl', mode='a', if_sheet_exists="overlay") as writer:
                try:
                    res_saved = pd.read_excel(res_path)
                    res_save = pd.concat([res_saved, res_save], ignore_index=True)
                    res_save.to_excel(writer, index=False, header=True)
                except ValueError:
                    res_save.to_excel(writer, index=False, header=True)
        except FileNotFoundError:
            res_save.to_excel(res_path, index=False, header=True)

    @staticmethod
    def energy_analysis(in_ma, prod_ma, left_ma, in_q, out_q, react_q):
        """
        perform enthalpy and exergy analysis of system
        :param in_ma: input material, dict(T=Tin, Pp=Pin, x=feed)
        :param prod_ma: product material
        :param left_ma: unreacted material
        :param in_q: input heat, [enthalpy, exergy]
        :param out_q: output heat
        :param react_q: reaction heat, float, minus for exothermic
        :return:
        """
        metric = pd.Series(index=['H_in', 'E_in', 'H_o', 'E_o', 'H_l', 'E_l',
                                  'H_qin', 'E_qin', 'H_qo', 'E_qo', 'H_r', "E_r"])
        eos = VLEThermo(np.array(["CO2", "H2", "Methanol", "H2O", "carbon monoxide"]), ref='ev')

        metric.loc[['H_o', 'E_o']] = 0
        for prod in prod_ma:
            metric.loc[['H_o', 'E_o']] += np.array([eos.cal_H(**prod), eos.cal_E(**prod)])
        metric.loc[['H_in', 'E_in']] = [eos.cal_H(**in_ma), eos.cal_E(**in_ma)]
        metric.loc[['H_l', 'E_l']] = [eos.cal_H(**left_ma), eos.cal_E(**left_ma)]
        metric.loc[['H_qin', 'E_qin', 'H_qo', 'E_qo']] = in_q + out_q
        metric.loc[['H_r', "E_r"]] = [react_q, react_q * (1 - 298.15 / in_ma['T'])]
        # thermodynamic efficiency
        metric['eta_e'] = (metric.loc['E_o'] + metric.loc['E_l'] + metric.loc['E_qo'] - metric.loc["E_r"]) / \
                          (metric.loc['E_in'] + metric.loc['E_qin'])
        metric['H_dev'] = ((metric.loc['H_o'] + metric.loc['H_l'] + metric.loc['H_qo'] - metric.loc['H_r']) -
                           (metric.loc['H_in'] + metric.loc['H_qin'])) / (metric.loc['H_in'] + metric.loc['H_qin'])
        return metric

    def series_metric(self, feed_paras, sep_paras, stage_res, process_res):
        """
        calculate the metric of series reactor with flash sep
        :param feed_paras: [Tin, Pin, feed] (dict)
        :param sep_paras: condensation temperature (dict)
        :param stage_res: sim res of each stage(pd.Dataframe)
        :param process_res: conversion and selectivity of the whole process (pd.Series)
        :return:
        """
        yield_ch3oh = pd.Series(stage_res["Methanol_l"].sum(), index=['y_ch3oh'])
        r_CH3OH_H2O = pd.Series(stage_res["Methanol_l"].sum() / stage_res["H2O_l"].sum(), index=['r_CH3OH_H2O'])
        # metric = pd.concat((energy_metric, yield_ch3oh, r_CH3OH_H2O, process_res))

        # calculate the exergy
        # generate the info of each material flow, including T,Pp,molar flow
        Tin, Pin, feed = feed_paras['Tin'], feed_paras['Pin'], feed_paras['comp'].to_numpy()
        Tc = sep_paras['Tc']
        feed = dict(T=Tin, P=Pin, x=feed)
        p_left = dict(T=Tin, P=Pin, x=stage_res.iloc[-1].filter(like='_g').to_numpy())
        p_cond = []
        for i in stage_res.filter(like='_l').to_numpy():
            p_cond.append(dict(T=Tc, P=Pin, x=i))
            # p_cond.append(dict(T=Tc, Pp=Pin, x=[0,0,i[2],i[3],0]))
        # p_cond = [dict(T=Tc, Pp=Pin, x=np.array([0, 0, stage_res["Methanol_l"].sum(), stage_res["H2O_l"].sum(), 0]))]
        # np.array([0, 0, stage_res["Methanol_l"].sum(), stage_res["H2O_l"].sum(), 0])

        # perform energy analysis
        reaction_energy = stage_res['Hrs'].sum()
        energy_metric = self.energy_analysis(feed, p_cond, p_left,
                                             stage_res[['H_h', 'E_h']].sum().abs().tolist(),
                                             stage_res[['H_c', 'E_c']].sum().abs().tolist(),
                                             reaction_energy)
        metric = pd.concat((energy_metric, yield_ch3oh, r_CH3OH_H2O, process_res))
        return metric

    def ideal_sep(self, feed_paras, sep_paras):
        """
        perform separation through a flash can
        :param feed_paras: [Tin, Pin, feed] (list)
        :param sep_paras: separation ratio (dict)
        """
        [Tin, Pin, feed] = feed_paras
        sp = sep_paras['sp']
        prod_out = feed.copy() * sp
        prod_out[[0, 1, -1]] = 0
        feed_out = feed - prod_out

        # metric of sep
        # released energy during sep, kW
        H = self.cal_H(Tin, Pin, prod_out) + self.cal_H(Tin, Pin, feed_out) - self.cal_H(Tin, Pin, feed)
        # released exergy during sep, kW
        E = self.cal_E(Tin, Pin, prod_out) + self.cal_E(Tin, Pin, feed_out) - self.cal_E(Tin, Pin, feed)

        res = dict(gas=feed_out, liq=prod_out, H_c=H, E_c=E)
        return res

    def flash_sep(self, feed_paras, sep_paras):
        """
        perform separation through a flash can
        :param feed_paras: [Tin, Pin, feed] (list)
        :param sep_paras: condensation temperature (dict)
        """
        [Tin, Pin, feed] = feed_paras
        Tc = sep_paras['Tc']
        gas, liq, vf = self.flash(Tc, Pin, feed)
        gas_out = np.array(gas) * vf * np.sum(feed)
        liq_out = np.array(liq) * (1 - vf) * np.sum(feed)

        # metric of sep
        H_c = self.cal_H(Tc, Pin, feed) - self.cal_H(Tin, Pin, feed)  # released energy during sep, kW
        E_c = self.cal_E(Tc, Pin, feed) - self.cal_E(Tin, Pin, feed)  # released exergy during sep, kW
        H_h = self.cal_H(Tin, Pin, gas_out) - self.cal_H(Tc, Pin, gas_out)  # required energy during sep, kW
        E_h = self.cal_E(Tin, Pin, gas_out) - self.cal_E(Tc, Pin, gas_out)  # required exergy during sep, kW

        res = dict(gas=gas_out, liq=liq_out, H_c=H_c, E_c=E_c, E_h=E_h, H_h=H_h)
        return res

    def cond_sep(self, feed_paras, sep_paras):
        """
        condenser separator
        :param feed_paras: [T K, Pp bar, feed_flow mol/s]
        :param sep_paras: separator paras location=0, Din=Din, thick=Dd, Tc=Tc, heater=0
        :return:
        """
        [Tin, Pin, feed] = feed_paras
        Tc = sep_paras['Tc']

        res = self.diffusion(feed_paras, **sep_paras)  # W/m mol/s m
        m_diff = res['mflux']
        q_sen, q_lat = res['hflux'], res['hlg']
        l = min(feed[2] / m_diff[2], feed[3] / m_diff[3], key=abs)

        N_CH3OH = l * m_diff[2]  # separated CH3OH, mol/s
        N_H2O = l * m_diff[3]  # separated H2O, mol/s
        # print(N_H2O, N_CH3OH)
        # diffused sensible heat, kW
        # r_F =

        Q_sen = q_sen / m_diff[2] * N_CH3OH / 1000  # * (np.sum(feed[[0, 1, -1]] / np.sum(feed)))
        # Q_lat = q_lat / m_diff[2] * N_CH3OH / 1000
        Q_lat = (PropsSI('HMOLAR', 'T', Tc, 'Q', 1, 'water') -
                 PropsSI('HMOLAR', 'T', Tc, 'Q', 0, 'water')) * N_H2O / 1000 + \
                (PropsSI('HMOLAR', 'T', Tc, 'Q', 1, 'Methanol') -
                 PropsSI('HMOLAR', 'T', Tc, 'Q', 0, 'Methanol')) * N_CH3OH / 1000  # kW

        # print(Q_sen)
        # print(Q_lat)
        # diffused latent heat, kW
        liq_out = np.array([0, 0, N_CH3OH, N_H2O, 0])
        gas_out = feed - liq_out
        gas_out[np.abs(gas_out) < 1e-10] = 0
        Q_mix = self.cal_H(Tc, Pin, liq_out) - \
                self.cal_H(Tc, Pin, np.array([0, 0, N_CH3OH, 0, 0])) - \
                self.cal_H(Tc, Pin, np.array([0, 0, 0, N_H2O, 0]))
        # print(Q_lat, Q_mix)
        # Q_lat = Q_lat - Q_mix
        # metric of sep

        H_c = Q_lat + Q_sen  # released energy during sep, kW
        E_c = (Q_lat + Q_sen) * (1 - 298.15 / Tc)  # released exergy during sep, kW

        # eos = VLEThermo(comp=self.comp, ref='ev')
        # H_h = eos.cal_H(Tin, Pin, feed) - eos.cal_H(Tin, Pin, gas_out) - eos.cal_H(Tc, Pin, liq_out) - Q_lat
        # H_h = eos.cal_H(Tin, Pin, liq_out) - eos.cal_H(Tc, Pin, liq_out) - Q_lat

        # print(eos.cal_H(Tin, Pin, Ff), eos.cal_H(Tin, Pin, Fr), eos.cal_H(Tc, Pin, Fp), Q_lat)
        # print(H_h, Q_sen)
        # Q_sen_p = self.cal_H(Tin, Pin, Fp) - self.cal_H(Tc, Pin, Fp) - Q_lat
        # print(Q_sen_p)
        H_h = Q_sen  # required energy during sep, kW
        E_h = H_h * (1 - 298.15 / Tin)  # required exergy during sep, kW

        res = dict(gas=gas_out, liq=liq_out, H_c=H_c, E_c=E_c, E_h=E_h, H_h=H_h)

        if np.any(liq_out) >= 0:
            return res
        else:
            raise ValueError('Gas flow should be positive!')

    def cond_sep_along(self, feed_paras, sep_paras, r_target):
        """
        condenser separator
        :param feed_paras: [T K, Pp bar, feed_flow mol/s]
        :param sep_paras: separator paras location=0, Din=Din, thick=Dd, Tc=Tc, heater=0
        :return:
        """
        feed_comp = feed_paras[2]  # mol/s
        reactor_feed = pd.Series(feed_comp.copy(), index=self.comp)
        q, n = self.diffusion(feed_paras, **sep_paras)  # W/m mol/s m
        L0 = feed_comp[2] * 0.9 / n[2]  # separated CH3OH, mol/s
        r_sim_CH3OH, r_sim_H2O = 0, 0
        L = abs(L0)
        while r_sim_CH3OH < r_target and r_sim_H2O <= 1:
            res_along = self.diff_along(feed_paras, L, **sep_paras)
            # print(res_along[:, -1])
            r_sim_CH3OH = (res_along[3, 0] - res_along[3, -1]) / res_along[3, 0]  # separation ratio of CH3OH
            r_sim_H2O = (res_along[4, 0] - res_along[4, -1]) / res_along[4, 0]  # separation ratio of H2O
            Q_diff = res_along[-2, -1] / 1000  # kW

            if r_sim_CH3OH < 0.6:
                L += 0.5
            elif 0.6 < r_sim_CH3OH < 0.8:
                L += 0.4
            elif 0.8 < r_sim_CH3OH < 0.9:
                L += 0.3
            else:
                L += 0.1
        reactor_feed.loc[:] = res_along[1:-2, -1]
        return reactor_feed, Q_diff

    def sep(self, Fin, Tin, Pin, Tos, Pos, sf=1, Fos=None):
        """
        calculate the metric of separator
        F1 is split into F2 and F3
        :param F1: molar flow rate (pd.Series), mol/s
        :param T: K
        :param Pp: bar
        :return:
        """
        if Fos is None:
            Fo_prod = pd.Series(0, index=Fin.index)
            Fo_prod.loc['Methanol'] = sf * Fin.loc['Methanol']
            Fo_prod.loc['H2O'] = min(sf * 1.25, 1) * Fin.loc['H2O']
            Fo_feed = Fin - Fo_prod
        Ws = self.cal_G(Tin, Pin, Fin.values) - \
             self.cal_G(Tos[0], Pos[0], Fo_feed.values) - self.cal_G(Tos[1], Pos[1], Fo_prod.values)
        Win = self.cal_G(Tin, Pin, Fin.values) - \
              self.cal_G(Tos[0], Pos[0], Fo_feed.values) - self.cal_G(Tos[1], Pos[1], Fo_prod.values)

    @staticmethod
    def diffusion(feed, location, Din, thick, Tc):
        """

        :param feed: [T0, P0, F0: mol/s(ndarray)]
        :param location:
        :param Din:
        :param thick:
        :param Tc:
        :param heater:
        :return:
        """
        [T0, P0, F0] = feed
        Do = Din + thick * 2
        # property_feed = mixture_property(T0, xi_gas=pd.Series(F0 / np.sum(F0), index=subs), Pt=P0)
        insula_sim = Insulation(Do, Din, 1, location)
        res = insula_sim.flux(T0, P0, F0, Tc)  # mol/(s m) W/m
        h_diff, m_diff = res['hflux'], res['mflux']
        # r_h_m = h_diff / 1e3 / (m_diff[2])  # kJ/mol CH4O
        return res  # h_diff, m_diff  # W/m, mol/(s m)

    @staticmethod
    def diff_along(feed, L, location, Din, thick, Tc, heater):

        [T0, P0, F0] = feed
        Do = Din + thick * 2
        insula_sim = Insulation(Do, Din, 1, location)

        def model(z, y):
            # y= [F_CO2, F_H2, F_CH3OH, F_H2O, F_CO
            # Tr]
            F = np.array(y[:5])
            Tr = y[-1]
            res_diff = insula_sim.flux(Tr, P0, F, Tc)  # mol/(s m) W/m
            dF_dz = res_diff["mflux"]
            dq_dz = res_diff['hflux']
            dTr_dz = 0  # res["Tvar"] + q_heater / heat_cap  # res_react['tc']
            res_dz = np.hstack((dF_dz, np.array([dq_dz, dTr_dz])))
            return res_dz

        # property_feed = mixture_property(T0, xi_gas=pd.Series(F0 / np.sum(F0), index=subs), Pt=P0)
        z_span = [0, L]
        ic = np.hstack((F0, np.array([0, T0])))
        res_sim = solve_ivp(model, z_span, ic, method='BDF', t_eval=np.linspace(0, L, 1000))  # LSODA BDF
        res = np.vstack((np.linspace(0, L, 1000), res_sim.y))
        return res


def find_best_cond(feed, T_range, Din_range, Dd_range):
    """
    find the lowest energy consumption with different conditions
    :param feed: molar flow rate (pd.Series), mol/s
    :param T_range:
    :param Din_range:
    :param Dd_range:
    :return:
    """
    sims_res = pd.DataFrame(columns=['Tc', 'Din', 'Dd', 'r', 's', 'Q', 'rq_in_diff'],
                            index=np.zeros(len(T_range) * len(Din_range) * len(Dd_range)))

    i = 0
    for Tc in T_range:
        for Din in Din_range:
            for Dd in Dd_range:
                diffusor_para = dict(location=0, Din=Din, thick=Dd, Tc=Tc, heater=0)
                gibbs_cal = Gibbs(feed.index)
                res = gibbs_cal.series_reactor(feed, 503, 70, r_target=0.95, sp_paras=diffusor_para)
                sims_res.iloc[i] = np.hstack((np.array([Tc, Din, Dd]), res))
                print(sims_res.iloc[i])
                i += 1
    path = f'res_Gibbs/eq_diff_{min(T_range)}_{max(T_range)}_' \
           f'{min(Din_range):.2f}_{max(Din_range):.2f}_' \
           f'{min(Dd_range):.3f}_{max(Dd_range):.3f}.xlsx'
    sims_res.to_excel(path, index=True, header=True, sheet_name='conversion')


def metric_single(feed, r, sp_para):
    gibbs_cal = Gibbs(feed['comp'].index, ref='ch')
    res = gibbs_cal.series_reactor(feed, r_target=r, sp_paras=sp_para, sp_type='diff')
    return res


if __name__ == '__main__':
    in_comp = pd.Series([7.18E-04,1.64E-04,0,0,	4.57E-04], index=["CO2", "H2", "Methanol", "H2O", "carbon monoxide"])
    Ts = np.arange(483, 553, 10)
    Ps = np.arange(30,80,10)
    gibbs_cal = Gibbs(in_comp.index)
    gibbs_cal.solve_over_range(Ts,Ps,in_comp,True)
    # 0.008154456, 0.024463369, 0, 0, 0
    # 0.00522348 0.01617213 0.00268013 0.00293098 0.00025085
    # Tcs = np.arange(323, 403, 5)
    Dins = np.arange(0.04, 0.09, 0.01)
    thick = np.arange(0.002, 0.014, 0.002)

    # find_best_cond(in_gas, Tcs, Dins, thick)
    # in_gas = dict(comp=in_comp, Tin=503, Pin=70)
    # # sp_para = dict(Tc=353)  # Din=0.08, Dd=0.01
    # sp_para = dict(Tc=353, location=0, Din=0.05, thick=0.005)
    # sim_res = metric_single(in_gas, r=0.95, sp_para=sp_para)
    # print(sim_res)
    # Gibbs.save_res(sim_res, in_gas, sp_para, 'diff')

    # metric_single()
    # Tins, Pins = np.arange(483, 513, 5), np.arange(30, 80, 10)
    # Tcs = np.arange(323, 393, 10)
    # for Tc in Tcs:
    #     for Din in Dins:
    #         for Dd in thick:
    #             try:
    #                 in_gas = dict(comp=in_comp, Tin=503, Pin=70)
    #                 # sp_para = dict(Tc=Tc)  # Din=0.08, Dd=0.01
    #                 sp_para = dict(Tc=Tc, location=1, Din=Din, thick=Dd)
    #                 sim_res = metric_single(in_gas, r=0.95, sp_para=sp_para)
    #                 print(sim_res)
    #                 # Gibbs.save_res(sim_res, in_gas, sp_para, 'diff')
    #             except AttributeError:
    #                 pass

    # calculate the metric at std
    # cal = VLEThermo(comp=["CO2", "H2", "Methanol", "H2O", "carbon monoxide"], ref='ev')
    # feed = np.array([1, 3, 0, 0, 0])
    # product = np.array([0, 0, 1, 1, 0])
    # T, Pp = 503, 70
    # Hr = cal.cal_H(T, Pp, product) - cal.cal_H(T, Pp, feed)
    # Er = cal.cal_E(T, Pp, product) - cal.cal_E(T, Pp, feed)
    # eta_E = cal.cal_E(T, Pp, product)/cal.cal_E(T, Pp, feed)
    # print(Hr, Er, eta_E)
