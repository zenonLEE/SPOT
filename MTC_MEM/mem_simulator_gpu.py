import os.path
import time

import scipy

from prop_calculator import VLEThermo
from CoolProp.CoolProp import PropsSI
from utility import heater, compressor, flasher, spliter, HeatExchanger, mixer, mult_compressor, \
    adsorber, valve, Distiller, DistillerOpt, multi_comp_opt, HX4Water, heater_duty
import numpy as np
import pandas as pd
from datetime import datetime
from reactor import Reaction
from reactor_HMT import HMT
import torch
from torchdiffeq import odeint


# from gibbs import Gibbs
RED = '\033[91m'
ENDC = '\033[0m'


class Simulation:
    def __init__(self, reactors_para, chem_para, feed_para, mem_para, eos, drop):

        # basic info
        self.comp_list = ["CO2", "H2", "Methanol", "H2O", "CO", 'N2']
        self.R = 8.314
        self.eos = eos  # 0 for ideal 1 for SRK
        self.drop = drop  # 1 for ergun 0 for zero drop

        # reactor insulator feed para
        self.reactors_para = reactors_para  #

        self.reactor_para = self._reactor_para()  # rewrite the reactor para
        self.chem_para = chem_para
        self.chem_para = chem_para
        self.feed_para = feed_para
        self.feed0_para = self._feed_para()
        self.Fr0 = self.feed0_para['F0']
        F_in = self.feed0_para.loc['F0']
        self.comp_num = len(F_in)
        # if len(F_in) == 6 and F_in[-1] == 0:
        #     self.F0 = self.F0[:-1]
        #     self.feed0_para.loc['F0'] = F_in[:-1]
        #     self.comp_num -= 1
        #     self.comp_list = self.comp_list[:-1]
        self.pro_cal = VLEThermo(self.comp_list, ref='ev')
        self.membranes = mem_para
        self.membrane = self._mem_para()
        self.cond_record = pd.DataFrame(index=self.comp_list + ["T", "P"], dtype='float64')
        self.utility_record = pd.DataFrame(index=['Q', 'E', 'W'], dtype='float64')
        self.rf_record = None
        self.iso_record = None

    def _reactor_para(self):
        """
        generate para for each reactor in series (pd.Dataframe)
        [[stage1, paras...], [stage2, paras...]]
        """
        self.stage = self.reactors_para['stage']
        self.Uc = self.reactors_para['Uc']
        self.recycle = self.reactors_para['recycle']
        Dt_name = [f'Dt{n + 1}' for n in range(self.stage)]
        L_name = [f'L{n + 1}' for n in range(self.stage)]
        Din_name = [f'Dc{n + 1}' for n in range(self.stage)]
        self.L, self.Dt = self.reactors_para[L_name], self.reactors_para[Dt_name]
        self.Din = self.reactors_para[Din_name]
        self.nrt = self.reactors_para['nrt']
        reactors = pd.DataFrame(index=np.arange(self.stage),
                                columns=['L', 'Dt', 'Dc'] + self.reactors_para.index[self.stage * 3:].tolist())
        for n in range(self.stage):
            reactors.loc[n, ['L', 'Dt', 'Dc']] = [self.L[n], self.Dt[n], self.Din[n]]
            reactors.iloc[n, 3:] = self.reactors_para[3 * self.stage:]
        return reactors

    def _feed_para(self):
        """
        generate feed para
        """
        self.P0, self.T0 = self.feed_para["P"], self.feed_para["T"]  # P0 bar, T0 K
        # self.T_feed = self.feed_para["T_feed"]
        try:
            self.T0_CO2, self.T0_H2 = self.feed_para["T0_CO2"], self.feed_para["T0_H2"]
            self.P0_CO2, self.P0_H2 = self.feed_para["P0_CO2"], self.feed_para["P0_H2"]
        except KeyError:
            self.T0_CO2, self.T0_H2 = self.T0, self.T0
            self.P0_CO2, self.P0_H2 = self.P0, self.P0
        self.vle_cal = VLEThermo(self.comp_list)

        if self.feed_para["fresh"] == 1:  # the feed to the plant is fresh stream
            self.r_inert = self.feed_para['inert']
            self.F0 = np.zeros(len(self.comp_list))  # component of feed gas, mol/s; ndarray
            # volumetric flux per tube from space velocity
            if self.feed_para["H2"] == 0:
                self.sv = self.feed_para["Sv"]
                # volumetric flux per tube under input temperature and pressure, m3/s
                self.v0 = self.sv * self.L[0] * np.pi * self.Dt[0] ** 2 / 4 / 3600 / self.nrt
                self.Ft0 = self.P0 * 1e5 * self.v0 / self.R / self.T0  # total flux of feed,mol/s
                self.F0[0] = 1 / (1 + 1 * self.feed_para["H2/CO2"] + self.feed_para['CO/CO2']) * \
                             self.Ft0 * (1 - self.r_inert)
                self.F0[4] = self.F0[0] * self.feed_para['CO/CO2']
                self.F0[1] = self.Ft0 * (1 - self.r_inert) - self.F0[0] - self.F0[4]
                x0 = self.F0 / self.Ft0
                z0 = self.vle_cal.z(T=self.T0, P=self.P0, x=self.F0)
                self.Ft0 = self.Ft0 / z0
                self.F0 = x0 * self.Ft0
                self.H2 = self.F0[1] * 8.314 * 273.15 / 1e5 * 3600  # Nm3/h
            else:
                self.H2 = self.feed_para["H2"]  # Nm3/h
                self.F0[1] = self.H2 / 3600 * 1e5 / self.R / 273.15  # mol/s
                self.F0[0] = self.F0[1] / self.feed_para["H2/CO2"]
                self.F0[4] = self.F0[0] * self.feed_para['CO/CO2']
                self.F0[-1] = np.sum(self.F0[:-1]) * self.r_inert / (1 - self.r_inert)
                self.Ft0 = np.sum(self.F0)
                z0 = self.vle_cal.z(T=self.T0, P=self.P0, x=self.F0)
                self.v0 = self.Ft0 * self.R * self.T0 * z0 / (self.P0 * 1e5)
                self.sv = self.v0 * self.nrt * 3600 * 4 / self.L[0] / np.pi / self.Dt[0] ** 2
            # print(self.sv, self.H2)
        else:
            self.F0 = self.feed_para[self.comp_list].to_numpy()
            self.Ft0 = np.sum(self.F0)
            z0 = self.vle_cal.z(T=self.T0, P=self.P0, x=self.F0)
            self.v0 = self.Ft0 * self.R * self.T0 * z0 / (self.P0 * 1e5)
            self.sv = self.v0 * self.nrt * 3600 * 4 / self.L[0] / np.pi / self.Dt[0] ** 2
            self.H2 = self.F0[1] * self.R * 273.15 / 1E5
        feed = pd.Series([self.T0, self.P0, self.F0], index=['T0', 'P0', 'F0'])
        return feed

    def _mem_para(self):
        """
        reconstruct insulator para
        :return: [[stage1, paras...], [stage2, paras...]] pd.Dataframe
        """
        paras_stage = ['Din', 'Thick', 'Tc', 'qmc', 'Th', 'qmh', 'Pp', 'per_H2O', 'per_S']
        paras_array_name = {'status_name': [f'status{n + 1}' for n in range(self.stage)],
                            'pattern_name': [f'pattern{n + 1}' for n in range(self.stage)],
                            'pattern_h_name': [f'pattern_h{n + 1}' for n in range(self.stage)]}
        # pattern -1 mean counter-flow
        for n in range(self.stage):
            for para_stage in paras_stage:
                paras_array_name[f'{para_stage}_name'] = [f'{para_stage}{n + 1}' for n in range(self.stage)]

        self.status = self.membranes[paras_array_name['status_name']].values
        self.pattern = self.membranes[paras_array_name['pattern_name']].values
        self.pattern_h = self.membranes[paras_array_name['pattern_h_name']].values
        self.Tc = self.membranes[paras_array_name['Tc_name']].values
        self.qmc = self.membranes[paras_array_name['qmc_name']].values
        self.Th = self.membranes[paras_array_name['Th_name']].values
        self.qmh = self.membranes[paras_array_name['qmh_name']].values
        self.Thick = self.membranes[paras_array_name['Thick_name']].values
        self.Pp = self.membranes[paras_array_name['Pp_name']].values
        self.Per_H2O = self.membranes[paras_array_name['per_H2O_name']].values
        self.Per_S = self.membranes[paras_array_name['per_S_name']].values

        self.per = self.membranes['per']
        self.nit = self.membranes["nit"]  # tube number of the insulator
        self.location = self.membranes["io"]
        self.Fp = self.membranes["Fp"]
        self.Fp = self.Fp * self.Fr0

        if self.location == 0:  # reaction occurs in the tube side
            self.Din = self.Dt.values  # self.insulator_para['Din']
            # self.heater = max(self.heater, (523 - 333) / self.insulator_para['Thick'] * 0.2 * np.pi * self.Dt / 3)
        else:
            self.Din = self.membranes[paras_array_name['Din_name']].values

        self.Do = self.Din + self.membranes[paras_array_name['Thick_name']].values * 2
        membranes = pd.DataFrame(index=np.arange(self.stage),
                                 columns=['status', 'pattern', 'pattern_h'] +
                                         paras_stage + ['nit', 'location', 'Fp', 'per'])

        for n in range(self.stage):
            membranes.loc[n, ['status', 'pattern', 'pattern_h'] + paras_stage] = [self.status[n],
                                                                                  self.pattern[n], self.pattern_h[n],
                                                                                  self.Din[n],
                                                                                  self.Thick[n], self.Tc[n],
                                                                                  self.qmc[n],
                                                                                  self.Th[n], self.qmh[n],
                                                                                  self.Pp[n], self.Per_H2O[n],
                                                                                  self.Per_S[n]]
            membranes.iloc[n, -4:] = [self.nit, self.location, self.Fp, self.per]
        # pd.set_option('display.max_columns', None)
        # print(membranes)

        return membranes

    def reactor_metric(self, sim_res):
        """
        save the one-pass performance of reactor
        :param sim_res: molar flux of each component along the reactor
        """

        # reactor metric
        # y= z [F_CO2, F_H2, F_CH3OH, F_H2O, F_CO, F_Ar
        # [molar fraction]
        # react: F_CO2, F_H2, F_CH3OH, F_H2O, F_CO, F_Ar
        # diff: F_CO2, F_H2, F_CH3OH, F_H2O, F_CO, F_Ar
        # permeate: F_CO2, F_H2, F_CH3OH, F_H2O, F_CO, F_Ar
        # Tr, Tc, Th, Pp, q_react, q_diff, q_heater, h_diff]

        To_r = sim_res[self.comp_num * 5 + 1, -1]  # reactor output temperature
        r = (sim_res[1, 0] - sim_res[1, -1]) / sim_res[1, 0]  # CO2 conversion ratio
        # r_H2 = (sim_res[2, 0] - sim_res[2, -1]) / sim_res[2, 0]  # H2 conversion ratio
        r_c = (sim_res[1, 0] - sim_res[1, -1] + sim_res[5, 0] - sim_res[5, -1]) / (
                sim_res[1, 0] + sim_res[5, 0])  # conversion ratio of Carbon
        dF_react_rwgs = sim_res[5][-1] - sim_res[5][0]  # amount of reaction CO2 to CO
        dF_react_ch3oh = (sim_res[1, 0] - sim_res[1][-1]) - dF_react_rwgs  # amount of reaction CO2 to CH3OH
        dF_h2 = (sim_res[2, 0] - sim_res[2, -1])
        dF_react_h2 = -dF_react_ch3oh * 3 - dF_react_rwgs
        r_H2 = -dF_react_h2 / sim_res[2, 0]
        dF_react_h2o = dF_react_rwgs + dF_react_ch3oh  # amount of water produced
        s_react = dF_react_ch3oh / (sim_res[1, 0] - sim_res[1, -1])  # selectivity of reactions
        dH_react = sim_res[self.comp_num * 5 + 5, -1]
        dP = sim_res[self.comp_num * 5 + 4, -1] - sim_res[self.comp_num * 5 + 4, 0]

        # in-situ separation metric

        # dF_diff_ch3oh = dF_react_ch3oh - (sim_res[3, -1] - sim_res[3, 0])  # amount of CH3OH sp
        dF_diff_h2o = dF_react_h2o - (sim_res[4][-1] - sim_res[4][0])  # amount of H2O sp
        dF_diff_h2 = dF_h2 - (-dF_react_h2)
        # sp_ch3oh = dF_diff_ch3oh / dF_react_ch3oh  # separation ratio of CH3OH
        sp_h2o = dF_diff_h2o / dF_react_h2o  # separation ratio of H2O
        sp_h2 = dF_diff_h2 / sim_res[2, 0]
        # N_CH3OH_H2O = dF_diff_ch3oh / dF_diff_h2o
        q_diff = sim_res[self.comp_num * 5 + 6, -1]
        q_heater = sim_res[self.comp_num * 5 + 7, -1]
        dH_diff = sim_res[self.comp_num * 5 + 8, -1]
        yield_CH3OH = dF_react_ch3oh  # mol/s
        yield_CH3OH_C = self.cond_record['S17']['Methanol'] if self.recycle == 1 else yield_CH3OH  # sim_res[3, -1]  #

        # volume and catalyst utilization
        if self.status == 0:
            V = self.reactor_para.iloc[0]['L'] * self.reactor_para.iloc[0]['Dt'] ** 2 / 4 * np.pi
        else:
            V = self.reactor_para.iloc[0]['L'] * (self.reactor_para.iloc[0]['Dt'] ** 2 -
                                                  self.membrane.iloc[0]['Din'] ** 2) / 4 * np.pi
        mcat = self.reactor_para.iloc[0]['rhoc'] * self.reactor_para.iloc[0]['phi'] * V * 1000  # g
        GHSV = np.sum(sim_res[1:self.comp_num + 1, 0]) * 3600 * 273.15 * 8.314 / 1e5 / V
        FHSV = self.Ft0 * 3600 * 273.15 * 8.314 / 1e5 / V
        STY = yield_CH3OH_C * 3600 * 32.042 * 1000 / mcat  # mg_CH3OH/h/g_cat

        # , sp_ch3oh, sp_h2o, N_CH3OH_H2O
        res = pd.Series([r, r_c, r_H2, s_react, yield_CH3OH, yield_CH3OH_C, dP, To_r,
                         dH_react, q_diff, q_heater, dH_diff, sp_h2o, sp_h2, dF_diff_h2o, dF_diff_h2,
                         GHSV, STY, FHSV],
                        index=['conversion', "r_c", "r_H2", 'select', 'y_CH3OH', 'y_CH3OH_C', 'dP', 'To_r',
                               'q_react', 'q_diff', 'q_heater', "H_diff", "sp_H2O",
                               'sp_H2', 'N_sp_H2O', 'N_sp_H2', 'GHSV', 'STY', "FHSV"])
        print(res)
        return res

    def recycle_metric(self, F_recycle):
        """
        calculate the metric for the whole plant with recycled stream
        :param F_recycle: recycled gas
        :param product: condensed product
        :return:
        """
        # y= z [F_CO2, F_H2, F_CH3OH, F_H2O, F_CO, F_Ar
        # molar fraction
        # react: F_CO2, F_H2, F_CH3OH, F_H2O, F_CO, F_Ar
        # diff: F_CO2, F_H2, F_CH3OH, F_H2O, F_CO, F_Ar
        # permeate: F_CO2, F_H2, F_CH3OH, F_H2O, F_CO, F_Ar
        # Tr, Tc, Th, Pp, q_react, q_diff, q_heater, h_diff]

        # [L1, F_CO2, F_H2, F_CH3OH, F_H2O, F_CO, Tr, Tc, q_react, q_diff]
        # calculate metric for recycled reactor, recycled ratio and enthalpy
        ratio = F_recycle[self.comp_list].sum() / self.Ft0

        p_conversion = self.cond_record["S17"]['Methanol'] / self.F0[0]
        p_conversion_H2 = self.cond_record["S17"]['Methanol'] * 2 / self.F0[1]

        cal = VLEThermo(self.comp_list, ref='ev')

        # calculate the process exergy efficiency
        CO2_feed, H2_feed = self.cond_record['SCO2'], self.cond_record['SH2']
        E_feed = cal.cal_E(CO2_feed['T'], CO2_feed['P'], CO2_feed[self.comp_list].to_numpy()) + \
                 cal.cal_E(H2_feed['T'], H2_feed['P'], H2_feed[self.comp_list].to_numpy())  # kW
        E_per = 0  # exergy of sweep gas
        if self.status == 0:
            E_per_h = 0  # heat exergy of permeate coming from reactor
            q_per = 0
        else:
            try:
                if self.cond_record['S19']['T'] == self.cond_record['S20p']['T']:
                    E_per_h = self.utility_record['Q3']['E']
                    q_per = self.utility_record['Q3']['Q']
                else:
                    E_per_h = 0
                    q_per = 0
            except KeyError:
                E_per_h = self.utility_record['Q3']['E']
                q_per = self.utility_record['Q3']['Q']

            if abs(self.Fp[1]) > 1e-10:  # permeate in is H2
                H2_feed = self.cond_record['S18pH2']
                E_per += cal.cal_E(H2_feed['T'], H2_feed['P'], H2_feed[:-2])
            if abs(self.Fp[0]) > 1e-10:  # permeate in is CO2
                CO2_feed = self.cond_record['S18pCO2']
                E_per += cal.cal_E(CO2_feed['T'], CO2_feed['P'], CO2_feed[:-2])
        try:
            product = self.cond_record['S17']
        except KeyError:
            product = self.cond_record['S9']
        E_product = cal.cal_E(product['T'], product['P'], product[self.comp_list].to_numpy())
        print(self.utility_record)

        E_ht, Q_ht = 0, 0
        if self.membrane.iloc[0]['qmh'] > 0:
            if self.cond_record['CM1p']['T'] == self.cond_record['CM3']['T']:
                E_ht = self.utility_record['CMQ']['E']
                Q_ht = self.utility_record['CMQ']['Q']
        else:
            if self.cond_record['S3pr']['T'] == self.cond_record['S3p']['T']:
                E_ht = self.utility_record['Q1']['E']
                Q_ht = self.utility_record['Q1']['Q']
            else:
                E_ht = 0
                Q_ht = 0

        # if ht stream is not exploited and the heat duty in reboiler not equal to zero
        # then part of duty in reboiler will be compensated by ht stream
        utility_record_mod = self.utility_record.copy()
        if self.utility_record.loc['Q', 'RB'] != 0 and Q_ht != 0:
            utility_record_mod.loc['E', 'RB'] = (utility_record_mod.loc['Q', 'RB'] + Q_ht) / \
                                                utility_record_mod.loc['Q', 'RB'] * utility_record_mod.loc['E', 'RB']
            utility_record_mod.loc['Q', 'RB'] = utility_record_mod.loc['Q', 'RB'] + Q_ht
            utility_record_mod.loc['Q', 'RB'] = max(utility_record_mod.loc['Q', 'RB'], 0)
            utility_record_mod.loc['E', 'RB'] = max(utility_record_mod.loc['E', 'RB'], 0)
            E_ht, Q_ht = 0, 0

        e_duty_term = utility_record_mod.loc[:, utility_record_mod.loc['Q'] > 0]
        q_exergy_duty = e_duty_term.loc['E'].sum()  # self.utility_record.loc['E', 'RB']
        heat_duty = 0
        for col in e_duty_term.columns:
            if e_duty_term[col]['E'] > 0:
                heat_duty += e_duty_term[col]['Q']

        work_duty = np.sum(self.utility_record.loc['W'])  # energy["W"]

        eta = E_product / (E_feed + E_per + q_exergy_duty + work_duty)
        loss_streams = ['S8', 'S14', 'S16']
        exergy_loss = 0
        for s in loss_streams:
            try:
                exergy_loss += cal.cal_E(self.cond_record[s]['T'], self.cond_record[s]['P'],
                                         self.cond_record[s][self.comp_list].to_numpy())
            except KeyError:
                exergy_loss += 0
        eta_mod = (E_product + exergy_loss) / (E_feed + E_per + q_exergy_duty + work_duty)
        eta_mod_mem = (E_product + exergy_loss + abs(E_ht) + abs(E_per_h)) \
                      / (E_feed + E_per + q_exergy_duty + work_duty)
        eta_mod_heat = (E_product + abs(E_ht) + abs(E_per_h)) / (E_feed + E_per + q_exergy_duty + work_duty)
        p_metric = pd.Series([ratio, heat_duty, q_exergy_duty, work_duty, exergy_loss,
                              q_per, E_per_h, E_feed, E_product, eta,
                              eta_mod, eta_mod_mem, p_conversion, p_conversion_H2, E_ht, Q_ht, eta_mod_heat],
                             index=['ratio', 'heat', 'q_exergy', "work", 'e_loss',
                                    'q_per', 'E_per_h', "Ef", "Ep", "eff_exergy",
                                    'eff_e_mod', "eff_e_mod_mem", "p_conv", 'p_conv_H2',
                                    'E_ht', 'Q_ht', 'eta_mod_heat'])
        return p_metric

    def save_point(self):
        Tr = self.feed0_para['T0']
        Pr = self.feed0_para['P0']
        L = self.L[0]
        per_H2O, per_S = self.membrane.iloc[0][['per_H2O', 'per_S']]
        Fp_coe = self.membranes["Fp"]

        mode = 'MR' if self.status == 1 else 'TR'
        point_path = f'result/sim_point_{Tr}_{Pr}_{mode}_{L}_{datetime.now().date()}_{per_H2O}_' \
                     f'{per_S}_{Fp_coe[0]}_{Fp_coe[1]}.xlsx'

        print(self.cond_record)

        # Ensure the directory exists
        os.makedirs(os.path.dirname(point_path), exist_ok=True)

        # Save data to the Excel file
        if os.path.exists(point_path):
            # Append to existing file
            with pd.ExcelWriter(point_path, engine='openpyxl', mode='a', if_sheet_exists='new') as writer:
                self.cond_record.to_excel(writer, index=True, sheet_name='stream')
                self.utility_record.to_excel(writer, index=True, sheet_name='utility')
                if self.rf_record is not None:
                    self.rf_record.to_excel(writer, index=True, sheet_name='rf')
                if self.iso_record is not None:
                    self.iso_record.to_excel(writer, index=True, sheet_name='iso')
        else:
            # Create a new file
            with pd.ExcelWriter(point_path, engine='openpyxl') as writer:
                self.cond_record.to_excel(writer, index=True, sheet_name='stream')
                self.utility_record.to_excel(writer, index=True, sheet_name='utility')
                if self.rf_record is not None:
                    self.rf_record.to_excel(writer, index=True, sheet_name='rf')
                if self.iso_record is not None:
                    self.iso_record.to_excel(writer, index=True, sheet_name='iso')

    @staticmethod
    def drop_tiny(F):
        """
        drop the tiny value in component flow rate
        """
        processed_F = F.copy()
        processed_F = np.where(np.abs(processed_F) < 1e-15, 0, processed_F)
        return processed_F

    def guess_initial(self, F, status, Fp):
        # recycle gas is the mix of unreacted comps and recycled sweep gas
        F_re0 = np.zeros(len(F))
        Ft0 = np.sum(F)
        r_guess = 4 if status == 1 else 7  # recycle ratio
        Fr_coe = np.array([0.19, 0.75, 0.005, 0.001]) if r_guess > 5 else [0.16, 0.78, 0.006, 0.0005]
        if status == 0:
            if self.membrane.iloc[0]['qmh'] != 0:
                if self.P0 == 30:
                    r_coe = [-1.19536776e-04, 1.84139527e-01, -9.44823877e+01, 1.61555924e+04]
                elif self.P0 == 40:
                    r_coe = [-9.52258888e-05, 1.48041635e-01, -7.66870353e+01, 1.32424676e+04]
                elif self.P0 == 50:
                    r_coe = [20.88822191, 14.05505311, 9.332621707, 6.828864878, 5.596684168, 5.100384427, 5.136366892]
                elif self.P0 == 70:
                    r_coe = [20.88822191, 14.05505311, 9.332621707, 6.828864878, 5.596684168, 5.100384427, 5.136366892]
                    Fr_coe = [0.1575, 0.8037, 0.0030, 0.0006]
                r_guess = 6  # np.polyval(r_coe, self.T0)*70/50

            else:
                if self.P0 == 30:
                    r_coe = [-1.15516438e-04, 1.78173857e-01, -9.15314384e+01, 1.56689929e+04]
                elif self.P0 == 40:
                    r_coe = [-3.23078612e-04, 4.81851841e-01, -2.39371700e+02, 3.96154707e+04]
                elif self.P0 == 50:
                    r_coe = [-3.59323727e-04, 5.36645675e-01, -2.66976214e+02, 4.42487746e+04]

                r_guess = 5  # np.polyval(r_coe, self.T0) * 40 / self.P0
        if self.status == 1:
            # r_guess = 6
            Ft0 = Ft0 - np.sum(self.Fp)
            r_coe = [-1.16392203e-06, 2.15661721e-03, -1.23391465e+00, 2.28532416e+02]
            if self.membrane.iloc[0]['qmh'] != 0:
                r_coe = [-1.71450000e-04, 2.62021115e-01, - 1.33457811e+02, 2.26585619e+04]
                if self.P0 == 30:
                    # r_guess = 28
                    r_coe = [-8.12703174e-05, 1.28701438e-01, -6.78187614e+01, 1.18984140e+04]
                    Fr_coe = [0.2188, 0.7433, 0.0065, 0.0012] if self.T0 < 483.15 else [0.1784, 0.7564, 0.0068, 0.0008]
                elif self.P0 == 40:
                    # r_guess = 6
                    Fr_coe = [0.2043, 0.7587, 0.0050, 0.0009] if self.T0 < 483.15 else [0.1503, 0.7921, 0.0058, 0.0005]
                    r_coe = [-6.26421267e-05, 1.01198217e-01, -5.43496627e+01, 9.70960209e+03]
                elif self.P0 == 50:
                    # r_guess = 20
                    r_coe = [-6.35548022e-05, 1.02890353e-01, -5.53772305e+01, 9.91363164e+03]
                    Fr_coe = [0.1935, 0.7700, 0.0043, 0.0007] if self.T0 < 483.15 else [0.1155, 0.8366, 0.0051, 0.0003]
                r_guess = np.polyval(r_coe, self.T0)  # * 40 / self.P0
                r_guess = min(28, r_guess)
            else:
                if self.P0 == 30:
                    # r_guess = 28
                    r_coe = [-1.73661246e-04, 2.61600429e-01, -1.31178015e+02, 2.19047731e+04]
                    Fr_coe = [0.1325, 0.7851, 0.0049, 0.0004]
                elif self.P0 == 40:
                    # r_guess = 6
                    Fr_coe = [0.1400, 0.7432, 0.0055, 0.0006]
                    r_coe = [-6.77568342e-06, 1.07767351e-02, -5.62681028e+00, 9.72488924e+02]
                elif self.P0 == 50:
                    # r_guess = 20
                    r_coe = [-1.16392203e-06, 2.15661721e-03, -1.23391465e+00, 2.28532416e+02]
                    # Fr_coe = [0.1213, 0.7615, 0.0050, 0.0003]
                    Fr_coe = [0.0401, 0.9113, 0.0042, 0.0001]
                r_guess = np.polyval(r_coe, self.T0)  # * 40 / self.P0
            # r_guess = 4
        r_conversion = 0.5 if status == 1 else 0.12  # reactor conversion
        r_guess = 3
        Fr_coe = np.array([0.05, 1, 0.001, 5e-5])
        # Fr_coe = [0.1888,0.7440,0.0061,0.0010]#[0.153, 0.788, 0.006, 0.001]
        for i in range(len(Fr_coe)):
            F_re0[i] = Ft0 * r_guess * Fr_coe[i]
        # print(r_guess)
        F_re0[4] = Ft0 * r_guess * 0.05 if status == 1 else F[0] * r_guess * r_conversion * 1.5
        F_re0[5] = self.Fp[-1] if status == 1 else 0

        return F_re0

    def one_pass(self, reactor, membrane, feed, diff_in=None):
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

        # 提取参数并转换为 PyTorch 张量
        [T_in, P_in, F_in_np] = feed.loc[['T0', 'P0', 'F0']].values
        F_in = torch.tensor(F_in_np, dtype=torch.float32, device=device)

        [L, Dt, Dc, nrt, phi, rhoc] = reactor.loc[['L', 'Dt', "Dc", 'nrt', 'phi', 'rhoc']].values
        [Din, thick, Tc_in, qmc, Th_in, qmh, Pp] = membrane.loc[['Din', 'Thick', 'Tc', 'qmc', 'Th', 'qmh', 'Pp']].values
        [status, nit, per, Fp_in_np] = membrane.loc[['status', 'nit', 'per', 'Fp']].values
        [per_H2O, per_S] = membrane.loc[['per_H2O', 'per_S']]
        per = torch.tensor([0, per_H2O / per_S, 0, per_H2O, 0, 0], device=device)

        # 转换 Fp_in 为 PyTorch 张量
        Fp_in_tensor = torch.tensor(Fp_in_np, dtype=torch.float32, device=device)

        try:
            Tc_in = self.feed_para['T0_H2'] if Fp_in_np[0] < 1E-10 else self.feed_para['T0_CO2']
        except KeyError:
            Tc_in = T_in

        Th_in = T_in
        if self.location == 0:
            Dt = Din
        else:
            Dc = Din + thick * 2
        Dt = Din if self.location == 0 else Dt

        Do = Din + thick * 2
        if diff_in is None:
            q_react_in, q_diff_in, q_heater_in, h_diff_in = 0, 0, 0, 0
        else:
            [q_react_in, q_diff_in, q_heater_in, h_diff_in, F_diff_in] = diff_in

        comp = ["CO2", "H2", "Methanol", "H2O", "CO", 'N2'] if len(F_in_np) != 5 else ["CO2", "H2", "Methanol", "H2O",
                                                                                       "CO"]

        # 创建 Reaction 对象并移动到设备
        react_sim = Reaction(L, Dt, Dc, nrt, phi, rhoc, self.chem_para, T_in, P_in, F_in, comp, self.eos, qmh).to(
            device)

        # 创建 HMT 对象并移动到设备
        if status == 0:
            cal_HMT = HMT(self.comp_list, phi, ds=2e-3, Dt=Dt).to(device)
        else:
            cal_HMT = HMT(self.comp_list, phi, ds=2e-3, Dt=Dt, Dm=Din).to(device)

        # 预计算 PropsSI 的值
        # 假设 cph 是恒定的，或者根据您的需求预计算
        # 如果 cph 是需要随 Tr, P 变化的，您需要找到其他方法处理
        cph = PropsSI('CPMOLAR', 'T', 523, 'P', 48E5, 'water')  # 使用 numpy 计算
        cph_tensor = torch.tensor(cph, dtype=torch.float32, device=device)

        # calculate the HMT performance
        if status == 0:
            cal_HMT = HMT(self.comp_list, phi, ds=2e-3, Dt=Dt)
        else:
            cal_HMT = HMT(self.comp_list, phi, ds=2e-3, Dt=Dt, Dm=Din)
        # 定义 ODE 函数
        def model(y, z):
            F = y[:self.comp_num]
            Fp = y[self.comp_num * 3:self.comp_num * 4]
            Tr, Tc = y[self.comp_num * 4], y[self.comp_num * 4 + 1]
            Th = y[self.comp_num * 4 + 2]
            P = y[self.comp_num * 4 + 3]

            F = torch.where(F > 0, F, torch.zeros_like(F))
            Fp = torch.where(Fp > 0, Fp, torch.zeros_like(Fp))

            res_react = react_sim(Tr, P, F)

            U_mr = 2.4  # W/m2 K
            U_rh = cal_HMT.htc4(Tr, P, F) if qmh > 0 else torch.tensor(0.0, device=device)
            dl2dw = torch.pi * ((Dt ** 2) / 4) * rhoc * phi
            dP_dz = torch.tensor(0.0, device=device) if self.drop == 0 else cal_HMT.ergun(Tr, P, F) * 1e-5

            if status == 1:
                r_v_ins_v_react = Do ** 2 / Dt ** 2 if self.location == 1 else 0
                Ap = torch.pi * Do
                Pi_r = F / torch.sum(F) * P
                Pi_p = Fp / torch.sum(Fp) * Pp
                Pi_dp = Pi_r - Pi_p
                mflux = Pi_dp * 1e5 * per * Ap  # mol/m s
                mflux_r2p = torch.where(mflux > 0, mflux, torch.zeros_like(mflux))
                mflux_p2r = torch.where(mflux < 0, -mflux, torch.zeros_like(mflux))

                H_per_delta_r2p = self.pro_cal.cal_H(Tr, P, mflux_r2p) - self.pro_cal.cal_H(Tc, Pp, mflux_r2p)
                H_per_delta_p2r = torch.where(torch.sum(mflux_p2r) < 1e-10,
                                              torch.tensor(0.0, device=device),
                                              self.pro_cal.cal_H(Tc, Pp, mflux_p2r) - self.pro_cal.cal_H(Tr, P,
                                                                                                         mflux_p2r))
                H_per_delta = H_per_delta_r2p + H_per_delta_p2r

                qc = U_mr * (Tr - Tc) * torch.pi * Do
            else:
                r_v_ins_v_react = 0
                mflux = torch.zeros_like(F)
                H_per_delta = torch.tensor(0.0, device=device)
                qc = torch.tensor(0.0, device=device)

            dF_react_dz = res_react['mflux'] * dl2dw * (1 - r_v_ins_v_react * nit) * self.nrt
            dF_diff_dz = -mflux * nit
            dF_dz = dF_react_dz + dF_diff_dz
            dFp_dz = mflux

            cpr = self.pro_cal.cal_cp(Tr, P, F)
            heat_cap_r = cpr * torch.sum(F)

            cpp = self.pro_cal.cal_cp(Tc, P, Fp)
            heat_cap_p = cpp * torch.sum(Fp)
            # 使用预计算的 cph 值
            cph_current = cph_tensor

            q_heater = U_rh * (Tr - Th) * torch.pi * Dt if qmh != 0 else torch.tensor(0.0, device=device)
            q_react = res_react["hflux"] * dl2dw * (1 - r_v_ins_v_react) * nrt

            dTr_dz = (q_react - qc * nit - q_heater - H_per_delta * nit) / heat_cap_r
            dTc_dz = (qc * nit + H_per_delta * nit) / heat_cap_p
            dTh_dz = q_heater / qmh / cph_current if qmh != 0 else torch.tensor(0.0, device=device)
            dq_rea_dz = res_react["hflux"] * dl2dw * (1 - r_v_ins_v_react) * nrt
            dq_dif_dz = -q_heater - qc * nit
            dh_dif_dz = H_per_delta * nit

            res_dz = torch.cat((
                dF_dz,
                dF_react_dz,
                dF_diff_dz,
                dFp_dz,
                torch.tensor([dTr_dz, dTc_dz, dTh_dz, dP_dz, dq_rea_dz, dq_dif_dz, q_heater, dh_dif_dz],
                             device=device)
            ))

            return res_dz

        # JIT 编译模型函数
        model_jit = torch.jit.script(model)

        # 定义积分点
        z_span = torch.linspace(0, L, 1000, device=device)
        Tc_ini = Tc_in  # if status == 1 else self.T_feed
        Th_ini = Th_in
        # 初始条件
        ic = torch.cat((
            F_in,
            torch.zeros(self.comp_num, device=device),
            torch.zeros(self.comp_num, device=device),
            Fp_in_tensor,
            torch.tensor([T_in, Tc_ini, Th_ini, P_in, q_react_in, q_diff_in, q_heater_in, h_diff_in], device=device)
        ))

        # 使用 torchdiffeq 的 odeint 进行 ODE 求解
        res_sim = odeint(model_jit, ic, z_span)

        # 将结果转换回 CPU 和 NumPy 以便后续处理
        sim_profile = res_sim.cpu().detach().numpy().T

        # 处理微小值
        sim_profile[0:self.comp_num, :] = self.drop_tiny(sim_profile[0:self.comp_num, :])
        sim_profile[self.comp_num * 3:self.comp_num * 4, :] = self.drop_tiny(
            sim_profile[self.comp_num * 3:self.comp_num * 4, :])

        res = np.vstack((z_span.cpu().numpy(), sim_profile))
        if qmh != 0:
            # cal_HMT.htc4 需要支持 NumPy 数组或预计算
            U_rh_in = cal_HMT.htc4(T_in, P_in, F_in.cpu().numpy())
            T_o, P_o = res[self.comp_num * 4 + 1, -1], res[self.comp_num * 4 + 4, -1]
            Fo = res[1:self.comp_num + 1, -1]
            U_rh_o = cal_HMT.htc4(T_o, P_o, Fo)
            self.iso_record = pd.Series([U_rh_in, U_rh_o], index=['Uin', 'Uo'])
        return res

    def multi_reactor(self):
        res = {}
        feed_para, diff_para = self.feed0_para, None

        for n in range(self.stage):
            mem_para = self.membrane.iloc[n]
            reactor_para = self.reactor_para.iloc[n]
            res[f'{n}'] = self.one_pass(reactor_para, mem_para, feed_para, diff_para)
            F_out = res[f'{n}'][1:(self.comp_num + 1), -1].copy()
            r_temp = (feed_para['F0'][0] - F_out[0]) / feed_para['F0'][0]

            Tr_out = res[f'{n}'][self.comp_num * 4 + 1, -1].copy()

            P_out = res[f'{n}'][self.comp_num * 4 + 4, -1].copy()
            F_diff = res[f'{n}'][self.comp_num * 2:self.comp_num * 3, -1].copy()
            diff_para = res[f'{n}'][self.comp_num * 4 + 5:self.comp_num * 4 + 9, -1].copy().tolist()
            diff_para.append(F_diff)
            feed_para = pd.Series([Tr_out, P_out, F_out], index=['T0', 'P0', 'F0'])

            print(Tr_out, r_temp)
        if self.stage >= 1:
            res_profile = res['0']

            custom_array = np.linspace(0, self.L[0], 1000).tolist()
            for i in np.arange(1, self.stage):
                res_profile = np.hstack((res_profile, res[f'{i}']))
                custom_array += list(np.linspace(0, self.L[i], 1000) + sum(self.L[:i]))
            res_profile[0] = custom_array
        xi = res_profile[1:self.comp_num + 1] / np.sum(res_profile[1:self.comp_num + 1], axis=0)
        res_profile = np.insert(res_profile, self.comp_num + 1, xi, axis=0)
        return res_profile

    def adjust_feed(self):
        """
        adjust feed comp if membrane is used
        :return:
        """
        self.feed0_para['F0'] = self.feed0_para['F0'] - self.Fp

    def process_feed(self, Tr, Pr):

        # construct the feed info
        F_feed = self.feed0_para['F0']
        feed_pd = pd.Series(F_feed, index=self.comp_list)

        T_CO2, T_H2 = self.feed_para['T0_CO2'], self.feed_para['T0_H2']
        P_CO2, P_H2 = self.feed_para['P0_CO2'], self.feed_para['P0_H2']
        CO2_pd = pd.concat([feed_pd, pd.Series([T_CO2, P_CO2], index=['T', 'P'])])
        CO2_pd["H2"] = 0
        H2_pd = pd.concat([feed_pd, pd.Series([T_H2, P_H2], index=['T', 'P'])])
        H2_pd["CO2"] = 0
        self.cond_record['SCO2'] = CO2_pd.copy()
        self.cond_record['SH2'] = H2_pd.copy()

        # compress them
        CO2_comp_res = multi_comp_opt(CO2_pd, P2=Pr)
        if Pr > H2_pd['P']:
            H2_comp_res = multi_comp_opt(H2_pd, P2=Pr)
            H2_comp = H2_comp_res.iloc[-1, :-2]
            W_comp_H2 = H2_comp_res['W'].sum()
            Q_comp_H2 = H2_comp_res['Q'].sum()
            self.utility_record['C2'] = [Q_comp_H2, W_comp_H2, W_comp_H2]
        else:
            H2_comp = H2_pd.copy()
        CO2_comp = CO2_comp_res.iloc[-1, :-2]
        W_comp_CO2 = CO2_comp_res['W'].sum()
        Q_comp_CO2 = CO2_comp_res['Q'].sum()

        self.utility_record['C1'] = [Q_comp_CO2, W_comp_CO2, W_comp_CO2]

        # mix them
        feed_mix_res = mixer(CO2_comp, H2_comp)
        # preheat them
        # feed_heat_res = heater(feed_mix_res['Fo'], T_out=453.15)
        # self.utility_record['Q1p'] = [feed_heat_res['Q'], feed_heat_res['E'], 0]
        return feed_mix_res  # feed_heat_res

    def valve_permeate(self, feed, p):
        res = valve(feed, P=p)
        return res['Fo']

    def post_recovery_ad(self):
        """
        heat recovery for adiabatic loop
        """
        # post process
        if self.cond_record['S3p']['T'] > self.cond_record['S3']['T']:
            try:
                if self.status == 0:
                    self.post_process(feed=self.cond_record['S9'], lt_stream=self.cond_record['S11'],
                                      ht_stream=self.cond_record['S3p'])
                else:
                    self.post_process(feed=self.cond_record['S9'], lt_stream=self.cond_record['S11'],
                                      ht_stream=self.cond_record['S3p'], permeate=self.cond_record['S19'])
            except ValueError:
                self.post_process(feed=self.cond_record['S9'], lt_stream=self.cond_record['S11'],
                                  ht_stream=None)
        else:
            self.post_process(feed=self.cond_record['S9'], lt_stream=self.cond_record['S11'],
                              ht_stream=None)
        # save status

    def post_recovery_iso(self, profile):
        """
        heat recovery for isothermal loop
        """

        Th_out = profile[self.comp_num * 4 + 3, -1]
        # Th_out_pd = pd.Series([Th_out, 48, self.membrane.iloc[0]['qmh']], index=['T', 'P', 'qm'])
        Th_out_record = pd.Series(np.zeros(self.comp_num + 2), index=self.comp_list + ['T', 'P'])
        Th_out_record['H2O'] = self.membrane.iloc[0]['qmh']
        Th_out_record['T'] = Th_out
        Th_out_record['P'] = 48
        F_fresh, T_r = self.feed0_para['F0'], self.feed0_para['T0']

        if self.cond_record['S3p']['T'] > self.cond_record['S3']['T']:
            print('reactor inlet temperature higher than prescribed value')
            exit(0)
            self.post_process(feed=self.cond_record['S9'], lt_stream=self.cond_record['S11'],
                              ht_stream=self.cond_record['S3p'])
        else:
            # heat the reactor feed using water
            hx_ht = HeatExchanger(Th_out_record, self.cond_record['S3p'], U=500)
            hx_ht_res = hx_ht.fixed_T_cold(Tout=T_r)
            hx_ht_ht_out = hx_ht_res['Fh_o']
            self.cond_record['CM2'] = Th_out_record.copy()
            self.cond_record['CM3'] = hx_ht_ht_out.copy()
            if self.cond_record['CM3']['T'] >= self.T0:
                self.utility_record.loc['E', 'Q1'] = 0
            else:
                hx_ht_res = hx_ht.fixed_T_hot(Tout=T_r)
                hx_ht_ht_out = hx_ht_res['Fh_o']
                self.cond_record['CM3'] = hx_ht_ht_out.copy()
                self.utility_record['CMQ'] = [0, 0, 0]
                # introduce external source to heat further reactor feed
                feed_pre_heated = hx_ht_res['Fc_o']
                self.cond_record['S3ph'] = feed_pre_heated.copy()
                feed_fur_heater = heater(fluid=feed_pre_heated, T_out=self.T0)
                self.utility_record.loc['E', 'Q1'] = 0
                self.utility_record.loc['Q', 'Q1h'] = feed_fur_heater['Q']
                self.utility_record.loc['E', 'Q1h'] = feed_fur_heater['E']
                # print('fail in heating the reactor feed')
                # exit(0)
            # the extra heat in water is fed to reboiler
            if self.status == 0:
                self.post_process(feed=self.cond_record['S9'], lt_stream=self.cond_record['S11'],
                                  ht_stream=hx_ht_ht_out.copy())
            else:
                if hx_ht_ht_out['T'] > self.T0:
                    self.post_process(feed=self.cond_record['S9'], lt_stream=self.cond_record['S11'],
                                      ht_stream=hx_ht_ht_out.copy(), permeate=self.cond_record['S19'])
                else:
                    self.post_process(feed=self.cond_record['S9'], lt_stream=self.cond_record['S11'],
                                      ht_stream=None, permeate=self.cond_record['S19'])

    def pre_process_permeate(self, permeate):

        permeate_feed = pd.Series(self.Fp, index=self.comp_list)
        Pp = self.membrane.iloc[0]['Pp']
        self.utility_record['Cper'] = [0, 0, 0]
        if abs(self.Fp[1]) > 1e-10:  # permeate in is H2
            permeate_H2_cond = pd.Series([self.feed_para['T0_H2'], self.feed_para['P0_H2']], index=['T', 'P'])
            permeate_H2_raw = pd.concat([permeate_feed, permeate_H2_cond])
            permeate_H2_raw['CO2'] = 0
            if permeate_H2_raw['P'] > Pp:
                permeate_H2_process = self.valve_permeate(permeate_H2_raw, Pp)
            else:
                permeate_H2_process = permeate_H2_raw
        else:
            permeate_H2_raw = None
            permeate_H2_process = None

        if abs(self.Fp[0]) > 1e-10:  # permeate in is CO2
            permeate_CO2_cond = pd.Series([self.feed_para['T0_CO2'], self.feed_para['P0_CO2']], index=['T', 'P'])
            permeate_CO2_raw = pd.concat([permeate_feed, permeate_CO2_cond])
            permeate_CO2_raw['H2'] = 0
            if permeate_CO2_raw['P'] < Pp:
                permeate_CO2_process_res = multi_comp_opt(permeate_CO2_raw, Pp)
                permeate_CO2_process = permeate_CO2_process_res.iloc[-1, :-2]
                W_comp = permeate_CO2_process_res['W'].sum()
                Q_comp = permeate_CO2_process_res['Q'].sum()
                self.utility_record['Cper'] = [W_comp, Q_comp, 0]
            else:
                permeate_CO2_process = permeate_CO2_raw
        else:
            permeate_CO2_raw = None
            permeate_CO2_process = None

        if abs(self.Fp[1]) > 1e-10 and abs(self.Fp[0]) > 1e-10:
            permeate_in_mix = mixer(permeate_H2_process, permeate_CO2_process)
            self.cond_record['S18'] = permeate_in_mix['Fo']
            self.cond_record['S18pCO2'] = permeate_CO2_raw
            self.cond_record['S18pH2'] = permeate_H2_raw
        elif abs(self.Fp[1]) > 1e-10 and abs(self.Fp[0]) < 1e-10:
            self.cond_record['S18'] = permeate_H2_process
            self.cond_record['S18pH2'] = permeate_H2_raw
        elif abs(self.Fp[1]) < 1e-10 and abs(self.Fp[0]) > 1e-10:
            self.cond_record['S18'] = permeate_CO2_process
            self.cond_record['S18pCO2'] = permeate_CO2_raw

        self.cond_record['S19'] = permeate

    def save_ht(self, sp_res, permeate_heater):
        if self.membrane.iloc[0]['qmh'] > 0:
            self.cond_record['CM1p'] = sp_res['ht']
            if sp_res['ht']['T'] > self.T0:
                cooler = heater(sp_res['ht'].copy(), T_out=self.T0)
                self.cond_record['CM1'] = cooler['Fo']
                self.utility_record['CMQ'] = [cooler['Q'], cooler['E'], 0]
            else:
                self.cond_record['CM1'] = self.cond_record['CM1p'].copy()
        else:
            self.cond_record['S3pr'] = sp_res['ht']
            if sp_res['ht']['T'] > self.T0:
                cooler = heater(sp_res['ht'].copy(), T_out=self.T0)
                self.cond_record['S3'] = cooler['Fo']
                self.utility_record['Q1'] = [cooler['Q'], cooler['E'], 0]
            else:
                self.cond_record['S3'] = sp_res['ht'].copy()
                self.utility_record['Q1'] = [0, 0, 0]
        if permeate_heater is not None:
            self.cond_record['S20p'] = permeate_heater['Fh_o']
            self.utility_record['Q3p'] = [permeate_heater['Q'], 0, 0]
            self.cond_record['S15pa'] = permeate_heater['Fc_o']
            cooler = heater(permeate_heater['Fh_o'], T_out=308.15)
            self.cond_record['S20'] = cooler['Fo']
            self.utility_record['Q3'] = [cooler['Q'], cooler['E'], 0]

    def save_ht1(self, sp_res, permeate):
        if not permeate:
            if self.membrane.iloc[0]['qmh'] > 0:
                self.cond_record['CM1p'] = sp_res['ht']
                if sp_res['ht']['T'] > self.T0:
                    cooler = heater(sp_res['ht'].copy(), T_out=self.T0)
                    self.cond_record['CM1'] = cooler['Fo']
                    self.utility_record['CMQ'] = [cooler['Q'], cooler['E'], 0]
                else:
                    self.cond_record['CM1'] = self.cond_record['CM1p'].copy()
            else:
                self.cond_record['S3pr'] = sp_res['ht']
                if sp_res['ht']['T'] > self.T0:
                    cooler = heater(sp_res['ht'].copy(), T_out=self.T0)
                    self.cond_record['S3'] = cooler['Fo']
                    self.utility_record['Q1'] = [cooler['Q'], cooler['E'], 0]
                else:
                    self.cond_record['S3'] = sp_res['ht'].copy()
                    self.utility_record['Q1'] = [0, 0, 0]
        else:
            self.cond_record['S20p'] = sp_res['ht']
            if sp_res['ht']['T'] > 308.15:
                cooler = heater(sp_res['ht'].copy(), T_out=308.15)
                self.cond_record['S20'] = cooler['Fo']
                self.utility_record['Q3'] = [cooler['Q'], cooler['E'], 0]
            else:
                self.cond_record['S20'] = sp_res['ht'].copy()
                self.utility_record['Q3'] = [0, 0, 0]

    def save_ht_no_rec(self, sp_heater=None):
        if self.membrane.iloc[0]['qmh'] > 0:
            self.cond_record['CM1p'] = self.cond_record['CM3'].copy()
            if abs(self.cond_record['CM1p']['T'] - self.T0) > 0.01:
                cooler = heater(self.cond_record['CM1p'], T_out=self.T0)
                self.cond_record['CM1'] = cooler['Fo']
                self.utility_record['CMQ'] = [cooler['Q'], cooler['E'], 0]
            else:
                self.cond_record['CM1'] = self.cond_record['CM3'].copy()
                self.cond_record['CM1']['T'] = self.T0
                self.utility_record['CMQ'] = [0, 0, 0]
        else:
            self.cond_record['S3pr'] = self.cond_record['S3p'].copy()
            cooler = heater(self.cond_record['S3p'], T_out=self.T0)
            self.utility_record['Q1'] = [cooler['Q'], cooler['E'], 0]

        if self.status == 1 and sp_heater is None:
            cooler = heater(self.cond_record['S19'], T_out=308.15)
            self.cond_record['S20'] = cooler['Fo']
            self.cond_record['S20p'] = self.cond_record['S19']
            self.utility_record['Q3'] = [cooler['Q'], cooler['E'], 0]
        elif self.status == 1 and sp_heater is not None:
            self.cond_record['S20p'] = sp_heater['Fh_o']
            self.utility_record['Q3p'] = [sp_heater['Q'], 0, 0]
            self.cond_record['S15pa'] = sp_heater['Fc_o']
            cooler = heater(sp_heater['Fh_o'], T_out=308.15)
            self.cond_record['S20'] = cooler['Fo']
            self.utility_record['Q3'] = [cooler['Q'], cooler['E'], 0]

    def post_process_ht(self, feed, ht_stream, Tdew, permeate=None):

        # prepare for distiller
        apk_path = r"E:\project\lsj_y_CH3OH_C\MTC_MEM/DT_with_MR_opt_boiler.bkp"
        block_name = 'B3'
        feed_name = 'S18'
        heavy_name = 'S8'
        light_name = 'S7'

        # calculate the max heat can be recovered in HT stream
        Qr_max = heater(ht_stream, T_out=self.feed0_para['T0'])['Q']

        # calculate the max heat duty in reboiler
        sp = DistillerOpt(apk_path, block_name=block_name, feed_name=feed_name,
                          light_name=light_name, heavy_name=heavy_name)
        sp_in = feed
        sp_res_check = sp.run_RF(stream=sp_in, valid_comp=self.comp_list, ht_stream=None)
        rb_hd = sp_res_check['block']['HD']
        ht_recovery = 0  # ht stream recovery in reboiler failed
        sp_res = None
        sp_in_heater_res = None

        if rb_hd < abs(Qr_max):
            # the waste heat in HT stream can supply the separation
            sp = DistillerOpt(apk_path, block_name=block_name, feed_name=feed_name,
                              light_name=light_name, heavy_name=heavy_name)

            sp_res = sp.run_RF(stream=sp_in, valid_comp=self.comp_list, ht_stream=ht_stream)
            ht_recovery = 1
            # self.save_ht(sp_res, permeate)

        elif abs(Qr_max / rb_hd) > 0.05:
            preheat_max_heater = heater(sp_in, T_out=Tdew)
            preheat_max = preheat_max_heater['Q']  # distiller feed cannot be preheated more than T_pew
            q_diff = rb_hd + Qr_max  # energy gap with HT stream heat recovery
            if permeate is None:
                # the waste heat in HT stream cannot supply the separation directly
                # but there is enough waste heat should be recovered
                print(rb_hd, Qr_max)
                q_diff_coes = np.arange(0.6, 2.2, 0.2)
                for q_diff_coe in q_diff_coes:
                    # print(q_diff_coe)
                    q_duty = min(preheat_max, q_diff_coe * q_diff)  # extra stream to preheat product
                    lt_pre_heater_res = heater_duty(sp_in, duty=q_duty)
                    sp_in_temp = lt_pre_heater_res['Fo']
                    # print(sp_in_temp)
                    if sp_in_temp['T'] > Tdew + 20:
                        print('over-heated product stream')
                        break
                    sp = DistillerOpt(apk_path, block_name=block_name, feed_name=feed_name,
                                      light_name=light_name, heavy_name=heavy_name)
                    sp_res_check2 = sp.run_RF(stream=sp_in_temp, valid_comp=self.comp_list, ht_stream=ht_stream)
                    rb_hd_cal = sp_res_check2['block']['HD']
                    time.sleep(1)
                    if rb_hd_cal < abs(Qr_max):
                        # water heat of water can supply the boiler
                        sp_res = sp_res_check2
                        sp_in = sp_in_temp
                        self.utility_record['QpD'] = [lt_pre_heater_res['Q'], lt_pre_heater_res['E'], 0]
                        # self.save_ht(sp_res, permeate)
                        ht_recovery = 1
                        print('success recovery in RB')
                        break
                    if preheat_max < min(q_diff_coes) * q_diff:
                        # waste heat cannot supply the boiler though the product is heated to Tdew
                        break
                if ht_recovery == 0:
                    sp = DistillerOpt(apk_path, block_name=block_name, feed_name=feed_name,
                                      light_name=light_name, heavy_name=heavy_name)
                    sp_res = sp.run_RF(stream=sp_in, valid_comp=self.comp_list, ht_stream=None)
            else:
                Qpr_max = heater(permeate, T_out=Tdew)['Q']
                if abs(Qpr_max + Qr_max) > rb_hd and permeate['T'] > (feed['T'] + 10):
                    # 渗透流及高温流或可以支撑产物预热及蒸馏塔热荷
                    sp_in_heater = HeatExchanger(permeate, sp_in)

                    delta_T_start = 1
                    for i in range(5):
                        delta_T = delta_T_start + i
                        try:
                            sp_in_heater_res = sp_in_heater.fixed_delta_hc(delta_T=delta_T)
                            if sp_in_heater_res['Fh_o']['T'] > permeate['T'] - 10:
                                # heat recovery fails
                                raise ValueError
                        except ValueError:
                            continue
                    sp_in_temp = sp_in_heater_res['Fc_o']

                    sp = DistillerOpt(apk_path, block_name=block_name, feed_name=feed_name,
                                      light_name=light_name, heavy_name=heavy_name)

                    sp_res_temp = sp.run_RF(stream=sp_in_temp, valid_comp=self.comp_list, ht_stream=ht_stream)

                    rb_hd_cal = sp_res_temp['block']['HD']
                    print('1st')
                    print(f"rb_hd_cal:{rb_hd_cal}")
                    print(f"Qr_max:{Qr_max}")
                    if rb_hd_cal < abs(Qr_max):
                        # ht can supply the boiler
                        sp_res = sp_res_temp
                        sp_in = sp_in_temp
                        # self.save_ht(sp_res, permeate)
                        ht_recovery = 1
                        print('success recovery in RB')
                if ht_recovery == 0 and permeate['T'] > (feed['T'] + 10):
                    # the combination of ht stream and permeate stream cannot supply the product separation
                    print('try to introduce external preheater besides permeate')
                    sp_in_heater = HeatExchanger(permeate, sp_in)
                    delta_T_start = 1
                    for i in range(5):
                        delta_T = delta_T_start + i
                        try:
                            sp_in_heater_res = sp_in_heater.fixed_delta_hc(delta_T=delta_T)
                            if sp_in_heater_res['Fh_o']['T'] > permeate['T'] - 10:
                                # heat recovery fails
                                raise ValueError
                        except ValueError:
                            continue
                    Qpr_max = sp_in_heater_res['Q']
                    q_diff = rb_hd + Qr_max - abs(Qpr_max)
                    q_diff_coes = np.arange(0.6, 2.2, 0.2)
                    print(f'sp_in_heater: {sp_in_heater_res}')
                    print(Qpr_max, Qr_max, rb_hd)
                    # try to introduce an external heat to preheat the product further
                    sp_in_pre_temp = sp_in_heater_res['Fc_o']
                    for q_diff_coe in q_diff_coes:
                        print(f'q_diff_coe:{q_diff_coe}')

                        q_duty = q_diff_coe * q_diff  # extra stream to preheat product
                        lt_pre_heater_res = heater_duty(sp_in_pre_temp, duty=q_duty)
                        sp_in_temp = lt_pre_heater_res['Fo']
                        print(f'sp_in_pre_temp[T]:{lt_pre_heater_res}')
                        if sp_in_temp['T'] > Tdew + 20:
                            break

                        sp = DistillerOpt(apk_path, block_name=block_name, feed_name=feed_name,
                                          light_name=light_name, heavy_name=heavy_name)
                        sp_res_check2 = sp.run_RF(stream=sp_in_temp, valid_comp=self.comp_list, ht_stream=ht_stream)
                        rb_hd_cal = sp_res_check2['block']['HD']
                        time.sleep(1)
                        if rb_hd_cal < abs(Qr_max):
                            # water heat of water can supply the boiler
                            sp_res = sp_res_check2
                            sp_in = sp_in_temp
                            self.utility_record['QpD'] = [lt_pre_heater_res['Q'], lt_pre_heater_res['E'], 0]
                            ht_recovery = 1
                            print('success recovery in RB with additional preheater')
                            break
                if ht_recovery == 0 and permeate['T'] > (feed['T'] + 10):
                    # try to use permeate to preheat the spin then use an external hs to supple reboiler
                    sp_in_heater = HeatExchanger(permeate, sp_in)
                    if sp_in['T'] < Tdew:
                        delta_T_start = 1
                        for i in range(5):
                            delta_T = delta_T_start + i
                            try:
                                sp_in_heater_res = sp_in_heater.fixed_delta_hc(delta_T=delta_T)
                                if sp_in_heater_res['Fh_o']['T'] > permeate['T'] - 10:
                                    # heat recovery fails
                                    raise ValueError
                            except ValueError:
                                continue
                        sp_in = sp_in_heater_res['Fc_o']
                    else:
                        sp_in_heater_res = None
                    print("try to use permeate to preheat the spin then use an external hs to supple reboiler")
                    print(sp_in_heater_res)
                    sp = DistillerOpt(apk_path, block_name=block_name, feed_name=feed_name,
                                      light_name=light_name, heavy_name=heavy_name)
                    sp_res = sp.run_RF(stream=sp_in, valid_comp=self.comp_list, ht_stream=None)
        else:
            if permeate is None:
                sp_in_heater = heater(sp_in, T_out=Tdew)
                sp_in = sp_in_heater['Fo']
                self.utility_record['QpD'] = [sp_in_heater['Q'], sp_in_heater['E'], 0]
            else:
                if permeate['T'] > (feed['T'] + 10):
                    sp_in_heater = HeatExchanger(permeate, sp_in)
                    delta_T_start = 1
                    for i in range(5):
                        delta_T = delta_T_start + i
                        try:
                            sp_in_heater_res = sp_in_heater.fixed_delta_hc(delta_T=delta_T)
                            if sp_in_heater_res['Fh_o']['T'] > permeate['T'] - 10:
                                # heat recovery fails
                                raise ValueError
                        except ValueError:
                            continue
                    sp_in = sp_in_heater_res['Fc_o']
                else:
                    ht_recovery = 0
                    sp_in = feed
                    sp_in_heater_res = None
            sp = DistillerOpt(apk_path, block_name=block_name, feed_name=feed_name,
                              light_name=light_name, heavy_name=heavy_name)
            sp_res = sp.run_RF(stream=sp_in, valid_comp=self.comp_list, ht_stream=None)

        self.cond_record['S15'] = sp_in
        if permeate is None:  # or ht_recovery == 0
            sp_in_heater_res = None
        return ht_recovery, sp_res, sp_in_heater_res

    def post_process_no_rc(self, feed, T_dew):
        # prepare for distiller
        apk_path = r"D:\study\00课题\06多联产系统\甲醇单元\反应器比较\膜反应器\DT_with_MR_opt_boiler.bkp"
        block_name = 'B3'
        feed_name = 'S18'
        heavy_name = 'S8'
        light_name = 'S7'

        # preheat using extra heat source
        if feed['T'] < T_dew - 10:
            pre_heater = heater(feed, T_out=T_dew - 10)
            self.utility_record['QpD'] = [pre_heater['Q'], pre_heater['E'], 0]
            sp_in = pre_heater['Fo']
        else:
            sp_in = feed
        sp = DistillerOpt(apk_path, block_name=block_name, feed_name=feed_name,
                          light_name=light_name, heavy_name=heavy_name)
        sp_res = sp.run_RF(stream=sp_in, valid_comp=self.comp_list, ht_stream=None)
        self.cond_record['S15'] = sp_in
        return sp_res

    def post_process(self, feed, ht_stream, lt_stream, permeate=None):

        # throttle
        valve_res = valve(feed, P=1.2)
        liq_valve = valve_res['Fo']
        flash_res = flasher(liq_valve)
        flash_liq, flash_gas = flash_res['Fl_o'], flash_res['Fg_o']

        # distill
        # heat liquid before entering column using reacted flow
        lt_hx = HeatExchanger(lt_stream, flash_liq, U=500)

        # calculate the dew temperature of product
        try:
            T_dew_in = self.vle_cal.T_dew(flash_liq['P'], flash_liq[self.comp_list].values)
        except ValueError:
            T_dew_in = 365.15

        lt_hx_res = lt_hx.fixed_delta_cold(delta_T=20)
        lt_recovery = lt_hx_res['Fh_o']
        lt_cooler = heater(lt_recovery, T_out=308.15)
        liq_heated_lt = lt_hx_res['Fc_o']
        if liq_heated_lt['T'] >= T_dew_in + 10:
            lt_hx_res = lt_hx.fixed_T_cold(Tout=T_dew_in + 10)
            lt_recovery = lt_hx_res['Fh_o']
            lt_cooler = heater(lt_recovery, T_out=308.15)
            liq_heated_lt = lt_hx_res['Fc_o']

        self.cond_record['S15p'] = liq_heated_lt.copy()
        sp_in = liq_heated_lt

        if ht_stream is not None:
            ht_recover, sp_res, sp_heater = self.post_process_ht(feed=sp_in, ht_stream=ht_stream,
                                                                 Tdew=T_dew_in, permeate=permeate)
            if ht_recover == 0:
                self.save_ht_no_rec(sp_heater=sp_heater)

            else:
                self.save_ht(sp_res, permeate_heater=sp_heater)
        else:
            sp_res = self.post_process_no_rc(sp_in, T_dew_in)
            ht_recover = 0
            self.save_ht_no_rec()

        product = sp_res['light'].copy()
        heavy_product = sp_res['heavy'].copy()
        reboiler_hd = sp_res['block']['HD']
        reboiler_T = sp_res['block']['HT']
        reboiler_hd_e = reboiler_hd * (1 - 298.15 / reboiler_T) * (1 - ht_recover)

        # record
        self.cond_record['S12'] = liq_valve.copy()
        self.cond_record['S13'] = flash_liq
        self.cond_record['S14'] = flash_gas
        self.cond_record['S5p'] = lt_recovery

        self.cond_record['S16'] = heavy_product
        self.cond_record['S17'] = product

        self.utility_record['RB'] = [reboiler_hd, reboiler_hd_e, 0]
        self.utility_record['Q2p'] = [lt_hx_res['Q'], 0, 0]
        self.utility_record['Q2'] = [lt_cooler['Q'], 0, 0]

        try:
            self.rf_record = sp_res['RF']
        except KeyError:
            self.rf_record = None
        return product

    def process_permeate(self, feed, Pr):
        """
        process permeate after reaction
        :param feed: permeate info
        :param Pr: reaction pressure
        :return:
        """
        # Cool and separate permeate stream
        permeate_cool = heater(feed, T_out=308.15)
        flash_res = flasher(permeate_cool['Fo'])
        flash_gas = flash_res['Fg_o']
        filter_res = adsorber(flash_gas)
        gas_permeate = filter_res["Fo"]

        # Further compress and heat permeate
        if Pr > self.membrane.iloc[0]['Pp']:
            compressed_result = multi_comp_opt(gas_permeate, P2=Pr)
            permeate_compressed = compressed_result.iloc[-1, :-2]
            W_comp = compressed_result['W'].sum()
            Q_comp = compressed_result['Q'].sum()
        else:
            permeate_compressed = gas_permeate
            W_comp = 0
            Q_comp = 0

        # record data
        self.cond_record['S20'] = permeate_cool['Fo'].copy()
        self.cond_record['S21'] = flash_gas.copy()
        self.cond_record['S22'] = flash_res['Fl_o'].copy()
        self.cond_record['S23'] = gas_permeate
        self.cond_record['S24'] = permeate_compressed
        self.utility_record['C3'] = [Q_comp, W_comp, W_comp]
        self.utility_record['Q3'] = [permeate_cool['Q'], permeate_cool['E'], 0]

        return permeate_compressed

    def recycle_loop(self, split_ratio, T_cool, max_iter=100, tol=1e-4):

        if self.status == 1:
            self.adjust_feed()

        F_fresh, T_r = self.feed0_para['F0'], self.feed0_para['T0']
        Pin = self.feed0_para['P0']
        # process the feed
        preprocess = self.process_feed(T_r, Pin)
        mix_fresh = preprocess['Fo']

        F_fresh_pd = mix_fresh.copy()
        self.cond_record['S1'] = F_fresh_pd.copy()
        Pp = self.membrane.iloc[0]['Pp']  # self.membranes.loc["Pp1"]
        Fp = self.membranes.loc["Fp"]
        recycle = self.guess_initial(self.Fr0, self.status, Fp)

        T_recycle = T_r
        recycle_pd = pd.concat([pd.Series(recycle, index=self.comp_list),
                                pd.Series([T_recycle, Pin], index=["T", "P"])])

        for i in range(max_iter):

            print(f"Iteration {i}")
            # print(f"Recycle: {recycle_pd}")
            # Mix fresh feed and recycle streams
            recycle_mixer = mixer(F_fresh_pd, recycle_pd)

            # heat the mix if necessary
            # if recycle_mixer['Fo']['T'] > T_r + 100:
            #     reactor_in = recycle_mixer['Fo']
            # else:
            mix_heater = heater(recycle_mixer['Fo'], T_out=T_r)
            self.utility_record['Q1'] = [mix_heater['Q'], mix_heater['E'], 0]
            reactor_in = mix_heater['Fo']
            # print(f"reactor_in: {reactor_in}")
            # Run the reactor simulation
            feed = pd.Series([reactor_in['T'], Pin, reactor_in.iloc[:-2].values],
                             index=['T0', 'P0', 'F0'])

            reactor_out = self.one_pass(self.reactor_para.iloc[0], self.membrane.iloc[0], feed)
            # print(f"reactor_out: {reactor_out}")
            # Split the reactor output into reactant and permeate streams
            reactant_stream = reactor_out[self.comp_num * 0 + 1:self.comp_num * 1 + 1, -1]
            T_reactant = reactor_out[self.comp_num * 4 + 1, -1]
            P_reactant = reactor_out[self.comp_num * 4 + 4, -1]

            # Cool and separate the reactant stream
            reactant_pd = pd.concat([pd.Series(reactant_stream, index=self.comp_list),
                                     pd.Series([T_reactant, P_reactant], index=["T", "P"])])

            reactant_cooler = heater(reactant_pd, T_cool)
            cool_reactant = reactant_cooler["Fo"]
            # print(f"cool_reactant: {cool_reactant}")
            flash_res = flasher(cool_reactant, T_cool=T_cool, P_cool=P_reactant)
            flash_gas, flash_liq = flash_res['Fg_o'], flash_res['Fl_o']

            split_res = spliter(flash_gas, split_ratio)
            purge_stream, recycled_reactant = split_res['purge'], split_res['recycle']

            # compress recycled_reactant if possible
            if P_reactant < Pin:
                recycle_comp_res = compressor(recycled_reactant, Pin / P_reactant)
                recycle_comped = recycle_comp_res['Fo']
                recycle_comped['P'] = np.round(recycle_comped['P'], 2)
            else:
                recycle_comped = recycled_reactant.copy()

            # heat recovery at reactor outlet
            hx = HeatExchanger(reactant_pd, recycle_comped, U=500)
            hx_res = hx.fixed_delta_cold(delta_T=10)  # hx.fixed_T_cold(Tout=T_r)
            reactant_recovery = hx_res['Fh_o']
            if reactant_recovery['T'] > T_cool:
                reactant_recovery_cooler = heater(reactant_recovery, T_out=T_cool)
            else:
                reactant_recovery_cooler = hx_res
            heater_res = hx_res['Fc_o']

            # # cooled reactant will be heated to fresh feed temperature
            # heater_res = heater(recycled_reactant, T_r)

            if self.status == 1:
                # Membrane is on
                permeate_stream = reactor_out[self.comp_num * 3 + 1:self.comp_num * 4 + 1, -1]
                T_per = reactor_out[self.comp_num * 4 + 2, -1]
                permeate_pd = pd.concat([pd.Series(permeate_stream, index=self.comp_list),
                                         pd.Series([T_per, Pp], index=["T", "P"])])

                # separate water and increase pressure
                try:
                    permeate_compressed = self.process_permeate(permeate_pd, Pin)
                except TypeError:
                    print(permeate_pd)
                    exit(0)

                # mix recycled reactant and recycled permeate
                mix_res = mixer(permeate_compressed, heater_res)  # heater_res['Fo']
                new_mix_T = mix_res['Fo']['T']
                new_mix = mix_res['Fo']
                # print(f"new_mix_F: {new_mix}")
            else:
                new_mix = heater_res  # ['Fo']
                new_mix_T = heater_res['T']  # T_r['Fo']

            new_recycle_F = new_mix.copy()
            # new_recycle_T = new_mix_T
            # Check for convergence
            # print(f"recycle_pd: {recycle_pd}")
            # print(f"new_recycle_F: {new_recycle_F}")
            diff = (new_recycle_F[self.comp_list] - recycle_pd[self.comp_list]).abs() / recycle_pd[self.comp_list]
            diff[np.isnan(diff)] = 0
            F_diff = np.max(diff)
            # T_diff = abs(new_recycle_T - T_recycle) / T_recycle
            if F_diff < tol:  # and T_diff < tol:
                break

            # Wegstein update
            if i > 0:
                wy = (new_recycle_F[self.comp_list] - recycle_pd[self.comp_list]) / \
                     (recycle_pd[self.comp_list] - recycle_prev[self.comp_list])
                wy[np.isnan(wy)] = 0
                qy = wy / (wy - 1)
                qy[np.isnan(qy)] = 0
                qy = np.clip(qy, -5, 0)
                qy["T"], qy["P"] = 0, 0
                # print(f"qy: {qy}")
                # try:
                #     wT = (new_recycle_T - T_recycle) / (T_recycle - T_recycle_prev)
                #     qT = wT / (wT - 1)
                #     qT = 0 if np.isnan(qT) else qT
                # except ZeroDivisionError:
                #     qT = 0
                # qT = -0 if qT > 0 else qT
                # qT = -4 if qT < -4 else qT
                # # print(f"qT:{qT}")
                recycle_prev = recycle_pd.copy()
                # T_recycle_prev = T_recycle

                recycle_pd = qy * recycle_pd.copy() + (1 - qy) * new_recycle_F.copy()
                # recycle = qy * recycle + (1 - qy) * new_recycle_F
                # T_recycle = qT * T_recycle + (1 - qT) * new_recycle_T
                recycle_pd['T'] = new_mix_T  # T_recycle

            if i == 0:
                recycle_prev = recycle_pd.copy()
                recycle_pd = new_recycle_F.copy()
                # T_recycle_prev = T_recycle
                # T_recycle = new_recycle_T
            recycle_pd.loc[recycle_pd < 0] = 1e-5
            print(f"diff: {diff}")
            # print(f"Tdiff: {T_diff}")
            print(f"recycle ratio: {np.sum(recycle_pd[self.comp_list]) / np.sum(F_fresh)}")
        # print(f"recycle ratio: {np.sum(recycle_pd[self.comp_list]) / np.sum(F_fresh)}")
        # record the condition for each point
        self.cond_record['S2'] = recycle_pd.copy()
        self.cond_record['S3'] = reactor_in.copy()  # recycle_mixer['Fo'].copy()
        self.cond_record['S3p'] = recycle_mixer['Fo']
        self.cond_record['S4'] = reactant_pd.copy()
        self.cond_record['S5'] = cool_reactant.copy()
        self.cond_record['S6'] = flash_gas.copy()
        self.cond_record['S7'] = recycled_reactant.copy()
        self.cond_record['S8'] = purge_stream.copy()
        self.cond_record['S9'] = flash_liq.copy()
        self.cond_record['S10'] = recycle_comped.copy()
        self.cond_record['S11'] = reactant_recovery.copy()

        # mix_heater = heater(recycle_mixer['Fo'], T_out=T_r)
        # self.utility_record['Q1'] = [mix_heater['Q'], mix_heater['E'], 0]
        self.utility_record['Q2t'] = [reactant_recovery_cooler['Q'], reactant_recovery_cooler['E'], 0]

        if self.cond_record['S7']['P'] != self.cond_record['S10']['P']:
            self.utility_record['CR'] = [0, recycle_comp_res['W'], recycle_comp_res['W']]

        if self.status == 1:
            self.pre_process_permeate(permeate_pd)
        if self.membrane.iloc[0]['qmh'] == 0:
            self.post_recovery_ad()
        else:
            self.post_recovery_iso(reactor_out)
        self.save_point()

        return reactor_out, recycle_pd  # recycle, combined_stream, purge_stream

    def sim(self, save_profile=0, loop='direct', rtol=0.05, r_target=None):
        kn_model = self.chem_para['kn_model']
        loop = "direct" if self.recycle == 0 else loop
        if self.recycle == 1:
            res_profile, F_recycle = self.recycle_loop(split_ratio=0.01, T_cool=35 + 273.15, tol=1e-3)
            xi = res_profile[1:self.comp_num + 1] / np.sum(res_profile[1:self.comp_num + 1], axis=0)
            res_profile = np.insert(res_profile, self.comp_num + 1, xi, axis=0)
            r_metric = self.reactor_metric(res_profile)
            # calculate metric for recycled reactor

            p_metric = self.recycle_metric(F_recycle)
            # p_metric = pd.concat([p_metric, p_conversion, p_conversion_H2])

            metric = pd.concat([p_metric, r_metric])
            res_path = 'result/sim_recycle_%s_%.4f_%s_log.xlsx' % (kn_model, rtol, datetime.now().date())
        else:
            res_profile = self.multi_reactor() if r_target is None else self.to_r(r_target)
            r_metric = self.reactor_metric(res_profile)
            metric = r_metric
            mode = 'ISO' if self.qmh != 0 else 'AIA'
            res_path = 'result/sim_one_pass_U_%s_%s_%s_%s_log.xlsx' % (
                kn_model, self.stage, mode, datetime.now().date())
        # calculate partial fugacity along the reactor

        print("*" * 10)
        print(metric)

        # concat the input para with the performance metrics
        feed_cond = self.feed0_para
        res = pd.concat([self.reactors_para, self.membranes, feed_cond, metric])
        res_save = pd.DataFrame(res.values.reshape(1, len(res.values)), columns=res.index)

        # save data to the Excel
        try:
            with pd.ExcelWriter(res_path, engine='openpyxl', mode='a', if_sheet_exists="overlay") as writer:
                try:
                    res_saved = pd.read_excel(res_path, sheet_name=loop)
                    res_save = pd.concat([res_saved, res_save], ignore_index=True)

                    res_save.to_excel(writer, index=False, header=True, sheet_name=loop)
                except ValueError:
                    res_save.to_excel(writer, index=False, header=True, sheet_name=loop)
        except FileNotFoundError:
            res_save.to_excel(res_path, index=False, header=True, sheet_name=loop)
        if save_profile == 1:
            save_data = pd.DataFrame(res_profile.T, columns=['z'] + self.comp_list +
                                                            ['xi_' + i for i in self.comp_list] +
                                                            ['dF_re_' + i for i in self.comp_list] +
                                                            ['dF_diff_' + i for i in self.comp_list] +
                                                            ['dF_per_' + i for i in self.comp_list] +
                                                            ['Tr', 'Tc', 'Th', 'dP',
                                                             'q_react', 'q_diff', 'q_heater', 'h_diff'])

            sim_path = 'result/sim_profile1_U_%s_%s_%s_%s_%s_%s_%s_%s.xlsx' \
                       % (self.stage, kn_model, self.Dt.values, self.L.values,
                          self.T0, self.P0, self.Tc, self.Th)

            sheet_name = 'U_%s_eos_%s_drop_%s' % (self.Uc, self.eos, self.drop) if self.recycle == 0 else \
                f'U_{self.Uc}_drop_{self.drop}_rtol_{rtol}'
            try:
                with pd.ExcelWriter(sim_path, engine='openpyxl', mode='a', if_sheet_exists='new') as writer:
                    save_data.to_excel(writer, index=False, sheet_name=sheet_name)
            except FileNotFoundError:
                save_data.to_excel(sim_path, index=False, sheet_name=sheet_name)
        return res_save
