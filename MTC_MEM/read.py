import os.path
import sys
import os
import json
import logging
import numpy as np
import pandas as pd

R = 8.314


class ReadData:
    """
    read data for simulation of conversion from CO2 to CH3OH
    """

    def __init__(self, kn_model='BU', in_path=None):
        current_dir = os.path.dirname(os.path.abspath(__file__))

        # chem data of reaction for kinetic model
        if kn_model == 'BU':
            self.chem = self.kn_bu
        elif kn_model == 'SL':
            self.chem = self.kn_sl
        elif kn_model == 'GR':
            self.chem = self.kn_gr
        else:
            print('kn_model should be one of BU, GR, and SL')
            sys.exit(1)
        self.chem['kn_model'] = kn_model

        # data path for reactor and feed gas
        self.root_path = sys.path[0]
        in_reactor_path = os.path.join(current_dir, 'in_reactor.json')
        feed_path = os.path.join(current_dir, 'in_feed.json')
        in_mem_path = os.path.join(current_dir, 'in_mem.json')
        if in_path is None:
            in_path = {'reactor': in_reactor_path, 'feed': feed_path,
                       'membrane': in_mem_path}
        in_data = dict()
        for key, values in in_path.items():
            try:
                file_path = os.path.join(self.root_path, values)
                with open(file_path) as f:
                    in_data[key] = json.load(f)
            except FileNotFoundError:
                with open(values) as f:
                    in_data[key] = json.load(f)

        # reactor parameters
        self.react_para = in_data["reactor"]["reactor"]

        # feed gas parameter
        self.feed_para = in_data['feed']

        # membrane parameter
        self.mem_para = in_data["membrane"]["membrane"]

    @property
    def feed_data(self):
        T0_CO2 = self.feed_para["condition"]["T0_CO2"]
        P0_CO2 = self.feed_para["condition"]["P0_CO2"]
        T0_H2 = self.feed_para["condition"]["T0_H2"]
        P0_H2 = self.feed_para["condition"]["P0_H2"]
        if self.feed_para["condition"]["fresh"] == 'on':
            feed = self.feed_para["condition"]
            H2_CO2 = feed["H2/CO2"]
            inert = feed["inert"]
            recycle = 1

            # feed para frame
            T0_array = self.data_array(feed["T"])
            P0_array = self.data_array(feed["P"])
            sv_array = self.data_array(feed["Sv"])
            H2_array = self.data_array(feed['H2'])
            CO_CO2_array = self.data_array(feed['CO/CO2'])
            feed_num = len(T0_array) * len(P0_array) * len(sv_array) * len(CO_CO2_array)
            feed_para = pd.DataFrame(index=np.arange(feed_num), columns=list(feed.keys()))
            i = 0
            for T in T0_array:
                for P in P0_array:
                    for sv in sv_array:
                        for CO_CO2 in CO_CO2_array:
                            for H2 in H2_array:
                                feed_para.iloc[i] = np.array([T, P, recycle, H2_CO2, CO_CO2, H2, sv, inert,
                                                              T0_CO2, T0_H2, P0_CO2, P0_H2])
                                i += 1
        else:
            feed = np.array([float(i) for i in self.feed_para["feed"].split('\t')])
            T = self.feed_para["condition"]['T'][0]
            P = self.feed_para["condition"]['P'][0]
            feed = np.append(np.array([T, P, 0]), feed)
            feed_para = pd.DataFrame(feed.reshape(1, len(feed)),
                                     columns=['T', 'P', 'fresh', "CO2", "H2", "Methanol", "H2O", "CO", "N2"])
        return feed_para

    @property
    def mem_data(self):
        stage = self.react_para['stage']

        paras_stage = ['Din', 'Thick', 'Tc', 'qmc', 'Th', 'qmh', 'Pp', 'per_H2O', 'per_S']
        paras_array = {'status': [0 if i == 'off' else 1 for i in self.mem_para['status']],
                       'pattern': self.mem_para["pattern"], 'pattern_h': self.mem_para["pattern_h"],
                       }
        paras_array_name = {'status_name': [f'status{n + 1}' for n in range(stage)],
                            'pattern_name': [f'pattern{n + 1}' for n in range(stage)],
                            'pattern_h_name': [f'pattern_h{n + 1}' for n in range(stage)]}
        for n in range(stage):
            for para_stage in paras_stage:
                paras_array_name[f'{para_stage}_name'] = [f'{para_stage}{n + 1}' for n in range(stage)]
                if n == 0:
                    paras_array[f'{para_stage}_array'] = [self.data_array(self.mem_para[f'{para_stage}'][n])]
                else:
                    paras_array[f'{para_stage}_array'].append(self.data_array(self.mem_para[f'{para_stage}'][n]))
        location = 0 if self.mem_para["io"] == 'in' else 1  # in:reaction in the insides of membrane
        nit = self.mem_para["nit"]  # tube number of the insulator
        q = self.mem_para['q']
        Fp = np.array([float(i) for i in self.mem_para['Fp'].split('\t')])
        per = self.mem_para['per']

        # insulator para frame
        # generate the combination of reactor length, diameter

        paras_array_comb = {}
        mem_num = 1
        for para_stage in paras_stage:
            if len(paras_array[f'{para_stage}_array']) >= stage:
                combinations = [[p] for p in paras_array[f'{para_stage}_array'][0]]
                for i in range(1, stage):
                    new_combinations = []
                    for c in combinations:
                        for p in paras_array[f'{para_stage}_array'][i]:
                            new_combinations.append(c + [p])
                    combinations = new_combinations
                paras_array_comb[f'{para_stage}'] = combinations
            else:
                paras_array_comb[f'{para_stage}'] = []  # 或者设置为空列表，根据需要处理
            mem_num *= len(paras_array_comb[f'{para_stage}'])
        column_name = []
        for para_name in paras_array_name.values():
            column_name += para_name
        column_name += ['io', 'nit', 'q', 'Fp', 'per']

        mem_para = pd.DataFrame(index=np.arange(mem_num), columns=column_name)
        i = 0
        for Din in paras_array_comb['Din']:
            for thick in paras_array_comb['Thick']:
                thick = [round(k, 3) for k in thick]
                for Tc in paras_array_comb['Tc']:
                    for qmc in paras_array_comb['qmc']:
                        for Th in paras_array_comb['Th']:
                            for qmh in paras_array_comb['qmh']:
                                for P in paras_array_comb["Pp"]:
                                    for per_H2O in paras_array_comb['per_H2O']:
                                        for per_S in paras_array_comb['per_S']:
                                            mem_para.iloc[i] = paras_array['status'] + \
                                                               paras_array['pattern'] + paras_array['pattern_h'] + \
                                                               Din + thick + \
                                                               Tc + qmc + Th + qmh + P + per_H2O + per_S + \
                                                               [location, nit, q, Fp, per]
                                            i += 1
        return mem_para

    @property
    def reactor_data(self):
        nrt = self.react_para['nrt']  # number of the reaction tube
        phi = self.react_para["phi"]  # void of fraction
        rhoc = self.react_para["rhoc"]  # density of catalyst, kg/m3
        stage = self.react_para['stage']
        # L2 = self.react_para['L2']  # length, m
        Uc = self.react_para["Uc"]  # total heat transfer coefficient of the reactor, W/m2 K
        recycle = 1 if self.react_para['recycle'] == "on" else 0

        # reactor para frame
        L_array, Dt_array = [], []
        Din_array = []
        for n in range(stage):
            L_array.append(self.data_array(self.react_para['L'][n]))
            Dt_array.append(self.data_array(self.react_para['Dt'][n]))
            Din_array.append(self.data_array(self.react_para['Dc'][n]))

        column = list(self.react_para.keys())
        column.remove('L')
        column.remove('Dt')
        column.remove('Dc')
        Dt_name = [f'Dt{n + 1}' for n in range(stage)]
        L_name = [f'L{n + 1}' for n in range(stage)]
        Din_name = [f'Dc{n + 1}' for n in range(stage)]
        # generate the combination of reactor length, diameter
        if stage > 0:
            Ls = [[L1] for L1 in L_array[0]]
            Dts = [[Dt1] for Dt1 in Dt_array[0]]
            Dins = [[Din1] for Din1 in Din_array[0]]
        if stage > 1:
            for i in range(1, stage):
                Ls = [[*L, L_i] for L in Ls for L_i in L_array[i]]
                Dts = [[*Dt, Dt_i] for Dt in Dts for Dt_i in Dt_array[i]]
                Dins = [[*Din, Din_i] for Din in Dins for Din_i in Din_array[i]]
        reactor_num = len(Ls) * len(Dts) * len(Dins)
        react_para = pd.DataFrame(index=np.arange(reactor_num), columns=L_name + Dt_name + Din_name + column)
        i = 0
        for L in Ls:
            for Dt in Dts:
                for Din in Dins:
                    react_para.iloc[i] = L + Dt + Din + [stage, nrt, rhoc, phi, recycle, Uc]
                    i += 1
        return react_para

    @staticmethod
    def data_array(in_data):
        try:
            data_lenth = len(in_data)
        except TypeError:
            print("Value in json should be list!")
            sys.exit(1)
        if data_lenth != 3:
            out_data = np.array(in_data)
        elif data_lenth == 3:
            out_data = np.linspace(in_data[0], in_data[1], in_data[2])
        return out_data

    @property
    def hr(self):
        """
        parameter for the calculation of reaction enthalpy, dH = aT^4+b J/kmol
        ref: Cui, 2020, Chemical Engineering Journal, 10.1016/j.cej.2020.124632
        :return: [a, b] for reactions
        """
        heat_reaction = {
            "1": [3.589e4, 4.0047e7],
            "2": [9.177e3, -4.4325e7]
        }
        return heat_reaction

    @property
    def keq(self):
        """
        parameter for the calculation of equilibrium constant, k = 10^(a/T+b)
        ref: Graaf, 1986, Chemical Engineering Science, 10.1016/0009-2509(86)80019-7
        :return: [a, b] for reactions
        """
        keq = {
            "1": [3066, -10.592],
            "2": [-2073, 2.029]
        }
        # keq = {
        #     "1": [3066, -10.92],
        #     "2": [-2073, 2.029]
        # }
        return keq

    @property
    def kn_sl(self):
        """
        reaction kinetic model proposed by Slotboom
        ref: Slotboom, 2020, Chemical Engineering Journal, 10.1016/j.cej.2020.124181
        :return: parameter dict
        """
        stoichiometry = {
            "1": [-1, -3, 1, 1, 0],
            "2": [-1, -1, 0, 1, 1]
        }
        kad = {
            "H2": [1.099, 0],
            "H2O": [126.4, 0]
        }
        kr = {
            "1": [7.414e14, -166000],
            "2": [1.111e19, -203700]
        }
        chem_data = {"stoichiometry": stoichiometry, "kad": kad, "kr": kr, "keq": self.keq, 'heat_reaction': self.hr}
        return chem_data

    @property
    def kn_bu(self):
        """
        reaction kinetic model proposed by Bussche
        ref: Bussche, 1996, Journal of Catalysis, 10.1006/jcat.1996.0156
        :return: parameter dict
        """
        stoichiometry = {
            "1": [-1, -3, 1, 1, 0],
            "2": [-1, -1, 0, 1, 1]
        }
        # kad = {
        #     "H2": [0.499, 17197],
        #     "H2O": [6.62e-11, 124119],
        #     "H2O/H2": [3453.38, 0]
        # }
        kad = {
            "H2": [0.7089, 8526],
            "H2O": [1.587e-7, 93410],
            "H2O/H2": [3322, 0]
        }
        # kr = {
        #     "1": [1.07, 36696],
        #     "2": [1.22e10, -94765]
        # }
        # kr = {
        #     "1": [1.07, 40000],
        #     "2": [1.22e10, -98084]
        # }
        kr = {
            "1": [0.1749, 44450],
            "2": [1.446e12, -122600]
        }
        # kad = {
        #     "H2": [0.7089, 8526],
        #     "H2O": [1.587e-7, 93410],
        #     "H2O/H2": [3322, 0]
        # }
        chem_data = {"stoichiometry": stoichiometry, "kad": kad, "kr": kr, "keq": self.keq, 'heat_reaction': self.hr}
        return chem_data

    @property
    def kn_gr(self):
        """
        reaction kinetic model proposed by graaf
        ref: Graaf, 1988, Chemical Engineering Science, 10.1016/0009-2509(88)85127-3
        :return: parameter dict
        """
        stoichiometry = {
            "1": [-1, -3, 1, 1, 0],
            "2": [-1, -1, 0, 1, 1],
            "3": [0, -2, 1, 0, -1]
        }
        kad = {
            "CO": [7.99e-7, 58100],
            "CO2": [1.02e-7, 67400],
            "H2O/H2": [4.13e-11, 104500]
        }
        kr = {
            "1": [436, -65200],
            "2": [7.31e8, -123400],
            "3": [2.69e7, -109900]
        }
        chem_data = {"stoichiometry": stoichiometry, "kad": kad, "kr": kr, "keq": self.keq, 'heat_reaction': self.hr}
        return chem_data

# data = ReadData(kn_model='BU')
# # print(data.feed_para.keys())
# # print(pd.DataFrame(columns=data.feed_para.keys()))
# print(data.feed_data)
# print(data.insulator_data)
# print(data.reactor_data)
# print()
# print(data.kn_bu["kr"])
# react_rate_constant = {}
#
# for key, value in data.kn_bu["kr"].items():
#     react_rate_constant[key] = value[0] * np.exp(value[1] / 503 / R)
#
# print(react_rate_constant)
