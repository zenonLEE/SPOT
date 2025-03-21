import time

import numpy as np
import math

import pandas as pd
import win32com.client
from CoolProp.CoolProp import PropsSI

from prop_calculator import VLEThermo


class HeatExchanger:
    def __init__(self, hot_fluid, cold_fluid, U=500):
        self.hot_fluid = hot_fluid
        self.cold_fluid = cold_fluid
        self.U = U

        # Create VLEThermo instances
        self.hot_species = hot_fluid.index.tolist()[:-2]
        self.cold_species = cold_fluid.index.tolist()[:-2]
        self.thermo_hot = VLEThermo(self.hot_species)
        self.thermo_cold = VLEThermo(self.cold_species)

    def fixed_delta_hot(self, delta_T):
        # Calculate enthalpies using VLEThermo
        h_hot_in = self.thermo_hot.cal_H(self.hot_fluid['T'], self.hot_fluid['P'],
                                         self.hot_fluid[self.hot_species].values)
        h_cold_in = self.thermo_cold.cal_H(self.cold_fluid['T'], self.cold_fluid['P'],
                                           self.cold_fluid[self.hot_species].values)

        # Determine hot fluid outlet temperature
        T_hot_out = self.cold_fluid["T"] + delta_T

        h_hot_out = self.thermo_hot.cal_H(T_hot_out, self.hot_fluid["P"], self.hot_fluid[self.hot_species].values)

        # Calculate heat duty (Q)
        Q = (h_hot_in - h_hot_out)

        # Determine cold fluid outlet enthalpy and temperature
        h_cold_out = h_cold_in + Q
        T_cold_out = self.thermo_cold.cal_T_from_H(h_cold_out, self.cold_fluid['P'],
                                                   self.cold_fluid[self.cold_species].values)

        # Logarithmic Mean Temperature Difference (LMTD)
        delta_T1 = self.hot_fluid["T"] - T_cold_out
        delta_T2 = T_hot_out - self.cold_fluid["T"]
        LMTD = (delta_T1 - delta_T2) / math.log(delta_T1 / delta_T2)

        # Surface area (A)
        A = Q / (self.U * LMTD)

        # Returning the results
        Fh_o, Fc_o = self.hot_fluid.copy(), self.cold_fluid.copy()
        Fh_o['T'] = T_hot_out
        Fc_o["T"] = T_cold_out

        return {
            'Fh_in': self.hot_fluid,
            'Fh_o': Fh_o,
            'Fc_in': self.cold_fluid,
            'Fc_o': Fc_o,
            'Q': Q,
            'A': A
        }

    def fixed_delta_cold(self, delta_T):
        # known the temperature diff between hot in and cold out
        # Calculate enthalpies using VLEThermo
        h_hot_in = self.thermo_hot.cal_H(self.hot_fluid['T'], self.hot_fluid['P'],
                                         self.hot_fluid[self.hot_species].values)
        h_cold_in = self.thermo_cold.cal_H(self.cold_fluid['T'], self.cold_fluid['P'],
                                           self.cold_fluid[self.hot_species].values)

        # Determine cold fluid outlet temperature
        T_cold_out = self.hot_fluid["T"] - delta_T
        h_cold_out = self.thermo_cold.cal_H(T_cold_out, self.cold_fluid["P"], self.cold_fluid[self.cold_species].values)

        # Calculate heat duty (Q)
        Q = (h_cold_in - h_cold_out)

        # Determine cold fluid outlet enthalpy and temperature
        h_hot_out = h_hot_in + Q
        T_hot_out = self.thermo_hot.cal_T_from_H(h_hot_out, self.hot_fluid['P'],
                                                 self.hot_fluid[self.hot_species].values)

        # Logarithmic Mean Temperature Difference (LMTD)
        delta_T1 = self.hot_fluid["T"] - T_cold_out
        delta_T2 = T_hot_out - self.cold_fluid["T"]
        LMTD = (delta_T1 - delta_T2) / math.log(delta_T1 / delta_T2)

        # Surface area (A)
        A = -Q / (self.U * LMTD)

        # Returning the results
        Fh_o, Fc_o = self.hot_fluid.copy(), self.cold_fluid.copy()
        Fh_o['T'] = T_hot_out
        Fc_o["T"] = T_cold_out

        return {
            'Fh_in': self.hot_fluid,
            'Fh_o': Fh_o,
            'Fc_in': self.cold_fluid,
            'Fc_o': Fc_o,
            'Q': Q,
            'A': A
        }

    def fixed_duty(self, duty):
        # Calculate enthalpies using VLEThermo
        h_hot_in = self.thermo_hot.cal_H(self.hot_fluid['T'], self.hot_fluid['P'],
                                         self.hot_fluid[self.hot_species].values)
        h_cold_in = self.thermo_cold.cal_H(self.cold_fluid['T'], self.cold_fluid['P'],
                                           self.cold_fluid[self.cold_species].values)

        # Determine cold fluid outlet enthalpy and temperature
        h_cold_out = h_cold_in + duty
        T_cold_out = self.thermo_cold.cal_T_from_H(h_cold_out, self.cold_fluid['P'],
                                                   self.cold_fluid[self.cold_species].values)

        # Determine hot fluid outlet enthalpy and temperature
        h_hot_out = h_hot_in - duty
        T_hot_out = self.thermo_hot.cal_T_from_H(h_hot_out, self.hot_fluid['P'],
                                                 self.hot_fluid[self.hot_species].values)

        # Logarithmic Mean Temperature Difference (LMTD)
        delta_T1 = self.hot_fluid["T"] - T_cold_out
        delta_T2 = T_hot_out - self.cold_fluid["T"]
        LMTD = (delta_T1 - delta_T2) / math.log(delta_T1 / delta_T2)

        # exergy input
        E_in = self.thermo_hot.cal_E(self.hot_fluid['T'], self.hot_fluid['P'], self.hot_fluid[self.hot_species].values) \
               - self.thermo_hot.cal_E(T_hot_out, self.hot_fluid['P'], self.hot_fluid[self.hot_species].values)

        # Surface area (A)
        A = duty / (self.U * LMTD)

        # Returning the results
        Fh_o, Fc_o = self.hot_fluid.copy(), self.cold_fluid.copy()
        Fh_o['T'] = T_hot_out
        Fc_o["T"] = T_cold_out

        return {
            'Fh_in': self.hot_fluid,
            'Fh_o': Fh_o,
            'Fc_in': self.cold_fluid,
            'Fc_o': Fc_o,
            'Q': duty,
            'E': E_in,
            'A': A
        }

    def fixed_delta_hc(self, delta_T):
        """
        specify the temperature difference between hot fluid output and cold fluid output
        :param delta_T:
        :return:
        """
        # Calculate enthalpies using VLEThermo
        h_hot_in = self.thermo_hot.cal_H(self.hot_fluid['T'], self.hot_fluid['P'],
                                         self.hot_fluid[self.hot_species].values)
        h_cold_in = self.thermo_cold.cal_H(self.cold_fluid['T'], self.cold_fluid['P'],
                                           self.cold_fluid[self.hot_species].values)
        if self.cold_fluid['T'] + 10 > self.hot_fluid['T']:
            print("heat integration fails in fixed_delta_hc")
            raise ValueError
        T_hot_outs = np.arange(self.cold_fluid['T'] + 1, self.hot_fluid['T'], 0.1)

        # Determine cold fluid outlet temperature
        for T_hot_out in T_hot_outs:
            h_hot_out = self.thermo_hot.cal_H(T_hot_out, self.hot_fluid["P"], self.hot_fluid[self.hot_species].values)

            # Calculate heat duty (Q)
            Q = (h_hot_in - h_hot_out)

            # Determine cold fluid outlet enthalpy and temperature
            h_cold_out = h_cold_in + Q
            try:
                T_cold_out = self.thermo_cold.cal_T_from_H(h_cold_out, self.cold_fluid['P'],
                                                           self.cold_fluid[self.cold_species].values)
            except TypeError:
                continue

            T_delta_hc = T_hot_out - T_cold_out

            if delta_T > T_delta_hc > 0:
                break

        # Logarithmic Mean Temperature Difference (LMTD)
        delta_T1 = self.hot_fluid["T"] - T_cold_out
        delta_T2 = T_hot_out - self.cold_fluid["T"]
        LMTD = (delta_T1 - delta_T2) / math.log(delta_T1 / delta_T2)

        # Surface area (A)
        A = Q / (self.U * LMTD)

        # Returning the results
        Fh_o, Fc_o = self.hot_fluid.copy(), self.cold_fluid.copy()
        Fh_o['T'] = T_hot_out
        Fc_o["T"] = T_cold_out

        return {
            'Fh_in': self.hot_fluid,
            'Fh_o': Fh_o,
            'Fc_in': self.cold_fluid,
            'Fc_o': Fc_o,
            'Q': Q,
            'A': A
        }

    def fixed_T_hot(self, Tout):
        # Calculate enthalpies using VLEThermo
        h_hot_in = self.thermo_hot.cal_H(self.hot_fluid['T'], self.hot_fluid['P'],
                                         self.hot_fluid[self.hot_species].values)
        h_cold_in = self.thermo_cold.cal_H(self.cold_fluid['T'], self.cold_fluid['P'],
                                           self.cold_fluid[self.hot_species].values)

        # Determine hot fluid outlet temperature
        T_hot_out = Tout
        h_hot_out = self.thermo_hot.cal_H(T_hot_out, self.hot_fluid['P'],
                                          self.hot_fluid[self.hot_species].values)

        # Calculate heat duty (Q)
        Q = (h_hot_in - h_hot_out)
        # Determine cold fluid outlet enthalpy and temperature
        h_cold_out = h_cold_in + Q

        T_cold_out = self.thermo_cold.cal_T_from_H(h_cold_out, self.cold_fluid['P'],
                                                   self.cold_fluid[self.cold_species].values)
        # Logarithmic Mean Temperature Difference (LMTD)
        delta_T1 = self.hot_fluid["T"] - T_cold_out
        delta_T2 = T_hot_out - self.cold_fluid["T"]
        LMTD = (delta_T1 - delta_T2) / math.log(delta_T1 / delta_T2)

        # Surface area (A)
        A = Q / (self.U * LMTD)

        # Returning the results
        Fh_o, Fc_o = self.hot_fluid.copy(), self.cold_fluid.copy()
        Fh_o['T'] = T_hot_out
        Fc_o["T"] = T_cold_out

        return {
            'Fh_in': self.hot_fluid,
            'Fh_o': Fh_o,
            'Fc_in': self.cold_fluid,
            'Fc_o': Fc_o,
            'Q': Q,
            'A': A
        }

    def fixed_T_cold(self, Tout):
        # Calculate enthalpies using VLEThermo
        h_hot_in = self.thermo_hot.cal_H(self.hot_fluid['T'], self.hot_fluid['P'],
                                         self.hot_fluid[self.hot_species].values)
        h_cold_in = self.thermo_cold.cal_H(self.cold_fluid['T'], self.cold_fluid['P'],
                                           self.cold_fluid[self.hot_species].values)

        # Determine cold fluid outlet temperature
        T_cold_out = Tout
        h_cold_out = self.thermo_cold.cal_H(T_cold_out, self.cold_fluid['P'],
                                            self.cold_fluid[self.hot_species].values)

        # Calculate heat duty (Q)
        Q = (h_cold_in - h_cold_out)

        # Determine cold fluid outlet enthalpy and temperature
        h_hot_out = h_hot_in + Q

        T_hot_out = self.thermo_hot.cal_T_from_H(h_hot_out, self.hot_fluid['P'],
                                                 self.hot_fluid[self.hot_species].values)

        # Logarithmic Mean Temperature Difference (LMTD)
        delta_T1 = self.hot_fluid["T"] - T_cold_out
        delta_T2 = T_hot_out - self.cold_fluid["T"]
        LMTD = (delta_T1 - delta_T2) / math.log(delta_T1 / delta_T2)

        # Surface area (A)
        A = -Q / (self.U * LMTD)

        # Returning the results
        Fh_o, Fc_o = self.hot_fluid.copy(), self.cold_fluid.copy()
        Fh_o['T'] = T_hot_out
        Fc_o["T"] = T_cold_out

        return {
            'Fh_in': self.hot_fluid,
            'Fh_o': Fh_o,
            'Fc_in': self.cold_fluid,
            'Fc_o': Fc_o,
            'Q': Q,
            'A': A
        }


class HX4Water:
    def __init__(self, hot_fluid, cold_fluid=None, U=500):
        self.hot_fluid = hot_fluid.copy()  # water
        self.cold_fluid = cold_fluid
        self.U = U

        # Create VLEThermo instances

        self.cold_species = cold_fluid.index.tolist()[:-2]
        self.thermo_cold = VLEThermo(self.cold_species)

    def fixed_delta_hot(self, delta_T):
        # Calculate enthalpies using VLEThermo

        h_hot_in = self.thermo_hot.cal_H(self.hot_fluid['T'], self.hot_fluid['P'],
                                         self.hot_fluid[self.hot_species].values)
        h_cold_in = self.thermo_cold.cal_H(self.cold_fluid['T'], self.cold_fluid['P'],
                                           self.cold_fluid[self.hot_species].values)

        # Determine hot fluid outlet temperature
        T_hot_out = self.cold_fluid["T"] + delta_T
        print(T_hot_out)
        h_hot_out = self.thermo_hot.cal_H(T_hot_out, self.hot_fluid["P"], self.hot_fluid[self.hot_species].values)

        # Calculate heat duty (Q)
        Q = (h_hot_in - h_hot_out)

        # Determine cold fluid outlet enthalpy and temperature
        h_cold_out = h_cold_in + Q
        T_cold_out = self.thermo_cold.cal_T_from_H(h_cold_out, self.cold_fluid['P'],
                                                   self.cold_fluid[self.cold_species].values)
        print(T_cold_out)

        # Logarithmic Mean Temperature Difference (LMTD)
        delta_T1 = self.hot_fluid["T"] - T_cold_out
        delta_T2 = T_hot_out - self.cold_fluid["T"]
        LMTD = (delta_T1 - delta_T2) / math.log(delta_T1 / delta_T2)

        # Surface area (A)
        A = Q / (self.U * LMTD)

        # Returning the results
        Fh_o, Fc_o = self.hot_fluid.copy(), self.cold_fluid.copy()
        Fh_o['T'] = T_hot_out
        Fc_o["T"] = T_cold_out

        return {
            'Fh_in': self.hot_fluid,
            'Fh_o': Fh_o,
            'Fc_in': self.cold_fluid,
            'Fc_o': Fc_o,
            'Q': Q,
            'A': A
        }

    def fixed_delta_cold(self, delta_T):
        # Calculate enthalpies using VLEThermo
        h_hot_in = self.thermo_hot.cal_H(self.hot_fluid['T'], self.hot_fluid['P'],
                                         self.hot_fluid[self.hot_species].values)
        h_cold_in = self.thermo_cold.cal_H(self.cold_fluid['T'], self.cold_fluid['P'],
                                           self.cold_fluid[self.hot_species].values)

        # Determine cold fluid outlet temperature
        T_cold_out = self.hot_fluid["T"] - delta_T
        h_cold_out = self.thermo_cold.cal_H(T_cold_out, self.cold_fluid["P"], self.cold_fluid[self.cold_species].values)

        # Calculate heat duty (Q)
        Q = (h_cold_in - h_cold_out)

        # Determine cold fluid outlet enthalpy and temperature
        h_hot_out = h_hot_in + Q
        T_hot_out = self.thermo_hot.cal_T_from_H(h_hot_out, self.hot_fluid['P'],
                                                 self.hot_fluid[self.hot_species].values)

        # Logarithmic Mean Temperature Difference (LMTD)
        delta_T1 = self.hot_fluid["T"] - T_cold_out
        delta_T2 = T_hot_out - self.cold_fluid["T"]
        LMTD = (delta_T1 - delta_T2) / math.log(delta_T1 / delta_T2)

        # Surface area (A)
        A = Q / (self.U * LMTD)

        # Returning the results
        Fh_o, Fc_o = self.hot_fluid.copy(), self.cold_fluid.copy()
        Fh_o['T'] = T_hot_out
        Fc_o["T"] = T_cold_out

        return {
            'Fh_in': self.hot_fluid,
            'Fh_o': Fh_o,
            'Fc_in': self.cold_fluid,
            'Fc_o': Fc_o,
            'Q': Q,
            'A': A
        }

    def fixed_T_hot(self, Tout):
        # Calculate enthalpies using VLEThermo
        cp_hot_in = PropsSI('CPMOLAR', 'T', self.hot_fluid['T'], 'P', self.hot_fluid['P'] * 1E5, 'water')
        # Determine hot fluid outlet temperature
        T_hot_out = Tout
        cp_hot_out = PropsSI('CPMOLAR', 'T', T_hot_out, 'P', self.hot_fluid['P'] * 1E5, 'water')
        # Calculate heat duty (Q)
        Q = (self.hot_fluid['T'] - T_hot_out) * (cp_hot_in + cp_hot_out) / 2 * self.hot_fluid['qm'] / 1000

        h_cold_in = self.thermo_cold.cal_H(self.cold_fluid['T'], self.cold_fluid['P'],
                                           self.cold_fluid[self.cold_species].values)

        # Determine cold fluid outlet enthalpy and temperature
        h_cold_out = h_cold_in + Q
        T_cold_out = self.thermo_cold.cal_T_from_H(h_cold_out, self.cold_fluid['P'],
                                                   self.cold_fluid[self.cold_species].values)
        print(T_cold_out)
        # Logarithmic Mean Temperature Difference (LMTD)
        delta_T1 = self.hot_fluid["T"] - T_cold_out
        delta_T2 = T_hot_out - self.cold_fluid["T"]
        LMTD = (delta_T1 - delta_T2) / math.log(delta_T1 / delta_T2)

        # Surface area (A)
        A = Q / (self.U * LMTD)

        # Returning the results
        Fh_o, Fc_o = self.hot_fluid.copy(), self.cold_fluid.copy()
        Fh_o['T'] = T_hot_out
        Fc_o["T"] = T_cold_out

        return {
            'Fh_in': self.hot_fluid,
            'Fh_o': Fh_o,
            'Fc_in': self.cold_fluid,
            'Fc_o': Fc_o,
            'Q': Q,
            'A': A
        }

    def fixed_T_cold(self, Tout):
        # Calculate enthalpies using VLEThermo
        h_cold_in = self.thermo_cold.cal_H(self.cold_fluid['T'], self.cold_fluid['P'],
                                           self.cold_fluid[self.cold_species].values)

        # Determine cold fluid outlet temperature
        T_cold_out = Tout
        h_cold_out = self.thermo_cold.cal_H(T_cold_out, self.cold_fluid['P'],
                                            self.cold_fluid[self.cold_species].values)

        # Calculate heat duty (Q)
        Q = (h_cold_in - h_cold_out)
        # Determine hot fluid outlet temperature
        T_hot_in = self.hot_fluid['T']
        cp_hot_in = PropsSI('CPMOLAR', 'T', self.hot_fluid['T'], 'P', self.hot_fluid['P'] * 1E5, 'water')
        Q_diff = 1000
        for T in np.arange(T_hot_in - 20, T_hot_in, 0.01):
            cp_hot_out = PropsSI('CPMOLAR', 'T', T, 'P', self.hot_fluid['P'] * 1E5, 'water')
            cp_hot_m = (cp_hot_out + cp_hot_in) / 2
            Q_hot_delta = cp_hot_m * (T_hot_in - T) * self.hot_fluid['qm']
            Q_diff_cal = abs(Q_hot_delta - Q) / abs(Q)
            if Q_diff_cal < Q_diff:
                Q_diff = Q_diff_cal
                T_hot_out = T
            if Q_diff_cal < 1E-5:
                break

        # Logarithmic Mean Temperature Difference (LMTD)
        delta_T1 = self.hot_fluid["T"] - T_cold_out
        delta_T2 = T_hot_out - self.cold_fluid["T"]
        LMTD = (delta_T1 - delta_T2) / math.log(delta_T1 / delta_T2)

        # Surface area (A)
        A = Q / (self.U * LMTD)

        # Returning the results
        Fh_o, Fc_o = self.hot_fluid.copy(), self.cold_fluid.copy()
        Fh_o['T'] = T_hot_out
        Fc_o["T"] = T_cold_out

        return {
            'Fh_in': self.hot_fluid,
            'Fh_o': Fh_o,
            'Fc_in': self.cold_fluid,
            'Fc_o': Fc_o,
            'Q': Q,
            'A': A
        }


class DistillerOpt:
    def __init__(self, apk_name, block_name, feed_name, light_name, heavy_name):

        self.block = block_name
        self.feed = feed_name
        self.light = light_name
        self.heavy = heavy_name
        self._init_apk(apk_name)

        self.apk_comp = self._get_comp()

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

    def _get_comp(self):
        try:
            stream_node = self.aspen.Tree.FindNode(f"\\Data\\Streams\\{self.feed}\\Output")
            CompoundLister = stream_node.Elements("STR_MAIN").Elements("MASSFLOW").Elements("MIXED").Elements
            CompoundNameList = []
            for compound in CompoundLister:
                Compoundname = compound.Name
                CompoundNameList.append(Compoundname)
            return CompoundNameList
        except Exception as e:
            print(f"An error occurred during get comps: {e}")
            self.aspen.Quit()
            exit(1)

    def set_feed_parameters(self, stream_name=None, stream_source=None):
        """
        Set the feed parameters for the specified stream in Aspen Plus.

        :param stream_name: The name of the stream to be modified.
        :param stream_source: The new temperature (in K) for the stream, if provided.
        """
        stream_name = self.feed if stream_name is None else stream_name

        try:
            # Access the feed stream node
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
            stream_node = self.aspen.Tree.FindNode(f"\\Data\\Streams\\{stream_name}\\Input")
            stream_node.Elements("TEMP").Elements("MIXED").Value = stream['T']
            stream_node.Elements("PRES").Elements("MIXED").Value = stream['P']
            for comp in self.apk_comp:

                comp_node = stream_node.Elements("FLOW").Elements("MIXED").Elements(comp)
                try:
                    comp_node.Value = stream[comp]
                except KeyError:
                    comp_node.Value = 0

        except Exception as e:
            print(f"Error setting feed parameters: {e}")
            self.aspen.Quit()
            exit(1)

    def set_block_parameters(self, block_para=None, block_name=None):
        """
        Set the feed parameters for the specified stream in Aspen Plus.

        :param block_name: The name of the stream to be modified.
        :param stream_source: The new temperature (in K) for the stream, if provided.
        """
        block_name = self.block if block_name is None else block_name
        block_para = self.opt_distiller() if block_para is None else block_para
        feed_name_rf = 'S6'
        try:
            # Access the feed stream node
            stream_node = self.aspen.Tree.FindNode(f"\\Data\\Blocks\\{block_name}\\Input")
            stream_node.Elements("NSTAGE").Value = block_para['nstage']
            stream_node.Elements("BASIS_RR").Value = block_para['RR']
            stream_node.Elements("FEED_STAGE").Elements(feed_name_rf).Value = block_para['FS']
            stream_node.Elements("D:F").Value = block_para['DF']

        except Exception as e:
            print(f"Error setting block parameters: {e}")
            self.aspen.Quit()
            exit(1)

    def access_to_rf(self, block_name=None):
        """
        Set the feed parameters for the specified stream in Aspen Plus.
        :param block_name: The name of the stream to be modified.
        """
        block_name = self.block if block_name is None else block_name
        try:
            # Access the feed stream node
            block_node = self.aspen.Tree.FindNode(f"\\Data\\Blocks\\{block_name}\\Output")
            CD = block_node.Elements("COND_DUTY").Value / 1000  # kW
            HD = block_node.Elements("REB_DUTY").Value / 1000  # kW
            HT = block_node.Elements("BOTTOM_TEMP").Value  # + 273.15  # Reboiler Temp

            return pd.Series([CD, HD, HT], index=['CD', 'HD', 'HT'])
        except Exception as e:
            print(f"Error accessing block parameters: {e}")
            self.aspen.Quit()
            exit(1)

    def access_to_st(self, stream=None, comp_list=None):
        stream = self.light if stream is None else stream
        comp_list = self.apk_comp if comp_list is None else comp_list
        new_comp = []
        for i in comp_list:
            if i == 'carbon monoxide':
                sub = 'CO'
            elif i == 'Methanol':
                sub = 'MEOH'
            elif i == 'water':
                sub = 'H2O'
            else:
                sub = i
            new_comp.append(sub)
        try:
            stream_node = self.aspen.Tree.FindNode(f"\\Data\\Streams\\{stream}\\Output")
            T = stream_node.Elements("STR_MAIN").Elements("TEMP").Elements("MIXED").Value  # + 273.15
            P = stream_node.Elements("STR_MAIN").Elements("PRES").Elements("MIXED").Value
            flow_cond = pd.Series([T, P], index=['T', 'P'])
            Flow = pd.Series(index=new_comp)
            for compound in new_comp:
                Flow[compound] = stream_node.Elements("STR_MAIN").Elements("MOLEFLOW"). \
                    Elements("MIXED").Elements(compound).Value
            res = pd.concat([Flow, flow_cond])
            new_index = []
            for i in res.index.tolist():
                # if i == 'CO':
                #     sub = 'carbon monoxide'
                if i == 'MEOH':
                    sub = 'Methanol'
                else:
                    sub = i
                new_index.append(sub)
            res.index = new_index
            return res
        except Exception as e:
            print(f"Error accessing stream parameters: {e}")
            self.aspen.Quit()
            exit(1)

    def run_apk(self, stream=None, valid_comp=None, quit_aspen=True):
        valid_comp = self.apk_comp if valid_comp is None else valid_comp
        self.set_feed_parameters(stream_source=stream)
        self.aspen.Engine.Run2()
        block_res = self.access_to_rf()
        light_res = self.access_to_st()
        if quit_aspen:
            self.aspen.Quit()
        return {'light': light_res[valid_comp + ['T', 'P']],
                'block': block_res}

    def opt_distiller(self):

        nstages = np.arange(10, 30, 1)
        Eva_min = 10000
        try:
            # determine the stage number, DF through DSTWU
            block_name = 'B1'
            block_node = self.aspen.Tree.FindNode(f"\\Data\\Blocks\\{block_name}\\Input")
            block_out_node = self.aspen.Tree.FindNode(f"\\Data\\Blocks\\{block_name}\\Output")
            block_node.Elements("OPT_NTRR").Value = 'NSTAGE'
            for nstage in nstages:
                block_node.Elements("NSTAGE").Value = nstage
                self.aspen.Engine.Run2()
                RR = self.aspen.Tree.Elements("Data").Elements("Blocks").Elements('B1').Elements("Output").Elements(
                    "ACT_REFLUX").Value  # block_out_node.Elements("ACT_REFLUX").Value
                FS = min(np.ceil(block_out_node.Elements("FEED_LOCATN").Value), nstage - 1)
                Eva = nstage * RR
                DF = block_out_node.Elements("DIST_VS_FEED").Value
                time.sleep(1)
                if Eva < Eva_min:
                    res = [nstage, RR, FS, DF]
                    Eva_min = Eva
            res = pd.Series(res, index=['nstage', 'RR', 'FS', 'DF'])
            self.set_block_parameters(block_para=res)

            # determine the feed stage
            block_name = 'B3'
            block_node = self.aspen.Tree.FindNode(f"\\Data\\Blocks\\{block_name}\\Input")
            block_out_node = self.aspen.Tree.FindNode(f"\\Data\\Blocks\\{block_name}\\Output")

            FSs = np.arange(np.round(res['nstage'] / 2), res['nstage'] - 1, 1)
            feed_name_rf = 'S6'
            duty = 1e9
            FS_opt = FS
            for FS in FSs:
                block_node.Elements("FEED_STAGE").Elements(feed_name_rf).Value = FS
                self.aspen.Engine.Run2()
                rf_res = self.access_to_rf()
                # rf_rr = self.aspen.Tree.Data.Blocks.B3.Output.RR
                rf_rr = block_out_node.Elements('RR').Value
                hd = rf_res['HD']
                if hd < duty:
                    FS_opt = FS
                    duty = hd
            res['FS'] = FS_opt
            res['RR'] = rf_rr
            return res
        except Exception as e:
            print(f"Error finding best para for distiller: {e}")
            self.aspen.Quit()
            exit(1)

    def run_RF(self, stream=None, valid_comp=None, opt=True, quit_aspen=True, ht_stream=None):
        valid_comp = self.apk_comp if valid_comp is None else valid_comp
        self.set_feed_parameters(stream_source=stream)
        if opt:
            opt_para = self.opt_distiller()
            self.set_block_parameters(block_para=opt_para)
        if ht_stream is not None:
            self.set_feed_parameters(stream_name='S11', stream_source=ht_stream)
        self.aspen.Engine.Run2()
        block_res = self.access_to_rf()
        light_res = self.access_to_st()
        heavy_res = self.access_to_st(stream='S8')
        ht_out = self.access_to_st(stream='S14')
        if quit_aspen:
            self.aspen.Quit()
        if opt:
            return {'light': light_res[valid_comp + ['T', 'P']],
                    'heavy': heavy_res[valid_comp + ['T', 'P']],
                    'block': block_res,
                    'RF': opt_para,
                    'ht': ht_out[valid_comp + ['T', 'P']]}
        else:
            return {'light': light_res[valid_comp + ['T', 'P']],
                    'heavy': heavy_res[valid_comp + ['T', 'P']],
                    'block': block_res,
                    'ht': ht_out[valid_comp + ['T', 'P']]}
    def close_aspen(self):
        if self.aspen is not None:
            try:
                self.aspen.Quit()
            except Exception as e:
                print(f"Warning while quitting Aspen: {e}")
            finally:
                self.aspen = None

class DistillerOpt_new:
    def __init__(self, apk_name, block_name, feed_name, light_name, heavy_name):
        self.block = block_name
        self.feed = feed_name
        self.light = light_name
        self.heavy = heavy_name
        self._init_apk(apk_name)
        self.apk_comp = self._get_comp()

    def _init_apk(self, apk):
        try:
            self.aspen = win32com.client.gencache.EnsureDispatch("Apwn.Document")
            self.aspen.InitFromArchive2(apk)
            self.aspen.Visible = False  # 若需要显示Aspen界面，可改为True
        except Exception as e:
            print(f"An error occurred: {e}")
            self.aspen.Quit()
            exit(1)

    def _get_comp(self):
        try:
            stream_node = self.aspen.Tree.FindNode(f"\\Data\\Streams\\{self.feed}\\Output")
            CompoundLister = stream_node.Elements("STR_MAIN").Elements("MASSFLOW").Elements("MIXED").Elements
            CompoundNameList = []
            for compound in CompoundLister:
                CompoundNameList.append(compound.Name)
            return CompoundNameList
        except Exception as e:
            print(f"An error occurred during get comps: {e}")
            self.aspen.Quit()
            exit(1)

    def set_feed_parameters(self, stream_name=None, stream_source=None):
        stream_name = self.feed if stream_name is None else stream_name
        if stream_source is None:
            raise ValueError("stream_source cannot be None")
        try:
            stream = stream_source.copy()
            new_index = []
            # 这里将一些常见别名统一，确保写入Aspen时名称符合要求
            for i in stream.index.tolist():
                if i.lower() in ['carbon monoxide', 'co']:
                    sub = 'CO'
                elif i.lower() in ['methanol', 'meoh']:
                    sub = 'MEOH'
                elif i.lower() in ['water', 'h2o']:
                    sub = 'H2O'
                else:
                    sub = i
                new_index.append(sub)
            stream.index = new_index
            stream_node = self.aspen.Tree.FindNode(f"\\Data\\Streams\\{stream_name}\\Input")
            stream_node.Elements("TEMP").Elements("MIXED").Value = stream['T']
            stream_node.Elements("PRES").Elements("MIXED").Value = stream['P']
            for comp in self.apk_comp:
                comp_node = stream_node.Elements("FLOW").Elements("MIXED").Elements(comp)
                try:
                    comp_node.Value = stream[comp]
                except KeyError:
                    comp_node.Value = 0.0
        except Exception as e:
            print(f"Error setting feed parameters for stream {stream_name}: {e}")
            self.aspen.Quit()
            exit(1)

    def set_block_parameters(self, block_para=None, block_name=None):
        block_name = self.block if block_name is None else block_name
        block_para = self.opt_distiller() if block_para is None else block_para
        feed_name_rf = 'S6'
        try:
            stream_node = self.aspen.Tree.FindNode(f"\\Data\\Blocks\\{block_name}\\Input")
            stream_node.Elements("NSTAGE").Value = block_para['nstage']
            stream_node.Elements("BASIS_RR").Value = block_para['RR']
            stream_node.Elements("FEED_STAGE").Elements(feed_name_rf).Value = block_para['FS']
            stream_node.Elements("D:F").Value = block_para['DF']
        except Exception as e:
            print(f"Error setting block parameters: {e}")
            self.aspen.Quit()
            exit(1)

    def access_to_rf(self, block_name=None):
        block_name = self.block if block_name is None else block_name
        try:
            block_node = self.aspen.Tree.FindNode(f"\\Data\\Blocks\\{block_name}\\Output")
            CD = block_node.Elements("COND_DUTY").Value / 1000.0  # kW
            HD = block_node.Elements("REB_DUTY").Value / 1000.0     # kW
            HT = block_node.Elements("BOTTOM_TEMP").Value
            return pd.Series([CD, HD, HT], index=['CD', 'HD', 'HT'])
        except Exception as e:
            print(f"Error accessing block parameters: {e}")
            self.aspen.Quit()
            exit(1)

    def access_to_st(self, stream=None, comp_list=None):
        stream = self.light if stream is None else stream
        comp_list = self.apk_comp if comp_list is None else comp_list
        new_comp = []
        for i in comp_list:
            # 统一转换（这里和set_feed_parameters相反，可根据需要修改）
            if i.lower() in ['carbon monoxide', 'co']:
                sub = 'CO'
            elif i.lower() in ['methanol', 'meoh']:
                sub = 'MEOH'
            elif i.lower() in ['water', 'h2o']:
                sub = 'H2O'
            else:
                sub = i
            new_comp.append(sub)
        try:
            stream_node = self.aspen.Tree.FindNode(f"\\Data\\Streams\\{stream}\\Output")
            T = stream_node.Elements("STR_MAIN").Elements("TEMP").Elements("MIXED").Value
            P = stream_node.Elements("STR_MAIN").Elements("PRES").Elements("MIXED").Value
            flow_cond = pd.Series([T, P], index=['T', 'P'])
            Flow = pd.Series(index=new_comp, dtype=float)
            for compound in new_comp:
                Flow[compound] = stream_node.Elements("STR_MAIN").Elements("MOLEFLOW").Elements("MIXED").Elements(compound).Value
            res = pd.concat([Flow, flow_cond])
            # 将“MEOH”转换回“Methanol”，如果你希望外部看到 Methanol
            rename_map = {'MEOH': 'Methanol', 'H2O': 'water', 'CO': 'carbon monoxide'}
            res.index = [rename_map.get(i, i) for i in res.index]
            return res
        except Exception as e:
            print(f"Error accessing stream parameters: {e}")
            self.aspen.Quit()
            exit(1)

    def run_apk(self, stream=None, valid_comp=None, quit_aspen=False):
        if stream is None:
            raise ValueError("stream cannot be None for run_apk.")
        valid_comp = self.apk_comp if valid_comp is None else valid_comp
        # 直接使用 valid_comp，不进行额外映射
        self.set_feed_parameters(stream_source=stream)
        self.aspen.Engine.Run2()
        block_res = self.access_to_rf()
        light_res = self.access_to_st()
        if quit_aspen:
            self.close_aspen()
        return {'light': light_res[valid_comp + ['T', 'P']], 'block': block_res}

    def opt_distiller(self):
        nstages = np.arange(10, 30, 1)
        Eva_min = 10000
        try:
            block_name = 'B1'
            block_node = self.aspen.Tree.FindNode(f"\\Data\\Blocks\\{block_name}\\Input")
            block_out_node = self.aspen.Tree.FindNode(f"\\Data\\Blocks\\{block_name}\\Output")
            block_node.Elements("OPT_NTRR").Value = 'NSTAGE'
            for nstage in nstages:
                block_node.Elements("NSTAGE").Value = nstage
                self.aspen.Engine.Run2()
                RR = self.aspen.Tree.Elements("Data").Elements("Blocks").Elements('B1').Elements("Output").Elements("ACT_REFLUX").Value
                FS = min(np.ceil(block_out_node.Elements("FEED_LOCATN").Value), nstage - 1)
                Eva = nstage * RR
                DF = block_out_node.Elements("DIST_VS_FEED").Value
                time.sleep(1)
                if Eva < Eva_min:
                    res = [nstage, RR, FS, DF]
                    Eva_min = Eva
            res = pd.Series(res, index=['nstage', 'RR', 'FS', 'DF'])
            self.set_block_parameters(block_para=res)
            block_name = 'B3'
            block_node = self.aspen.Tree.FindNode(f"\\Data\\Blocks\\{block_name}\\Input")
            block_out_node = self.aspen.Tree.FindNode(f"\\Data\\Blocks\\{block_name}\\Output")
            FSs = np.arange(np.round(res['nstage'] / 2), res['nstage'] - 1, 1)
            feed_name_rf = 'S6'
            duty = 1e9
            FS_opt = FS
            for FS in FSs:
                block_node.Elements("FEED_STAGE").Elements(feed_name_rf).Value = FS
                self.aspen.Engine.Run2()
                rf_res = self.access_to_rf()
                rf_rr = block_out_node.Elements('RR').Value
                hd = rf_res['HD']
                if hd < duty:
                    FS_opt = FS
                    duty = hd
            res['FS'] = FS_opt
            res['RR'] = rf_rr
            return res
        except Exception as e:
            print(f"Error finding best para for distiller: {e}")
            self.aspen.Quit()
            exit(1)

    def run_RF(self, stream=None, valid_comp=None, opt=True, quit_aspen=True, ht_stream=None):
        #valid_comp = self.apk_comp if valid_comp is None else valid_comp
        valid_comp = ['CO2', 'H2', 'Methanol', 'water', 'carbon monoxide', 'N2']
        # 直接使用 valid_comp，无需额外映射
        self.set_feed_parameters(stream_source=stream)
        if opt:
            opt_para = self.opt_distiller()
            self.set_block_parameters(block_para=opt_para)
        if ht_stream is not None:
            self.set_feed_parameters(stream_name='S11', stream_source=ht_stream)
        self.aspen.Engine.Run2()
        block_res = self.access_to_rf()
        light_res = self.access_to_st()
        heavy_res = self.access_to_st(stream=self.heavy)
        ht_out = self.access_to_st(stream='S14')
        if quit_aspen:
            self.close_aspen()
        if opt:
            return {
                'light': light_res[valid_comp + ['T', 'P']],
                'heavy': heavy_res[valid_comp + ['T', 'P']],
                'block': block_res,
                'RF': opt_para,
                'ht': ht_out[valid_comp + ['T', 'P']]
            }
        else:
            return {
                'light': light_res[valid_comp + ['T', 'P']],
                'heavy': heavy_res[valid_comp + ['T', 'P']],
                'block': block_res,
                'ht': ht_out[valid_comp + ['T', 'P']]
            }

    def close_aspen(self):
        if self.aspen is not None:
            try:
                self.aspen.Quit()
            except Exception as e:
                print(f"Warning while quitting Aspen: {e}")
            finally:
                self.aspen = None


class Distiller:
    def __init__(self, apk, block, feed, light, heavy):

        self.block = block
        self.feed = feed
        self.light = light
        self.heavy = heavy
        self._init_apk(apk)

        self.apk_comp = self._get_comp()

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

    def _get_comp(self):
        try:
            stream_node = self.aspen.Tree.FindNode(f"\\Data\\Streams\\{self.feed}\\Output")
            CompoundLister = stream_node.Elements("STR_MAIN").Elements("MASSFLOW").Elements("MIXED").Elements
            CompoundNameList = []
            for compound in CompoundLister:
                Compoundname = compound.Name
                CompoundNameList.append(Compoundname)
            return CompoundNameList
        except Exception as e:
            print(f"An error occurred during get comps: {e}")
            self.aspen.Quit()
            exit(1)

    def set_feed_parameters(self, stream_name=None, stream_source=None):
        """
        Set the feed parameters for the specified stream in Aspen Plus.

        :param stream_name: The name of the stream to be modified.
        :param stream_source: The new temperature (in K) for the stream, if provided.
        """
        stream_name = self.feed if stream_name is None else stream_name

        try:
            # Access the feed stream node
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
            stream_node = self.aspen.Tree.FindNode(f"\\Data\\Streams\\{stream_name}\\Input")
            stream_node.Elements("TEMP").Elements("MIXED").Value = stream['T']
            stream_node.Elements("PRES").Elements("MIXED").Value = stream['P']
            for comp in self.apk_comp:

                comp_node = stream_node.Elements("FLOW").Elements("MIXED").Elements(comp)
                try:
                    comp_node.Value = stream[comp]
                except KeyError:
                    comp_node.Value = 0

        except Exception as e:
            print(f"Error setting feed parameters: {e}")
            self.aspen.Quit()
            exit(1)

    def access_to_rf(self, block_name=None):
        """
        Set the feed parameters for the specified stream in Aspen Plus.
        :param block_name: The name of the stream to be modified.
        """
        block_name = self.block if block_name is None else block_name
        try:
            # Access the feed stream node
            block_node = self.aspen.Tree.FindNode(f"\\Data\\Blocks\\{block_name}\\Output")
            CD = block_node.Elements("COND_DUTY").Value / 1000  # kW
            HD = block_node.Elements("REB_DUTY").Value / 1000  # kW
            HT = block_node.Elements("BOTTOM_TEMP").Value  # Reboiler Temp

            return pd.Series([CD, HD, HT], index=['CD', 'HD', 'HT'])
        except Exception as e:
            print(f"Error accessing block parameters: {e}")
            self.aspen.Quit()
            exit(1)

    def access_to_st(self, stream=None, comp_list=None):
        stream = self.light if stream is None else stream
        comp_list = self.apk_comp if comp_list is None else comp_list
        new_comp = []
        for i in comp_list:
            if i == 'carbon monoxide':
                sub = 'CO'
            elif i == 'Methanol':
                sub = 'MEOH'
            elif i == 'water':
                sub = 'H2O'
            else:
                sub = i
            new_comp.append(sub)
        try:
            stream_node = self.aspen.Tree.FindNode(f"\\Data\\Streams\\{stream}\\Output")
            T = stream_node.Elements("STR_MAIN").Elements("TEMP").Elements("MIXED").Value
            P = stream_node.Elements("STR_MAIN").Elements("PRES").Elements("MIXED").Value
            flow_cond = pd.Series([T, P], index=['T', 'P'])
            Flow = pd.Series(index=new_comp)
            for compound in new_comp:
                Flow[compound] = stream_node.Elements("STR_MAIN").Elements("MOLEFLOW"). \
                    Elements("MIXED").Elements(compound).Value
            res = pd.concat([Flow, flow_cond])
            new_index = []
            for i in res.index.tolist():
                # if i == 'CO':
                #     sub = 'carbon monoxide'
                if i == 'MEOH':
                    sub = 'Methanol'
                else:
                    sub = i
                new_index.append(sub)
            res.index = new_index
            return res
        except Exception as e:
            print(f"Error accessing stream parameters: {e}")
            self.aspen.Quit()
            exit(1)

    def opt_distiller(self):
        opt_block = r"D:\study\00课题\06多联产系统\甲醇单元\反应器比较\膜反应器\DT_with_MR_opt.bkp"
        block_name = 'B1'
        block_node = self.aspen.Tree.FindNode(f"\\Data\\Blocks\\{block_name}\\Input")
        block_node.Elements("OPT_NTRR").Value = 'NSTAGE'
        # nstages = np.arange(9)
        block_node.Elements("NSTAGE").Value = -1.5

    def run_apk(self, stream=None, valid_comp=None):
        valid_comp = self.apk_comp if valid_comp is None else valid_comp
        self.set_feed_parameters(stream_source=stream)
        self.aspen.Engine.Run2()
        block_res = self.access_to_rf()
        light_res = self.access_to_st()
        self.aspen.Quit()
        return {'light': light_res[valid_comp + ['T', 'P']],
                'block': block_res}


def heater(fluid, T_out):
    """
    :param fluid: fluid infor, including molar flow rate, temperature, pressure
    :param T_out:
    :return:
    """
    # read infor
    T_in = fluid["T"]
    P_in = fluid["P"]
    # Create VLEThermo instance for the fluid
    thermo = VLEThermo(fluid.index.tolist()[:-2])

    # Calculate enthalpy at inlet and outlet temperatures
    h_in = thermo.cal_H(T_in, P_in, fluid.values[:-2])
    h_out = thermo.cal_H(T_out, P_in, fluid.values[:-2])

    # Calculate heat duty (Q)
    Q = (h_out - h_in)  # kW

    # Calculate exergy duty
    e_in = thermo.cal_E(T_in, P_in, fluid.values[:-2])
    e_out = thermo.cal_E(T_out, P_in, fluid.values[:-2])
    E = (e_out - e_in)

    # Returning the results
    fo_cond = pd.Series([T_out, P_in], index=['T', "P"])
    return {
        # 'fluid': pd.concat([fluid, f_cond]),
        'Fin': fluid,
        'Fo': pd.concat([fluid.iloc[:-2], fo_cond]),
        'Q': Q,
        "E": E
    }


def heater_duty(fluid, duty):
    """
    :param fluid: fluid infor, including molar flow rate, temperature, pressure
    :return:
    """
    # read infor
    T_in = fluid["T"]
    P_in = fluid["P"]
    # Create VLEThermo instance for the fluid
    thermo = VLEThermo(fluid.index.tolist()[:-2])

    # Calculate enthalpy at inlet and outlet temperatures
    h_in = thermo.cal_H(T_in, P_in, fluid.values[:-2])
    h_out = h_in + duty
    T_out = thermo.cal_T_from_H(h_out, P_in, fluid.values[:-2])

    # Calculate heat duty (Q)
    Q = (h_out - h_in)  # kW

    # Calculate exergy duty
    e_in = thermo.cal_E(T_in, P_in, fluid.values[:-2])
    e_out = thermo.cal_E(T_out, P_in, fluid.values[:-2])
    E = (e_out - e_in)

    # Returning the results
    fo_cond = pd.Series([T_out, P_in], index=['T', "P"])
    return {
        # 'fluid': pd.concat([fluid, f_cond]),
        'Fin': fluid,
        'Fo': pd.concat([fluid.iloc[:-2], fo_cond]),
        'Q': Q,
        "E": E
    }


def compressor(fluid, compression_ratio, eta_isentropic=0.78, eta_mechanical=0.98):
    """
    :param fluid: fluid info,
    :param compression_ratio:
    :param eta_isentropic:
    :param eta_mechanical:
    :return:
    """
    thermo = VLEThermo(fluid.index.tolist()[:-2])
    # Calculate outlet pressure using compression ratio
    T1, P1 = fluid["T"], fluid["P"]
    P2 = P1 * compression_ratio

    # Get inlet enthalpy and entropy
    h1 = thermo.cal_H(T1, P1, fluid.values[:-2])
    s1 = thermo.cal_S(T1, P1, fluid.values[:-2])
    # Get isentropic outlet temperature and enthalpy
    precision_start = 11
    for i in range(14):
        precision = precision_start + i
        try:
            T2s = thermo.cal_T_from_S(np.round(s1, precision), P2, fluid.values[:-2])
            break
        except TypeError:
            continue
    h2s = thermo.cal_H(T2s, P2, fluid.values[:-2])

    # Isentropic work
    work_isentropic = h2s - h1

    # Actual work and actual outlet enthalpy
    work_actual = work_isentropic / eta_isentropic
    h2_actual = h1 + work_actual

    # Get actual outlet temperature
    T2_actual = thermo.cal_T_from_H(h2_actual, P2, fluid.values[:-2])

    # Work input to the compressor shaft considering mechanical efficiency
    work_input_shaft = work_actual / eta_mechanical

    # Returning the results
    Fo_cond = pd.Series([T2_actual, np.round(P2, 2)], index=['T', "P"])

    return {
        "Fin": fluid,
        "Fo": pd.concat([fluid.iloc[:-2], Fo_cond]),
        'W': work_input_shaft  # kW
    }


def mult_compressor(fluid, P2, T_cool=308.15, r_max=3.5):
    """
    multi-stage compressor with interstage cooling
    :param fluid: fluid info, pd.Series
    :param P2: target pressure
    :param T_cool: interstage cooling temperature
    :param r_max: compression ratio of single compressor
    :return:
    """
    P1 = fluid["P"]

    wf = fluid.copy()
    r_target = P2 / P1
    n_stage = np.ceil(np.log(r_target) / np.log(r_max))
    res = pd.DataFrame(index=np.arange(n_stage) + 1, columns=fluid.index.tolist() + ["W", "Q"], dtype='float64')
    n = 0
    wf["T"] = T_cool
    while n < n_stage:
        r_comp = r_max if n < n_stage - 1 else P2 / wf["P"]
        temp = compressor(wf, r_comp)
        res.iloc[n, :-2] = temp['Fo']
        res.iloc[n, -2] = temp['W']
        res.iloc[n, -1] = heater(temp['Fo'], T_cool)["Q"]
        wf["T"] = T_cool
        wf["P"] = temp['Fo']['P']
        n += 1
    return res


def multi_comp_opt(fluid, P2, T_cool=308.15, r_max=3.5, r_min=2.5):
    P1 = fluid["P"]
    r_target = P2 / P1
    T_tar = 150 + 273.15
    duty_min = 1e10
    for r_comp in np.arange(r_min, r_max, 0.1):
        r_comp=np.round(r_comp,1)
        n_stage = np.ceil(np.log(r_target) / np.log(r_comp))
        res = pd.DataFrame(index=np.arange(n_stage) + 1, columns=fluid.index.tolist() + ["W", "Q"], dtype='float64')
        wf = fluid.copy()
        # wf["T"] = T_cool
        n = 0
        while n < n_stage:
            r_comp = r_comp if n < n_stage - 1 else P2 / wf["P"]
            # res.iloc[n, -1] = heater(temp['Fo'], T_cool)["Q"]
            res.iloc[n, -1] = heater(wf, T_cool)["Q"]
            wf["T"] = T_cool
            temp = compressor(wf, r_comp)
            res.iloc[n, :-2] = temp['Fo']
            res.iloc[n, -2] = temp['W']
            wf["T"] = temp['Fo']['T']
            wf["P"] = temp['Fo']['P']

            n += 1
        heat_res = heater(res.iloc[-1, :-2], T_out=T_tar)
        duty = res['W'].sum() + heat_res['E']
        if duty < duty_min:
            res_opt = res
            duty_min = duty
        if r_comp == r_min:
            res_opt = res
        if res_opt.iloc[-1, -4] > T_tar:
            break
    return res_opt


def spliter(fluid, split_ratio=0.01):
    """
    Split a stream into two parts based on a given split ratio.
    :param fluid: fluid info, pd.Series
    :param split_ratio: The split ratio (e.g., 0.01 for 1% purge)
    :return: Two streams (purge and recycle)
    """
    purge_stream, recycle_stream = fluid.copy(), fluid.copy()
    purge_stream.iloc[:-2] = purge_stream.iloc[:-2] * split_ratio
    recycle_stream.iloc[:-2] = recycle_stream.iloc[:-2] * (1 - split_ratio)
    return {
        'feed': fluid,
        'purge': purge_stream,
        'recycle': recycle_stream
    }


def adsorber(stream, sp_ratio=1 - 1e-4):
    """
    adsorber to drop the water in the stream
    :param stream: fluid info, pd.Series
    :param sp_ratio: The separation ratio (e.g., 0.99 for 99% water adsorbed)
    :return: Two streams (purge and recycle)
    """
    processed_steam = stream.copy()
    processed_steam["H2O"] = processed_steam["H2O"] * (1 - sp_ratio)
    return {
        'Fin': stream,
        'Fo': processed_steam
    }


def flasher(feed, T_cool=None, P_cool=None):
    """
    Cool permeate stream and separate methanol and water from it.
    :param feed: permeate stream (molar flow rates of components), T, P; pd.Series
    :param T_cool: The cooling temperature
    :param P_cool: The cooling pressure
    :return: Gas and liquid streams after cooling
    """
    # Create a VLEThermo instance for permeate stream
    thermo = VLEThermo(feed.index.tolist()[:-2])
    if T_cool is None:
        T_cool = feed["T"]
    if P_cool is None:
        P_cool = feed["P"]

    # Perform flash calculation at the cooling conditions
    Ft = np.sum(feed.values[:-2])
    gas, liquid, sf = thermo.flash(T=T_cool, P=P_cool, x=feed.values[:-2])
    Fl_o = pd.Series(np.array(liquid) * Ft * (1 - sf), index=feed.index.tolist()[:-2])
    Fg_o = pd.Series(np.array(gas) * Ft * sf, index=feed.index.tolist()[:-2])

    # Returning the results
    fo_cond = pd.Series([T_cool, P_cool], index=['T', "P"])
    return {
        'F_in': feed,
        'Fl_o': pd.concat([Fl_o, fo_cond]),
        'Fg_o': pd.concat([Fg_o, fo_cond])
    }

    # return np.array(gas) * Ft * sf, np.array(liquid) * Ft * (1 - sf)


def mixer(F1, F2):
    """
    real mixer
    ref: Modelling, Estimation and Optimization of the Methanol Synthesis with Catalyst Deactivation
    :param F1: fluid1 info; pd.Series
    :param F2: fluid1 info; pd.Series
    :return: molar flux of components, temperature
    """
    if F1['P'] != F2['P']:
        print(f"pressure not equal, {F1['P']} vs {F2['P']}")
        exit(0)
    species = F1.index.tolist()[:-2]
    cal = VLEThermo(species)
    P = F1["P"]
    T1, T2 = F1["T"], F2["T"]
    F_out = F1.copy()
    F_out[species] = F1[species] + F2[species]
    if abs(T1 - T2) > 0.1:
        H_in = cal.cal_H(T1, P, F1.iloc[:-2].values) + cal.cal_H(T2, P, F2.iloc[:-2].values)
        T_out = cal.cal_T_from_H(H_in, P, F_out[species].values)
    else:
        T_out = T1
    # Returning the results
    F_out['T'], F_out['P'] = T_out, P
    return {
        'Fin_1': F1,
        'Fin_2': F2,
        'Fo': F_out
    }


def valve(stream, P=1):
    """
    valve to tune the stream
    :param stream: fluid info, pd.Series
    :param P: the output pressure
    :return: Two streams (purge and recycle)
    """
    species = stream.index.tolist()[:-2]
    T1, P1 = stream['T'], stream['P']
    thermo = VLEThermo(stream.index.tolist()[:-2])
    h1 = thermo.cal_H(T1, P1, stream[species].values)
    T2 = thermo.cal_T_from_H(h1, P, stream[species].values)
    processed_steam = stream.copy()
    processed_steam['T'], processed_steam['P'] = T2, P
    return {
        'Fin': stream,
        'Fo': processed_steam
    }
