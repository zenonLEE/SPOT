import pandas as pd
import numpy as np
from utility import heater, valve, HeatExchanger, DistillerOpt, flasher
from prop_calculator import VLEThermo
from itertools import combinations

pd.set_option('display.max_columns', None)

feed = [0.000115664,
7.17007E-06,
0.007471147,
0.007680317,
8.58907E-07,
0,
308.15,
49.8563963
        ]
sub = ['CO2', 'H2', 'Methanol', 'H2O', 'CO', 'N2']

s3p = [0.040704173,
0.210747806,
0.000891264,
0.000174391,
0.018133671,
0,
531.448,
50
]
CM2 = [0.000000, 0.000000, 0.000000, 0.003000, 0.000000, 0.000000, 540.135718, 48.000000]
s11 = [0.033022004,
0.188105171,
0.008371497,
0.00785656,
0.018335607,
0,
369.218,
49.8563963]
s19 = [0.004077, 0.000297, 0.000000, 0.004338, 0.000000, 0.000000, 476.574369, 1.000000]

feed = pd.Series(feed, index=sub + ['T', 'P'], name='feed')
s3p = pd.Series(s3p, index=sub + ['T', 'P'], name='S3p')
CM2 = pd.Series(CM2, index=sub + ['T', 'P'], name='CM2')
s11 = pd.Series(s11, index=sub + ['T', 'P'], name='S11')
s19 = pd.Series(s19, index=sub + ['T', 'P'], name='S19')
Tr = 513.15
Pr = 50

cond_record = pd.DataFrame(index=sub + ["T", "P"], dtype='float64')
utility_record = pd.DataFrame(index=['Q', 'E', 'W'], dtype='float64')
vle_cal = VLEThermo(sub)


def ht_filter(df, rb_hd):
    """
    根据给定的DataFrame和rb_hd值，筛选出满足条件的组合并按abs(E_sum)排序。

    参数:
    df (pd.DataFrame): 包含元素的DataFrame，索引为元素名称，列为 'Q' 和 'E'。
    rb_hd (float): 用于筛选组合的Q值的绝对值阈值。

    返回:
    list: 按abs(E_sum)排序的组合列表。
    """

    # 初始化一个空的列表来记录满足条件的组合及其详细信息
    res = []

    # 获取所有可能的两两组合
    for comb in combinations(df.index, 2):
        # 计算两行的Q和E的和
        Q_sum = df.loc[comb[0], 'Q'] + df.loc[comb[1], 'Q']
        E_sum = df.loc[comb[0], 'E'] + df.loc[comb[1], 'E']

        # 检查是否满足 abs(Q_sum) > rb_hd
        if abs(Q_sum) > rb_hd:
            # 获取组合中的元素并按abs(E)进行排序
            elem1 = pd.Series(df.loc[comb[0]])
            elem2 = pd.Series(df.loc[comb[1]])

            if abs(elem1['E']) < abs(elem2['E']):
                elem1, elem2 = elem2, elem1

            combination_info = [
                elem1,  # 按照abs(E)较大的元素放在前面
                elem2,  # 较小的元素放在后面
                pd.Series({'Q_sum': Q_sum, 'E_sum': E_sum})  # 相加结果存为 pd.Series
            ]
            res.append(combination_info)

    # 根据 abs(E_sum) 对 res 进行排序
    res.sort(key=lambda x: abs(x[2]['E_sum']))

    return res


def post_process_ht_multi(feed, ht_streams, Tdew):
    # prepare for distiller
    apk_path = r"D:\study\00课题\06多联产系统\甲醇单元\反应器比较\膜反应器\DT_with_MR_opt_boiler.bkp"
    block_name = 'B3'
    feed_name = 'S18'
    heavy_name = 'S8'
    light_name = 'S7'

    # calculate the max heat can be recovered in HT stream
    # the minimum temperature for possible ht stream when heat recovery
    Tcs = pd.Series([Tr, Tr, Tr, Tdew - 10],
                    index=['S3p', 'CM2', 'CM3', 'S19'])

    # the cooled temperature for possible ht stream
    Tcs0 = pd.Series([Tr, Tr, Tr, 308.15],
                     index=['S3p', 'CM2', 'CM3', 'S19'])

    # the name of ht stream after heat recovery
    St_ro = pd.Series(['S3pr', 'CM3', 'CM3r', 'S19r'],
                      index=['S3p', 'CM2', 'CM3', 'S19'])

    # the name of ht stream after cooler
    St_co = pd.Series(['S3', 'CM3', 'CM1', 'S20'],
                      index=['S3p', 'CM2', 'CM3', 'S19'])

    ht_pd = pd.DataFrame(columns=sub + ['T', 'P'])
    ht_info_pd = pd.DataFrame(columns=['Q', 'E'])

    for ht_stream in ht_streams:
        ht_pd.loc[ht_stream.name] = ht_stream
        heater_res = heater(ht_stream, T_out=Tcs[ht_stream.name])
        ht_info_pd.loc[ht_stream.name] = [heater_res['Q'], heater_res['E']]
    ht_info_pd = ht_info_pd.iloc[(ht_info_pd['E'].abs()).argsort()]  # rank from smallest to biggest

    # record the use status of ht
    # 0 means no use; 1 means in rb; 2 means in preheater
    ht_use = pd.Series(0, index=ht_pd.index)

    # calculate the max heat duty in reboiler
    sp = DistillerOpt(apk_path, block_name=block_name, feed_name=feed_name,
                      light_name=light_name, heavy_name=heavy_name)
    sp_in = feed
    sp_res_check = sp.run_RF(stream=sp_in, valid_comp=sub, ht_stream=None)
    rb_hd = sp_res_check['block']['HD']
    ht_recovery = 0  # ht stream recovery failed
    rb_recovery = 0
    sp_res = None

    # first try if the heat duty can be satisfied by single stream
    ht_filters_sin = ht_info_pd[abs(ht_info_pd['Q']) > rb_hd]

    if not ht_filters_sin.empty:
        ht_filters_sin = ht_filters_sin.iloc[(ht_info_pd['E'].abs()).argsort()]
        candidate_num = len(ht_filters_sin.index)
        for i in range(candidate_num):
            # the waste heat in HT stream can supply the separation
            ht_filter_sin = ht_filters_sin.iloc[i]

            sp = DistillerOpt(apk_path, block_name=block_name, feed_name=feed_name,
                              light_name=light_name, heavy_name=heavy_name)

            sp_res_temp = sp.run_RF(stream=sp_in, valid_comp=sub,
                                    ht_stream=ht_pd.loc[ht_filter_sin.name])
            rb_hd_temp = sp_res_temp['block']['HD']

            if abs(ht_filter_sin['Q']) >= rb_hd_temp:
                rb_recovery = 1
                ht_recovery = 1
                sp_res = sp_res_temp
                ht_use[ht_filter_sin.name] = 1
                for st in ht_pd.index.tolist():
                    if st != ht_filter_sin.name:
                        ht_use[st] = 0
                        cooler = heater(ht_pd.loc[st], T_out=Tcs0[st])
                        cond_record[St_co[st]] = cooler['Fo']
                        utility_record['Q' + st] = [cooler['Q'], cooler['E'], 0]
                    else:
                        cond_record[St_ro[st.name]] = sp_res['ht']
                        cooler = heater(sp_res['ht'], T_out=Tcs[st])
                        cond_record[St_co[st.name]] = cooler['Fo']
                        utility_record['Q' + st] = [cooler['Q'], cooler['E'], 0]
                break

    # if single failed then try to use combination of two stream

    if ht_recovery == 0 and len(ht_streams) > 1:
        ht_combs = ht_filter(ht_info_pd, rb_hd)
        candidate_num = len(ht_combs)
        if candidate_num > 0:
            for i in range(candidate_num):
                ht_comb = ht_combs[i]

                lt_ht_source = ht_pd.loc[ht_comb[1].name]
                pre_heater = HeatExchanger(lt_ht_source, sp_in)
                pre_heater_res = pre_heater.fixed_T_hot(Tcs[lt_ht_source.name])
                sp_in_temp = pre_heater_res['Fc_o']
                if sp_in_temp['T'] > Tdew + 5:
                    pre_heater_res = pre_heater.fixed_T_cold(Tdew + 5)
                    sp_in_temp = pre_heater_res['Fc_o']
                sp = DistillerOpt(apk_path, block_name=block_name, feed_name=feed_name,
                                  light_name=light_name, heavy_name=heavy_name)

                ht_ht_source = ht_pd.loc[ht_comb[0].name]
                sp_res_temp = sp.run_RF(stream=sp_in_temp, valid_comp=sub,
                                        ht_stream=ht_pd.loc[ht_ht_source.name])
                rb_hd_temp = sp_res_temp['block']['HD']
                # verify if the heat in stream satisfy the requirement in RB
                if abs(ht_comb[0]['Q']) >= rb_hd_temp:
                    rb_recovery = 1
                    ht_recovery = 1
                    sp_res = sp_res_temp
                    ht_use[ht_comb[0].name] = 1
                    ht_use[ht_comb[1].name] = 2
                    sp_in = sp_in_temp
                    for st in ht_pd.index.tolist():
                        if st != lt_ht_source.name and st != ht_ht_source.name:
                            ht_use[st] = 0
                            cooler = heater(ht_pd.loc[st], T_out=Tcs0[st])
                            cond_record[St_co[st]] = cooler['Fo']
                            utility_record['Q' + st] = [cooler['Q'], cooler['E'], 0]
                    cond_record[St_ro[ht_ht_source.name]] = sp_res['ht']
                    cond_record[St_ro[lt_ht_source.name]] = pre_heater_res['Fh_o']
                    cooler = heater(sp_res['ht'], T_out=Tcs0[ht_ht_source.name])
                    utility_record['Q' + ht_ht_source.name] = [cooler['Q'], cooler['E'], 0]
                    cooler = heater(pre_heater_res['Fh_o'], T_out=Tcs0[lt_ht_source.name])
                    utility_record['Q' + lt_ht_source.name] = [cooler['Q'], cooler['E'], 0]
                    break

    # the combination of two stream still cannot meet the requirement
    if ht_recovery == 0:
        if len(ht_streams) > 2:
            # the combination of three ht stream is utilized
            lt_ht_source1 = ht_pd.loc[ht_info_pd.iloc[0].name]

            pre_heater1 = HeatExchanger(lt_ht_source1, sp_in)
            pre_heater1_res = pre_heater1.fixed_T_hot(Tout=Tcs[lt_ht_source1.name])
            sp_in_pre1 = pre_heater1_res['Fc_o']

            lt_ht_source2 = ht_pd.loc[ht_info_pd.iloc[1].name]
            pre_heater2 = HeatExchanger(lt_ht_source2, sp_in_pre1)
            pre_heater2_res = pre_heater2.fixed_T_hot(Tout=Tcs[lt_ht_source2.name])
            sp_in_pre2 = pre_heater2_res['Fc_o']

            # distiller feed temperature should not be higher than Tdew+5
            if sp_in_pre2['T'] > Tdew + 5:
                pre_heater2_res = pre_heater2.fixed_T_cold(Tout=Tdew + 5)
                sp_in_pre2 = pre_heater2_res['Fc_o']

            ht_ht_source = ht_pd.loc[ht_info_pd.iloc[2].name]

            sp = DistillerOpt(apk_path, block_name=block_name, feed_name=feed_name,
                              light_name=light_name, heavy_name=heavy_name)

            sp_res_temp = sp.run_RF(stream=sp_in_pre2, valid_comp=sub,
                                    ht_stream=ht_ht_source)
            rb_hd_temp = sp_res_temp['block']['HD']

            # verify if the heat in stream satisfy the requirement in RB
            if abs(ht_info_pd.iloc[2]['Q']) >= rb_hd_temp:
                rb_recovery = 1
                sp_res = sp_res_temp
                ht_use[ht_ht_source.name] = 1
                ht_use[lt_ht_source1.name] = 2
                ht_use[lt_ht_source2.name] = 2
                sp_in = sp_in_pre2

                cond_record[St_ro[lt_ht_source1.name]] = pre_heater1_res['Fh_o']
                cond_record[St_co[lt_ht_source1.name]] = pre_heater1_res['Fh_o']

                cond_record[St_ro[lt_ht_source2.name]] = pre_heater2_res['Fh_o']
                cond_record[St_co[lt_ht_source2.name]] = pre_heater2_res['Fh_o']
                if pre_heater2_res['Fh_o']['T'] > Tcs0[lt_ht_source2.name]:
                    cooler = heater(pre_heater2_res['Fh_o'], T_out=Tcs0[lt_ht_source2.name])
                    utility_record['Q' + lt_ht_source2.name] = [cooler['Q'], cooler['E'], 0]
                    cond_record[St_co[lt_ht_source2.name]] = cooler['Fo']

                cond_record[St_ro[ht_ht_source.name]] = sp_res['ht']
                cooler = heater(sp_res['ht'], T_out=Tcs0[ht_ht_source.name])
                cond_record[St_co[ht_ht_source.name]] = cooler['Fo']
                utility_record['Q' + ht_ht_source.name] = [cooler['Q'], cooler['E'], 0]

                cond_record['S15p2'] = pre_heater1_res['Fc_o']

            else:
                lt_ht_source2 = ht_pd.loc[ht_info_pd.iloc[1].name]
                pre_heater2 = HeatExchanger(lt_ht_source2, sp_in_pre1)
                pre_heater2_res = pre_heater2.fixed_T_hot(Tout=Tcs[lt_ht_source2.name])
                sp_in_pre2 = pre_heater2_res['Fc_o']
                # distiller feed temperature should not be higher than Tdew+5
                if sp_in_pre2['T'] > Tdew + 5:
                    pre_heater2_res = pre_heater2.fixed_T_cold(Tout=Tdew + 5)
                    sp_in_pre2 = pre_heater2_res['Fc_o']

                    # the ht stream cannot meet the requirement of reboiler though there is preheater
                    # hence use an external heat source work in reboiler
                    sp = DistillerOpt(apk_path, block_name=block_name, feed_name=feed_name,
                                      light_name=light_name, heavy_name=heavy_name)
                    sp_in = sp_in_pre2
                    sp_res = sp.run_RF(stream=sp_in, valid_comp=sub,
                                       ht_stream=None)
                    cond_record[St_ro[lt_ht_source1.name]] = pre_heater1_res['Fh_o']
                    cond_record[St_co[lt_ht_source1.name]] = pre_heater1_res['Fh_o']

                    cond_record[St_ro[lt_ht_source2.name]] = pre_heater2_res['Fh_o']
                    cond_record[St_co[lt_ht_source2.name]] = pre_heater2_res['Fh_o']
                    if pre_heater2_res['Fh_o']['T'] > Tcs0[lt_ht_source2.name]:
                        cooler = heater(pre_heater2_res['Fh_o'], T_out=Tcs0[lt_ht_source2.name])
                        utility_record['Q' + lt_ht_source2.name] = [cooler['Q'], cooler['E'], 0]
                        cond_record[St_co[lt_ht_source2.name]] = cooler['Fo']
                    ht_use[ht_ht_source.name] = 0
                    ht_use[lt_ht_source1.name] = 2
                    ht_use[lt_ht_source2.name] = 2
                    cond_record['S15p2'] = pre_heater1_res['Fc_o']
                else:
                    # try to introduce an additional preheater
                    pre_heater3_res = heater(sp_in_pre2, T_out=Tdew + 5)
                    sp_in_pre3 = pre_heater3_res['Fo']

                    # the ht stream cannot meet the requirement of reboiler though there is preheater
                    # hence use an external heat source work in reboiler
                    sp = DistillerOpt(apk_path, block_name=block_name, feed_name=feed_name,
                                      light_name=light_name, heavy_name=heavy_name)
                    sp_in = sp_in_pre3
                    sp_res = sp.run_RF(stream=sp_in, valid_comp=sub,
                                       ht_stream=ht_ht_source)

                    rb_hd_temp = sp_res_temp['block']['HD']
                    # verify if the heat in stream satisfy the requirement in RB
                    if abs(ht_info_pd.iloc[2]['Q']) >= rb_hd_temp:
                        cond_record[St_ro[lt_ht_source1.name]] = pre_heater1_res['Fh_o']
                        cond_record[St_co[lt_ht_source1.name]] = pre_heater1_res['Fh_o']
                        if pre_heater1_res['Fh_o']['T'] > Tcs0[lt_ht_source1.name]:
                            cooler = heater(pre_heater1_res['Fh_o'], T_out=Tcs0[lt_ht_source1.name])
                            cond_record[St_co[lt_ht_source1.name]] = cooler['Fo']
                            utility_record['Q' + St_ro[lt_ht_source1.name]] = [cooler['Q'], cooler['E'], 0]

                        cond_record[St_ro[lt_ht_source2.name]] = pre_heater2_res['Fh_o']
                        cond_record[St_co[lt_ht_source2.name]] = pre_heater2_res['Fh_o']

                        if pre_heater2_res['Fh_o']['T'] > Tcs0[lt_ht_source2.name]:
                            cooler = heater(pre_heater2_res['Fh_o'], T_out=Tcs0[lt_ht_source2.name])
                            cond_record[St_co[lt_ht_source2.name]] = cooler['Fo']
                            utility_record['Q' + St_ro[lt_ht_source2.name]] = [cooler['Q'], cooler['E'], 0]

                        utility_record['Q' + 'sp'] = [pre_heater3_res['Q'], pre_heater3_res['E'], 0]

                        cond_record[St_ro[ht_ht_source.name]] = sp_res['ht']
                        cooler = heater(sp_res['ht'], T_out=Tcs0[ht_ht_source.name])
                        cond_record[St_co[ht_ht_source.name]] = cooler['Fo']
                        utility_record['Q' + St_ro[ht_ht_source.name]] = [cooler['Q'], cooler['E'], 0]

                        rb_recovery = 1

                        ht_use[ht_ht_source.name] = 1
                        ht_use[lt_ht_source1.name] = 2
                        ht_use[lt_ht_source2.name] = 2
                        cond_record['S15p2'] = pre_heater1_res['Fc_o']
                        cond_record['S15p3'] = pre_heater2_res['Fc_o']

                    else:
                        pre_heater3 = HeatExchanger(ht_ht_source, sp_in_pre2)
                        pre_heater3_res = pre_heater3.fixed_T_cold(Tout=Tdew + 5)
                        if pre_heater3_res['Fh_o']['T'] < Tcs0[ht_ht_source.name]:
                            pre_heater3_res = pre_heater3.fixed_T_hot(Tout=Tcs[ht_ht_source.name])
                        # the ht stream cannot meet the requirement of reboiler though there is preheater
                        # hence use an external heat source work in reboiler
                        sp = DistillerOpt(apk_path, block_name=block_name, feed_name=feed_name,
                                          light_name=light_name, heavy_name=heavy_name)
                        sp_in = pre_heater3_res['Fc_o']
                        sp_res = sp.run_RF(stream=sp_in, valid_comp=sub,
                                           ht_stream=None)

                        cond_record[St_ro[lt_ht_source1.name]] = pre_heater1_res['Fh_o']
                        cond_record[St_co[lt_ht_source1.name]] = pre_heater1_res['Fh_o']
                        if pre_heater1_res['Fh_o']['T'] > Tcs0[lt_ht_source1.name]:
                            cooler = heater(pre_heater1_res['Fh_o'], T_out=Tcs0[lt_ht_source1.name])
                            cond_record[St_co[lt_ht_source1.name]] = cooler['Fo']
                            utility_record['Q' + St_ro[lt_ht_source1.name]] = [cooler['Q'], cooler['E'], 0]

                        cond_record[St_ro[lt_ht_source2.name]] = pre_heater2_res['Fh_o']
                        cond_record[St_co[lt_ht_source2.name]] = pre_heater2_res['Fh_o']

                        if pre_heater2_res['Fh_o']['T'] > Tcs0[lt_ht_source2.name]:
                            cooler = heater(pre_heater2_res['Fh_o'], T_out=Tcs0[lt_ht_source2.name])
                            cond_record[St_co[lt_ht_source2.name]] = cooler['Fo']
                            utility_record['Q' + St_ro[lt_ht_source2.name]] = [cooler['Q'], cooler['E'], 0]

                        cond_record[St_ro[ht_ht_source.name]] = pre_heater3_res['Fh_o']
                        cooler = heater(sp_res['ht'], T_out=Tcs0[ht_ht_source.name])
                        cond_record[St_co[ht_ht_source.name]] = cooler['Fo']
                        utility_record['Q' + St_ro[ht_ht_source.name]] = [cooler['Q'], cooler['E'], 0]

                        rb_recovery = 0
                        ht_use[ht_ht_source.name] = 2
                        ht_use[lt_ht_source1.name] = 2
                        ht_use[lt_ht_source2.name] = 2
                        cond_record['S15p2'] = pre_heater1_res['Fc_o']
                        cond_record['S15p3'] = pre_heater2_res['Fc_o']


        elif len(ht_streams) > 1:
            # there are only two ht streams. external ht stream is needed
            lt_ht_source1 = ht_pd.loc[ht_info_pd.iloc[0].name]
            pre_heater1 = HeatExchanger(lt_ht_source1, sp_in)
            pre_heater1_res = pre_heater1.fixed_T_hot(Tout=Tcs[lt_ht_source1.name])
            sp_in_pre1 = pre_heater1_res['Fc_o']
            ht_ht_source = ht_pd.loc[ht_info_pd.iloc[1].name]

            if sp_in_pre1['T'] > Tdew + 5:
                # distiller feed temperature should not be higher than Tdew+5
                # the ht stream cannot meet the requirement of reboiler though there is preheater
                # hence use an external heat source work in reboiler
                pre_heater1 = HeatExchanger(lt_ht_source1, sp_in)
                pre_heater1_res = pre_heater1.fixed_T_cold(Tout=Tdew + 5)
                sp_in = pre_heater1_res['Fc_o']
                sp = DistillerOpt(apk_path, block_name=block_name, feed_name=feed_name,
                                  light_name=light_name, heavy_name=heavy_name)

                sp_res = sp.run_RF(stream=sp_in, valid_comp=sub,
                                   ht_stream=None)

                cond_record[St_ro[lt_ht_source1.name]] = pre_heater1_res['Fh_o']
                cond_record[St_co[lt_ht_source1.name]] = pre_heater1_res['Fh_o']
                if pre_heater1_res['Fh_o']['T'] > Tcs0[lt_ht_source1.name]:
                    cooler = heater(pre_heater1_res['Fh_o'], T_out=Tcs0[lt_ht_source1.name])
                    utility_record['Q' + lt_ht_source1.name] = [cooler['Q'], cooler['E'], 0]
                    cond_record[St_co[lt_ht_source1.name]] = cooler['Fo']

                cooler = heater(ht_ht_source, Tcs0[ht_ht_source.name])
                utility_record['Q' + ht_ht_source.name] = [cooler['Q'], cooler['E'], 0]
                cond_record[St_co[ht_ht_source.name]] = cooler['Fo']
                cond_record[St_ro[ht_ht_source.name]] = cooler['Fo']

                ht_use[ht_ht_source.name] = 0
                ht_use[lt_ht_source1.name] = 2

            else:
                # try to introduce another heat source to preheat the distiller feed further
                pre_heater2_res = heater(sp_in_pre1, T_out=Tdew + 5)
                sp_in_temp = pre_heater2_res['Fo']
                sp = DistillerOpt(apk_path, block_name=block_name, feed_name=feed_name,
                                  light_name=light_name, heavy_name=heavy_name)

                sp_res_temp = sp.run_RF(stream=sp_in_temp, valid_comp=sub,
                                        ht_stream=ht_ht_source)
                rb_hd_temp = sp_res_temp['block']['HD']
                # verify if the heat in stream satisfy the requirement in RB
                if abs(ht_info_pd.iloc[1]['Q']) >= rb_hd_temp:
                    # the ht stream can meet the requirement of reboiler with an additional preheater
                    rb_recovery = 1
                    sp_res = sp_res_temp
                    sp_in = sp_in_pre1

                    utility_record['Q' + 'sp'] = [pre_heater2_res['Q'], pre_heater2_res['E'], 0]

                    cond_record[St_ro[lt_ht_source1.name]] = pre_heater1_res['Fh_o']
                    cond_record[St_co[lt_ht_source1.name]] = pre_heater1_res['Fh_o']
                    if pre_heater1_res['Fh_o']['T'] > Tcs0[lt_ht_source1.name]:
                        cooler = heater(pre_heater1_res['Fh_o'], T_out=Tcs0[lt_ht_source1.name])
                        utility_record['Q' + lt_ht_source1.name] = [cooler['Q'], cooler['E'], 0]
                        cond_record[St_co[lt_ht_source1.name]] = cooler['Fo']

                    cond_record[St_ro[ht_ht_source.name]] = sp_res['ht']
                    cooler = heater(sp_res['ht'], Tcs0[ht_ht_source.name])
                    utility_record['Q' + ht_ht_source.name] = [cooler['Q'], cooler['E'], 0]
                    cond_record[St_co[ht_ht_source.name]] = cooler['Fo']

                    ht_use[ht_ht_source.name] = 1
                    ht_use[lt_ht_source1.name] = 2
                    cond_record['S15p2'] = pre_heater1_res['Fc_o']

                else:
                    # still fails with an additional preheater
                    # hence introduce a heat source working in boiler
                    sp = DistillerOpt(apk_path, block_name=block_name, feed_name=feed_name,
                                      light_name=light_name, heavy_name=heavy_name)

                    sp_res = sp.run_RF(stream=sp_in_pre1, valid_comp=sub,
                                       ht_stream=None)

                    cond_record[St_ro[lt_ht_source1.name]] = pre_heater1_res['Fh_o']
                    cond_record[St_co[lt_ht_source1.name]] = pre_heater1_res['Fh_o']
                    utility_record['Q' + lt_ht_source1.name] = [pre_heater2_res['Q'], pre_heater2_res['E'], 0]

                    cooler = heater(ht_ht_source, Tcs0[ht_ht_source.name])
                    utility_record['Q' + ht_ht_source.name] = [cooler['Q'], cooler['E'], 0]
                    cond_record[St_co[ht_ht_source.name]] = cooler['Fo']
                    cond_record[St_ro[ht_ht_source.name]] = cooler['Fo']

                    ht_use[ht_ht_source.name] = 0
                    ht_use[lt_ht_source1.name] = 2

        else:
            lt_ht_source1 = ht_pd.loc[ht_info_pd.iloc[0].name]
            pre_heater_res = heater(sp_in, T_out=Tdew + 5)

            sp_in_temp = pre_heater_res['Fo']
            sp = DistillerOpt(apk_path, block_name=block_name, feed_name=feed_name,
                              light_name=light_name, heavy_name=heavy_name)

            sp_res_temp = sp.run_RF(stream=sp_in_temp, valid_comp=sub, ht_stream=lt_ht_source1)
            rb_hd_temp = sp_res_temp['block']['HD']
            # verify if the heat in stream satisfy the requirement in RB
            print('heat energy in ht')
            print(ht_info_pd)
            print(rb_hd_temp)
            if abs(ht_info_pd.iloc[0]['Q']) >= rb_hd_temp:

                # the ht stream can meet the requirement of reboiler with an additional preheater
                rb_recovery = 1
                sp_res = sp_res_temp
                sp_in = sp_in_temp

                cond_record[St_ro[lt_ht_source1.name]] = sp_res['ht']
                cooler = heater(sp_res['ht'], Tcs0[lt_ht_source1.name])
                cond_record[St_co[lt_ht_source1.name]] = cooler['Fo']
                utility_record['Q' + St_ro[lt_ht_source1.name]] = [cooler['Q'], cooler['E'], 0]
                utility_record['Q' + 'sp'] = [pre_heater_res['Q'], pre_heater_res['E'], 0]

                ht_use[lt_ht_source1.name] = 1
            else:
                pre_heater = HeatExchanger(lt_ht_source1, sp_in)
                pre_heater_res = pre_heater.fixed_T_cold(Tout=Tdew + 5)

                if pre_heater_res['Fh_o']['T'] > Tcs[lt_ht_source1.name]:
                    sp_in = pre_heater_res['Fc_o']
                    sp = DistillerOpt(apk_path, block_name=block_name, feed_name=feed_name,
                                      light_name=light_name, heavy_name=heavy_name)

                    sp_res = sp.run_RF(stream=sp_in, valid_comp=sub, ht_stream=None)
                    cond_record[St_ro[lt_ht_source1.name]] = pre_heater_res['Fh_o']
                    cooler = heater(pre_heater_res['Fh_o'], Tcs0[lt_ht_source1.name])
                    cond_record[St_co[lt_ht_source1.name]] = cooler['Fo']
                    utility_record['Q' + St_ro[lt_ht_source1.name]] = [cooler['Q'], cooler['E'], 0]
                    ht_use[lt_ht_source1.name] = 2
                else:
                    pre_heater_res = pre_heater.fixed_T_hot(Tout=Tcs[lt_ht_source1.name])
                    sp_in = pre_heater_res['Fc_o']
                    sp = DistillerOpt(apk_path, block_name=block_name, feed_name=feed_name,
                                      light_name=light_name, heavy_name=heavy_name)

                    sp_res = sp.run_RF(stream=sp_in, valid_comp=sub, ht_stream=None)
                    cond_record[St_ro[lt_ht_source1.name]] = pre_heater_res['Fh_o']
                    cond_record[St_co[lt_ht_source1.name]] = pre_heater_res['Fh_o']
                    ht_use[lt_ht_source1.name] = 2

    cond_record['S15'] = sp_in

    print(ht_use)
    return rb_recovery, sp_res


def post_process(feed, ht_streams, lt_stream):
    # print(f"feed: {feed}")
    # print(f"ht: {ht_streams}")
    # print(f"lt: {lt_stream}")
    # throttle
    valve_res = valve(feed, P=1.2)
    liq_valve = valve_res['Fo']
    flash_res = flasher(liq_valve)
    flash_liq, flash_gas = flash_res['Fl_o'], flash_res['Fg_o']

    # distill
    # heat liquid before entering column using reacted flow
    lt_hx = HeatExchanger(lt_stream, flash_liq, U=500)
    lt_hx_res = lt_hx.fixed_delta_cold(delta_T=20)
    lt_recovery = lt_hx_res['Fh_o']
    lt_cooler = heater(lt_recovery, T_out=308.15)

    liq_heated_lt = lt_hx_res['Fc_o']
    cond_record['S15p'] = liq_heated_lt.copy()
    sp_in = liq_heated_lt

    # calculate the dew temperature of product
    try:
        T_dew_in = vle_cal.T_dew(flash_liq['P'], flash_liq[sub].values)
    except ValueError:
        T_dew_in = 368.15
    # print(T_dew_in)
    # print(f"feed: {feed}")
    if ht_streams is not None:
        rb_recover, sp_res = post_process_ht_multi(feed=sp_in, ht_streams=ht_streams, Tdew=T_dew_in)

    # else:
    #     sp_res = post_process_no_rc(sp_in, T_dew_in)
    #     ht_recover = 0
    #     save_ht_no_rec()

    product = sp_res['light'].copy()
    heavy_product = sp_res['heavy'].copy()
    reboiler_hd = sp_res['block']['HD']
    reboiler_T = sp_res['block']['HT']
    reboiler_hd_e = reboiler_hd * (1 - 298.15 / reboiler_T) * (1 - rb_recover)

    # record
    cond_record['S12'] = liq_valve.copy()
    cond_record['S13'] = flash_liq
    cond_record['S14'] = flash_gas
    cond_record['S5p'] = lt_recovery

    cond_record['S16'] = heavy_product
    cond_record['S17'] = product

    utility_record['RB'] = [reboiler_hd, reboiler_hd_e, 0]
    utility_record['Q2p'] = [lt_hx_res['Q'], 0, 0]
    utility_record['Q2'] = [lt_cooler['Q'], 0, 0]

    try:
        rf_record = sp_res['RF']
    except KeyError:
        rf_record = None
    return product


res = post_process(feed=feed, ht_streams=[s3p], lt_stream=s11)
print(cond_record)
print(utility_record)
