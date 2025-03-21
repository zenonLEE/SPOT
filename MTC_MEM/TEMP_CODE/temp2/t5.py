import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from itertools import combinations

d = [-1, -0.5]
dn = 'S1'
c = [-1, -0.6]
cn = 'S2'

e = [-0.5, -0.2]
en = 'S3'

rb_hd = 1.2
df_in = pd.DataFrame(columns=['Q','E'])
df_in.loc[dn] = d
df_in.loc[cn] = c
df_in.loc[en] = e

def filter(df, rb_hd):
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


res = filter(df_in, rb_hd)
print(res[0][0])
print(res[0][1])
print(res[0][1].name)
# for i in res:
#     print(i)

# ht_filters = pd_info[abs(pd_info['Q']) > rb_hd]
# print(ht_filters)
# print(ht_filters.empty)
# ht_filter = ht_filters.loc[abs(ht_filters['E']).idxmin()]
# print(ht_filter)
# a.loc['S', b.name] =[0, 0]
# a.loc['Q', b.name] = 0
# a.loc['E', b.name] = 1
#
# a.loc['S', c.name] = c.values
# a.loc['Q', c.name] = 10
# a.loc['E', c.name] = 11


# e_duty_term = a.loc[:, a.loc['Q'] > 0]
# q_exergy_duty = e_duty_term.loc['E'].sum()  # self.utility_record.loc['E', 'RB']
# heat_duty = 0
# for col in e_duty_term.columns:
#     if e_duty_term[col]['E'] >= 0:
#         heat_duty += e_duty_term[col]['Q']
#
# print(q_exergy_duty)
# print(heat_duty)
