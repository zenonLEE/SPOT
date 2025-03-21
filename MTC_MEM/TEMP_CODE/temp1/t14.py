import numpy as np
import pandas as pd
from prop_calculator import VLEThermo


# k = np.array([[0, 0, 0],
#               [0, 0, 0],
#               [0, 0, 0]])
# eos = VLEThermo(["O2", "H2", "H2O"], ref='ev')
#
#
# g1 = np.array([0, 0, 1])
# g2 = np.array([0, 1, 0])
# g3 = np.array([0.5, 0, 0])
#
# T = 353.15
# Pp = 1
# a = eos.cal_E(T, Pp, g1)
# b = eos.cal_E(T, Pp, g2)
# c = eos.cal_E(T, Pp, g3)
# print(a, b, c)
# print(b + c - a)

# Hfgs=[-393474.0, 0.0, -200700.0, -241822.0, -110525.0]
# Hvap_298s=[5265.543624363681, 0.0, 37457.19357385385, 43987.4460555689, 0.0]
# Hf_STPs=[-393474.0, 0.0, -238400.0, -285825.0, -110525.0]
# Hvap_298s=[5265.543624363681, 0.0, 37457.19357385385, 43987.4460555689, 0.0]
# print(Hf_STPs[3]-Hfgs[3])
# eos = VLEThermo(np.array(["O2", "H2", "H2O"]), ref='ch', kij=k)
#
# g1 = np.array([0, 0, 1])
# g2 = np.array([0, 1, 0])
# g3 = np.array([0.5, 0, 0])
#
# T = 298.15
# Pp = 1
# a = eos.cal_H(T, Pp, g1)
# b = eos.cal_H(T, Pp, g2)
# c = eos.cal_H(T, Pp, g3)
# print(a, b, c)
# print(b + c - a)
# eos = VLEThermo(np.array(["CO2", "H2", "Methanol", "H2O", "carbon monoxide"]), ref='ev')
#
# g1 = np.array([1, 3, 0, 0, 0])
# g2 = np.array(([0.582, 1.746, 0.418, 0.418, 0]))
#
# T = 473
# Pp = 60
# a = eos.cal_H(T, Pp, g1)
# b = eos.cal_H(T, Pp, g2)
# print((b - a)/0.418)
#
# eos = VLEThermo(np.array(["CO2", "H2", "Methanol", "H2O", "carbon monoxide"]), ref='ch')
#
# g1 = np.array([1, 3, 0, 0, 0])
# g2 = np.array(([0.582, 1.746, 0.418, 0.418, 0]))
#
# T = 473
# Pp = 60
# a = eos.cal_H(T, Pp, g1)
# b = eos.cal_H(T, Pp, g2)
# print((b - a)/0.418)


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

data = np.array([[979.02143, 75.28314, 756.99933, 209.56516, 263.34064],
                 [807.67733, 63.06492, 691.89727, 85.11491, 40.91851],
                 [0.82498, 0.8377, 0.914, 0.40615, 0.15538]])

in_data = data[:, [0, 3]]
out_data = data[:, [1, 2, 4]]
labels = ['In', 'Left', 'Cond', 'Qin', 'Qo']

fig, ax = plt.subplots()

# 绘制 'In' 和 'Qin' 组的矩形
for i, (h_value, alpha_value, area_value) in enumerate(zip(in_data[0], in_data[2], in_data[1])):
    # 计算矩形左下角坐标
    x = sum(in_data[0][:i])  # 累积前面的H值
    y = 0  # 矩形的底部在y=0处

    # 创建矩形并添加到图中
    rect = Rectangle((x, y), h_value, alpha_value, linewidth=2, edgecolor='none', facecolor='blue', alpha=0.05)
    ax.add_patch(rect)

    # 在矩形上方添加边框
    ax.plot([x, x + h_value], [alpha_value, alpha_value], color='blue', linewidth=2)

    # 在矩形上方添加面积值
    ax.text(x + h_value / 2, alpha_value/2, f'{area_value:.2f}', ha='center', va='bottom',color='blue',)

# 绘制 'Left', 'Cond', 'Qo' 组的矩形
for i, (h_value, alpha_value, area_value) in enumerate(zip(out_data[0], out_data[2], out_data[1])):
    # 计算矩形左下角坐标
    x = sum(out_data[0][:i])   # 累积前面的H值，并考虑上一组矩形的宽度
    y = 0  # 矩形的底部在y=0处

    # 创建矩形并添加到图中
    rect = Rectangle((x, y), h_value, alpha_value, linewidth=1, edgecolor='none', facecolor='red', alpha=0.05)
    ax.add_patch(rect)

    # 在矩形上方添加边框
    ax.plot([x, x + h_value], [alpha_value, alpha_value], color='red', linewidth=2)

    # 在矩形上方添加面积值
    ax.text(x + h_value / 2, alpha_value/2, f'{area_value:.2f}', ha='center', va='bottom',color='red')

ax.set_xlim(0, np.sum(in_data[0])*1.05)  # 设置x轴范围
ax.set_ylim(0, 1)  # 设置y轴范围

ax.set_xlabel('H')
ax.set_ylabel('alpha')

plt.show()
