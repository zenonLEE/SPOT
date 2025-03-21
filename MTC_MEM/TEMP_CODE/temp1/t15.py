import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

data = np.array([[-146.73896, -52.71514, 0.359244334],
                 [-209.56516, -85.11491, 0.406150097],
                 [263.34064, 40.91851, 0.155382435],
                 ])
data = np.abs(data.T)
in_data = data[:, [0, 1]]
out_data = data[:, [2]]
labels = ['r', 'h', 'c']

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
    ax.text(x + h_value / 2, alpha_value / 2, f'{area_value:.2f}', ha='center', va='bottom', color='blue', )

# 绘制 'Left', 'Cond', 'Qo' 组的矩形
for i, (h_value, alpha_value, area_value) in enumerate(zip(out_data[0], out_data[2], out_data[1])):
    # 计算矩形左下角坐标
    x = sum(out_data[0][:i])  # 累积前面的H值，并考虑上一组矩形的宽度
    y = 0  # 矩形的底部在y=0处

    # 创建矩形并添加到图中
    rect = Rectangle((x, y), h_value, alpha_value, linewidth=1, edgecolor='none', facecolor='red', alpha=0.05)
    ax.add_patch(rect)

    # 在矩形上方添加边框
    ax.plot([x, x + h_value], [alpha_value, alpha_value], color='red', linewidth=2)

    # 在矩形上方添加面积值
    ax.text(x + h_value / 2, alpha_value / 2, f'{area_value:.2f}', ha='center', va='bottom', color='red')

ax.set_xlim(0, np.sum(in_data[0]) * 1.05)  # 设置x轴范围
ax.set_ylim(0, 1)  # 设置y轴范围

ax.set_xlabel('dH')
ax.set_ylabel('A')

plt.show()
