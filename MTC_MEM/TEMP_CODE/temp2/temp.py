import numpy as np
import matplotlib.pyplot as plt

# 您的数据
# 您的数据
T = [
463.15,
473.15,
483.15,
493.15,
503.15,
513.15,
523.15


]
P = [
30,
30,
30,
30,
30,
30,
40,
50,
40,
50,
40,
50,
40,
50,
40,
50,
40,
50,
40,
50
]  # 由于P的值是重复的，我们将忽略它
r = [
22.27954942,
14.13459643,
8.383439982,
4.446160172,
2.826392061,
2.595202623,
2.895134812



]

# 将列表转换为NumPy数组
T = np.array(T)
r = np.array(r)

# 选择多项式的度数，这里我们选择3次多项式作为示例
polynomial_degree = 3

# 执行多项式拟合
coefficients = np.polyfit(T, r, polynomial_degree)
print(coefficients)
# 使用多项式系数进行预测
T_new = np.linspace(min(T)-10, max(T)+10, 10)  # 创建一个新的温度范围用于绘制拟合曲线
r_fit = np.polyval(coefficients, T_new)

# 绘制原始数据点和拟合曲线
plt.scatter(T, r, label='Data points')
plt.plot(T_new, r_fit, label=f'Polynomial fit of degree {polynomial_degree}', color='red')

# 添加图例和标签
plt.legend()
plt.xlabel('Temperature (T)')
plt.ylabel('r')
plt.title('Polynomial Fit of Temperature vs r')
plt.show()