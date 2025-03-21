import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import r2_score
from sklearn.tree import DecisionTreeRegressor

# 生成模拟数据，假设你的数据存储在一个CSV文件中
path = r"D:\document\04Code\PycharmProjects\MTC\result\sim_one_pass_0.2_BU_2_2023-11-22-23_log.xlsx"
keys = ['Dt2', 'Thick2', 'Tc2','heater2','conversion','y_CH3OH','eff']
data = pd.read_excel(path)[keys]
data = data[data['eff']<10]

# 这里用随机数生成模拟数据，你需要替换为你的实际数据
np.random.seed(42)


# 划分数据集为训练集和测试集
X = data[keys[:-3]]
y = data['eff']
Xs = StandardScaler().fit_transform(X)
X_train, X_test, y_train, y_test = train_test_split(Xs, y, test_size=0.2, random_state=42)

# 建立线性回归模型
model = DecisionTreeRegressor(max_depth=4, random_state=0)
model.fit(X_train, y_train)

# 在测试集上进行预测
y_pred = model.predict(X_test)

# 评估模型性能
r2 = r2_score(y_test, y_pred)
print(f'R2: {r2}')

# 输出模型系数
# coefficients = model.coef_
# intercept = model.intercept_
# print(f'Coefficients: {coefficients}')
# print(f'Intercept: {intercept}')
importances = model.feature_importances_
print(f'Feature Importances: {importances}')

# 可视化预测值和实际值的比较
plt.scatter(y_test, y_pred)
plt.xlabel('Actual Performance')
plt.ylabel('Predicted Performance')
plt.title('Actual vs Predicted Performance')
plt.show()
