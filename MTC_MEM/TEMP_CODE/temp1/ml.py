import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.svm import SVR
from sklearn.preprocessing import StandardScaler
from sklearn.tree import DecisionTreeRegressor
from sklearn.metrics import r2_score
from sklearn.model_selection import cross_val_score

path = r"D:\document\04Code\PycharmProjects\MTC\result\sim_one_pass_0.2_BU_2_2023-11-22-23_log.xlsx"

data = pd.read_excel(path)
data = data.drop(data.columns[0], axis=1)

X = data.iloc[:, :3].values
Y = data.iloc[:, -1].values
Xs = StandardScaler().fit_transform(X)

Xtrain, Xtest, Ytrain, Ytest = train_test_split(Xs, Y, test_size=0.3, random_state=520)

# dtr = DecisionTreeRegressor(criterion="absolute_error", max_depth=3)
# dtr.fit(Xtrain, Ytrain)
# r2 = dtr.score(Xtest, Ytest)
# print(r2)

# 交叉验证/剪枝参数优化
# cv_list = []
#
# for i in range(1, 10):
#     dtr = DecisionTreeRegressor(criterion="squared_error" ,max_depth=i ,min_samples_leaf=10)
#   #  scoring='neg_mean_absolute_error'  # 交叉验证默认评价参数为 R2，可通过 scoring=“” 修改为其它参数，在这里不作修改
#     cv = cross_val_score(dtr, X, Y, cv=10)
#     cv_list.append(-cv.mean())
#     #print(f"i is {i}, the mean of cv is {cv.mean()}, the score is {dtr.score(Xtest, Ytest)}")

dtr_2 = DecisionTreeRegressor(max_depth=3, random_state=0)
dtr_2.fit(Xtrain, Ytrain)
Ypred = dtr_2.predict(Xtest)

# 绘制最大深度为2时回归树的拟合效果
fig, ax = plt.subplots(ncols=2, nrows=1, sharey=True, figsize=(10, 5))
ax[0].set_xlabel('Velocity')
ax[0].set_ylabel('Capacity')
ax[0].scatter(Xtest[:, 0], Ytest, c='r')
ax[0].scatter(Xtest[:, 0], Ypred)

# 将预测结果降序排序，以便plot
df_results = pd.DataFrame(data={'Xtest': Xtest[:, 0].ravel(),
                                'Ytest': Ytest,
                                'Ypred': Ypred})
df_results.sort_values('Ytest', inplace=True)

ax[1].scatter(df_results['Xtest'], df_results['Ytest'], c='r')
ax[1].plot(df_results['Xtest'], df_results['Ypred'])
plt.show()
