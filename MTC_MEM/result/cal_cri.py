import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def get_eta_mod_heat_range_via_interpolation(eta_threshold, excel_path, per_S1_value=100, P0_value=50):
    """
    使用线性插值确定 eta_mod_heat 超过指定阈值时的 per_H2O1 范围。

    参数:
    eta_threshold (float): eta_mod_heat 的阈值。
    excel_path (str): 包含数据的 Excel 文件路径。
    per_S1_value (int, optional): 用于过滤数据的 per_S1 值，默认为 100。
    P0_value (int, optional): 用于过滤数据的 P0 值，默认为 50。

    返回:
    tuple: 包含 per_H2O1 的最小和最大值，以及完整的 per_H2O1 和 eta_mod_heat 列表。
    """
    # 读取 Excel 文件
    df = pd.read_excel(excel_path, header=0, index_col=None)

    # 过滤出 per_S1 等于指定值且 P0 等于指定值的数据
    filtered_df = df[(df['per_S1'] == per_S1_value) & (df['P0'] == P0_value)]

    # 提取 per_H2O1 和 eta_mod_heat 列
    per_H2O1 = filtered_df['per_H2O1'].values
    eta_mod_heat = filtered_df['eta_mod_heat'].values

    # 初始化 per_H2O1 范围
    min_per_H2O1, max_per_H2O1 = None, None

    # 逐对检查相邻点，进行线性插值
    for i in range(len(per_H2O1) - 1):
        if eta_mod_heat[i] < eta_threshold and eta_mod_heat[i + 1] > eta_threshold:
            # 线性插值计算出低于阈值和高于阈值之间的 per_H2O1
            interp_per_H2O1 = np.interp(eta_threshold, [eta_mod_heat[i], eta_mod_heat[i + 1]], [per_H2O1[i], per_H2O1[i + 1]])
            if min_per_H2O1 is None:
                min_per_H2O1 = interp_per_H2O1
            max_per_H2O1 = interp_per_H2O1

        elif eta_mod_heat[i] > eta_threshold and eta_mod_heat[i + 1] < eta_threshold:
            # 线性插值计算出高于阈值和低于阈值之间的 per_H2O1
            interp_per_H2O1 = np.interp(eta_threshold, [eta_mod_heat[i], eta_mod_heat[i + 1]], [per_H2O1[i], per_H2O1[i + 1]])
            if min_per_H2O1 is None:
                min_per_H2O1 = interp_per_H2O1
            max_per_H2O1 = interp_per_H2O1

    return per_H2O1, eta_mod_heat, (min_per_H2O1, max_per_H2O1)

# 使用示例
excel_path = r'D:\study\00课题\06多联产系统\甲醇单元\反应器比较\膜反应器\模拟结果\MR\等温\iso_all\iso_all.xlsx'
eta_threshold = 0.8740615487

# 获取数据和符合条件的 per_H2O1 范围
per_H2O1_values, eta_mod_heat_values, per_H2O1_range = get_eta_mod_heat_range_via_interpolation(eta_threshold, excel_path)

# 输出 per_H2O1 范围
print(f"当 eta_mod_heat > {eta_threshold} 时，per_H2O1 的范围为: {per_H2O1_range}")

# 绘制拟合曲线
plt.figure(figsize=(10, 6))
plt.plot(per_H2O1_values, eta_mod_heat_values, 'o-', label='Data Points')
plt.axhline(y=eta_threshold, color='r', linestyle='--', label=f'Threshold = {eta_threshold}')
if per_H2O1_range[0] and per_H2O1_range[1]:
    plt.axvline(x=per_H2O1_range[0], color='g', linestyle='--', label=f'Min per_H2O1 = {per_H2O1_range[0]:.10f}')
    plt.axvline(x=per_H2O1_range[1], color='b', linestyle='--', label=f'Max per_H2O1 = {per_H2O1_range[1]:.10f}')
plt.xlabel('per_H2O1')
plt.ylabel('eta_mod_heat')
plt.title('eta_mod_heat vs. per_H2O1 with Threshold')
plt.legend()
plt.grid(True)
plt.show()
