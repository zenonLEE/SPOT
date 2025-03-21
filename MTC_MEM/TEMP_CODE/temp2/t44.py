import os
import re
import pandas as pd

# 文件夹路径
folder_path = r"D:\study\00课题\06多联产系统\甲醇单元\反应器比较\膜反应器\模拟结果\TR\绝热v3.1"

# 获取当前文件夹的名称
folder_name = os.path.basename(folder_path)

# 获取所有以'sim_point'开头的Excel文件名
files = [f for f in os.listdir(folder_path) if f.startswith("sim_point") and f.endswith(".xlsx")]

# 存储结果的列表
data = []

# 读取每个文件中的'stream' sheet，并提取所需的数据
for file in files:
    file_path = os.path.join(folder_path, file)
    # 读取sheet并指定第一列为index，第一行为列名
    df = pd.read_excel(file_path, sheet_name='stream', index_col=0)

    # 提取列名为'S15'，行名为'T'的数据
    if 'S15' in df.columns and 'T' in df.index:
        ts15_value = df.loc['T', 'S15']
    else:
        continue  # 如果数据不在，则跳过该文件

    # 从文件名中提取温度、压力和参数
    match = re.search(r"(\d+\.\d+)_(\d+\.\d+)_.+_(\d+e?[-+]?\d*)_(\d+\.?\d*)", file, re.IGNORECASE)
    if match:
        temperature = float(match.group(1))
        pressure = float(match.group(2))
        p1 = match.group(3)
        p2 = match.group(4)
        # 检查是否有更多参数（如0.5, 0.0）
        additional_params = file.split("_")[-2:]
    else:
        continue  # 如果文件名不符合预期格式，跳过

    # 将数据存储为一个字典，并添加到列表中
    data.append({
        'T': temperature,
        'P': pressure,
        'p1': p1,
        'p2': p2,
        'TS15': ts15_value,
        'additional_param_1': additional_params[0] if len(additional_params) > 0 else None,
        'additional_param_2': additional_params[1] if len(additional_params) > 1 else None
    })

# 将列表转换为DataFrame
result_df = pd.DataFrame(data)

# 输出结果
print(result_df)

# 构建输出文件名和路径
output_filename = f"{folder_name}_T_S15.csv"
output_path = os.path.join(folder_path, output_filename)

# 将结果保存为CSV文件
result_df.to_csv(output_path, index=False)

print(f"File saved as: {output_path}")
