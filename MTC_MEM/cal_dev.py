import pandas as pd
import matplotlib.pyplot as plt
# Load the Excel file
file_path = r'D:\study\00课题\06多联产系统\甲醇单元\反应器比较\膜反应器\模拟结果\MR\等温\50bar_mem/MR_ISO_50_BU_0.0001_2024-08-23_log.xlsx'
data = pd.read_excel(file_path)

# Group by 'per_H2O1' and calculate the variance of 'p_conv'
variance_results = data.groupby('per_H2O1')['p_conv'].var()

# Display the results
print(variance_results)

plt.plot(variance_results)
plt.show()