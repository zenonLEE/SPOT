import os

from prop_calculator import VLEThermo
import pandas as pd

path = r"D:\document\00Study\05多联产系统\甲醇单元\反应器设计\设计核算\recycle_483_70.xlsx"

data = pd.read_excel(path, sheet_name='recycle', index_col=0)
streams = data.columns.tolist()
metric = pd.DataFrame(index=['H', 'S', "G", 'E'], columns=data.columns)
cal = VLEThermo(["CO2", "H2", "Methanol", "H2O", "carbon monoxide"], ref='ev')
for stream in streams:
    [T, P] = data[stream][:2].values
    f = data.loc[['CO2', 'H2', 'MEOH', 'H2O', 'CO'], stream].values
    metric.loc['H', stream] = cal.cal_H(T, P, f)
    metric.loc['S', stream] = cal.cal_S(T, P, f)
    metric.loc['G', stream] = cal.cal_G(T, P, f)
    metric.loc['E', stream] = cal.cal_E(T, P, f)

res_save = pd.concat((data, metric))
file_name = os.path.basename(path)
new_file_name = os.path.splitext(file_name)[0] + '_E' + os.path.splitext(file_name)[1]
new_path = os.path.join(os.path.dirname(path), new_file_name)

res_save.to_excel(new_path)
