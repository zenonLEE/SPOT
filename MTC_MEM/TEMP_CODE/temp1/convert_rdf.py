import numpy as np
import scipy.optimize as opt

file = r"D:\document\00Study\05多联产系统\甲醇单元\反应器设计\物性计算\H2O-CO2\production\rdf_423_70_1.dat"
file_out = r"D:\document\00Study\05多联产系统\甲醇单元\反应器设计\物性计算\H2O-CO2\production\rdf_423_70_1_convert.dat"
with open(file=file, mode="r") as fb:
    data = fb.readlines()
data_save = np.zeros((len(data) - 1, 7))
i = 0
for line in data[1:]:
    r = float(line.strip("\n").split(" ")[0])
    temp = line.strip("\n").split(" ")[1]
    rdf = [float(i) for i in temp.strip('\t').split("\t")]
    rdf = [0 for i in rdf] if r < 0.8 else rdf
    data_save[i, 0] = r
    data_save[i, 1:] = rdf
    i += 1
dr = data_save[2, 0] - data_save[1, 0]
R = data_save[-1, 0]

rdf_corrected = np.zeros((len(data) - 1, 7))
rdf_corrected[:, :4] = data_save[:, [0, 2, 4, 6]]

for i in np.arange(1, 4, 1):
    rdf_corrected[:, i + 3] = dr * 4 * np.pi * (rdf_corrected[:, i] - 1) * \
                              (1 - 3 * rdf_corrected[:, 0] / 2 / R + rdf_corrected[:, 0] ** 3 / 16 / (R / 2) ** 3) * \
                              rdf_corrected[:, 0] ** 2

G = np.sum(rdf_corrected[:, 4:], axis=0)
ct = 1000 / 734241
x1, x2 = 0.9, 0.1
tau = 1 / (1 + ct * x1 * x2 * (G[0] + G[2] - 2 * G[1]))
print(tau)
np.savetxt(file_out, data_save, fmt="%.10e")
