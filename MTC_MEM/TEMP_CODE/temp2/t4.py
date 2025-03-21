import numpy as np
import pandas as pd

from utility import DistillerOpt, HeatExchanger, heater

ht_stream = pd.Series([466.821, 48, 0.3], index=['T', 'P', 'qm'])
lt_stream = [0.007601, 0.047227, 0.004233, 0.001386, 0.002391, 0]

Tlt = 369.546000
Plt = 39.984576

sub = ['CO2', 'H2', 'Methanol', 'H2O', 'CO', 'N2']
lt_pd = pd.Series(lt_stream + [Tlt, Plt], index=sub + ['T', 'P'])
ht_pd = pd.Series(np.zeros(8), index=lt_pd.index)
ht_pd[['T', 'P', 'H2O']] = [500.798527, 48, 0.3]
liq_heated_ult = pd.Series([1.649393e-05, 6.618428e-09, 6.422239e-03, 6.705136e-03, 4.825489e-10, 0, 3.648890e+02, 1.2],
                           index=sub + ['T', 'P'])
liq_heated_lt = pd.Series([1.667324e-05, 1.440638e-08, 3.877844e-03, 1.355306e-03, 1.784643e-09, 0, 3.485160e+02, 1.2],
                          index=sub + ['T', 'P'])

apk_path = r"D:\study\00课题\06多联产系统\甲醇单元\反应器比较\膜反应器\DT_with_MR_opt_boiler.bkp"
block_name = 'B3'
feed_name = 'S18'
heavy_name = 'S8'
light_name = 'S7'
sp = DistillerOpt(apk_path, block_name=block_name, feed_name=feed_name,
                  light_name=light_name, heavy_name=heavy_name)

sp_res = sp.run_RF(stream=liq_heated_ult, valid_comp=sub, ht_stream=ht_pd)
print(sp_res)
