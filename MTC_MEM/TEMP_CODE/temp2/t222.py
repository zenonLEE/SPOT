import numpy as np
import pandas as pd
from prop_calculator import VLEThermo
from utility import HX4Water, HeatExchanger, DistillerOpt, heater, heater_duty, flasher

Tlt = 369.546000
Plt = 39.984576
# prepare for distiller
apk_path = r"D:\study\00课题\06多联产系统\甲醇单元\反应器比较\膜反应器\DT_with_MR_opt_boiler.bkp"
block_name = 'B3'
feed_name = 'S18'
heavy_name = 'S8'
light_name = 'S7'

sub = ['CO2', 'H2', 'Methanol', 'H2O', 'CO', 'N2']
feed = pd.Series([0,
1.679659e-05,
3.148009e-07,
8.058405e-03,
1.531932e-05,
1.528122e-07,
4.280060e+02,
1.2
], index=sub + ['T', 'P'])

ht_stream = pd.Series([0,
0,
0,
0.5,
0,
0,
522.39,
48
], index=sub + ['T', 'P'])

lt_stream = pd.Series([0.005806,
0.113472,
0,
0.008096,
0,
0,
405.486940,
1
], index=sub + ['T', 'P'])

# calculate the max heat can be recovered in HT stream
Qr_max = heater(ht_stream, T_out=522.39)['Q']

# calculate the max heat duty in reboiler
sp = DistillerOpt(apk_path, block_name=block_name, feed_name=feed_name,
                  light_name=light_name, heavy_name=heavy_name)
sp_in = feed
sp_res_check = sp.run_RF(stream=sp_in, valid_comp=sub, ht_stream=None)
rb_hd = sp_res_check['block']['HD']
print()