import pandas as pd

from prop_calculator import VLEThermo

sub = ['CO2', 'H2', 'Methanol', 'H2O', 'CO']
H2 = 0.02446337
CO2 = 0.00815446/2
CO2_pd = pd.Series([CO2, 0, 0, 0, 0], index=sub)
H2_pd = pd.Series([0, H2, 0, 0, 0], index=sub)

# hs_pd = pd.Series(ls, index=sub + ['T', 'P'])
a = VLEThermo(sub,ref='ev')
H2_E = a.cal_E(T=353.15, P=30, x=H2_pd)
CO2_E = a.cal_E(T=383.15, P=1, x=CO2_pd)
print(H2_E, CO2_E, H2_E + CO2_E)