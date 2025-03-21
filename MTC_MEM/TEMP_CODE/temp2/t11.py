import pandas as pd

from utility import multi_comp_opt

f = [8.153639e-06, 1.200393e-03, 0.000000e+00, 5.517946e-09, 0.000000e+00, 0.000000e+00, 3.081500e+02, 1.000000e+00]
sub = ['CO2', 'H2', 'Methanol', 'H2O', 'CO', 'N2', 'T', 'P']
hs_pd = pd.Series(f, index=sub)

res = multi_comp_opt(hs_pd, P2=50)
print(res)
