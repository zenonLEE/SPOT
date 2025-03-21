import numpy as np
import pandas as pd

from gibbs import Gibbs
sub = ["CO2", "H2", "Methanol", "H2O", "carbon monoxide"]
a = Gibbs(sub)
feed = pd.Series(np.array([1,3,0,0,0]),index=sub)
T = np.arange(200,290,10)+273.15
p = np.array([20,30,40,50,60,70,80,90,100])
a.solve_over_range(T,p,feed,save=True)