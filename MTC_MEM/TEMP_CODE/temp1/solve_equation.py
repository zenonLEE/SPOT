import numpy as np
from scipy.optimize import fsolve

T = 513
kp = 0.3 #0.21 + 0.00015 * T
kf = 0.14233
# print(kf/kp*100)
vof = 0.33
Cf = 1.25*1
B = Cf * ((1 - vof) / vof) ** 1.11
kappa = kp / kf


# Define the function representing the system of equations
def equations(vars):
    x = vars
    eq1 = x / kf  # Equation 1
    eq2 = (1 - (1 - vof) ** 0.5) + (2 * (1 - vof) ** 0.5 * x / (1 - B / kappa)) * (
            B * (1 - 1 / kappa) / (1 - B / kappa) ** 2 * np.log(kappa / B) - (B - 1) / (1 - B / kappa) + (
            B + 1) / 2)
    return eq1 - eq2


# Initial guess for the solution
initial_guess = [1]

# Solve the system of equations
solution = fsolve(equations, initial_guess)

# Print the solution
print("Solution:", solution)
