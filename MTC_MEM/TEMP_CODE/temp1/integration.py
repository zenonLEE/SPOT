import sympy as sp

r = sp.symbols('r')
T = sp.Function('T')
ode = T(r).diff(r) + r * T(r).diff(r, 2)
res = sp.dsolve(ode)
print(res)
con = {T(0.01): 383, T(r).diff(r).subs(r, 0.025): 100}
res = sp.dsolve(ode, ics=con)
print(res)
