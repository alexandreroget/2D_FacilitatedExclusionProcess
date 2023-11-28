import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize

rho = []
rho_a = []
a = []

with open("data/results.txt", "r") as file:
    for line in file:
        parts = line.strip().split(";")
        if len(parts) == 3:
            value_1, value_2, value_3 = map(float, parts)
            rho.append(value_1)
            rho_a.append(value_2)
            a.append(value_3)

x = list(range(0, 300))

def power_law_regression(x,rho_c,C,beta):
    return rho_c + C * x**(1/beta)
    
rho_c,C,beta = optimize.curve_fit(power_law_regression, xdata=x,ydata=rho,p0=[0.3,1,0.5])[0]
rho_regression = power_law_regression(x,rho_c,C,beta)

print(f'rho_c = {rho_c:.4}, beta = {beta:.3}')

plt.figure(1)
plt.plot(x, rho, label='\u03C1')
plt.plot(x, rho_a, label='\u03C1\u2090')
plt.plot(x, rho_regression, label=f'{rho_c:.4} + x^{1/beta:.5}')
plt.xlabel('x')
plt.grid(True)
plt.legend(fontsize=14)
plt.savefig('active_sites.png')

plt.figure(2)
plt.plot(x, a)
plt.xlabel('x')
plt.ylabel('a')
plt.grid(True)
plt.savefig('activity.png')
