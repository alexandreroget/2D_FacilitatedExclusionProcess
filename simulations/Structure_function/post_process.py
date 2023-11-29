import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from scipy import optimize

# We want to fit S(k) = A (rho-rho_c)^gamma + B |k|^alpha
def Sfit(rhok, A, rho_c, gamma, B, alpha):
    rho,k = rhok
    return A*abs(rho-rho_c)**gamma + B*k**alpha

L = 100
files_array = np.arange(0,20)
rho = []
k = [] 
S = []

for i in range(np.size(files_array)):

  filename = "data/S" + str(files_array[i]) + ".txt"
  file = open(filename,'r')

  Lines = file.read().splitlines()

  for j in range(np.size(Lines)-1):
    line = Lines[j+1]
    splited = line.split(',')
    
    x = (float(splited[0]))
    k1 = 2*np.pi*x/L
    
    y = (float(splited[1]))
    k2 = 2*np.pi*y/L
    
    if 0.1**2 < k1**2 + k2**2 < 0.5**2:
        k.append((k1**2 + k2**2)**0.5)
        S.append(float(splited[2]))
        rho.append(float(Lines[0]))
    file.close()
    
A,rho_c,gamma,B,alpha = optimize.curve_fit(Sfit,xdata=(rho,k),ydata=S,p0=[1,0.336,1,1,0.45])[0]

print(f'S(k) = {A:.2f}(rho-{rho_c:.4f})^{gamma:.2f} + {B:.2f}|k|^{alpha:.2f}')
print(f'rho_c = {rho_c:.4f}, gamma = {gamma:.2f}, alpha = {alpha:.2f}, zeta = {(1-alpha/2):.2f}')


files_array = [5, 6, 8, 10, 19] # The graph looks nicer without the densities with very large xi

fig,ax = plt.subplots()
for i in range(np.size(files_array)):
  k=[]
  S=[]
    
  filename = "data/S" + str(files_array[i]) + ".txt"
  file = open(filename,'r')

  Lines = file.read().splitlines()
  rho = float(Lines[0])

  for j in range(np.size(Lines)-1):
    line = Lines[j+1]
    splited = line.split(',')
    
    x = (float(splited[0]))
    k1 = 2*np.pi*x/L
    
    y = (float(splited[1]))
    k2 = 2*np.pi*y/L
    
    if 0 < k1**2 + k2**2 < 1**2:
        k.append((k1**2 + k2**2)**0.5)
        S.append(float(splited[2]))
  
  chi = A*abs(rho-rho_c)**gamma
  xi = (B/A)**(1/alpha) * abs(rho-rho_c)**(-gamma/alpha)
  
  label = "\u03C1 = " + str(rho)
  ax.plot(xi*np.array(k),S/chi,'+',label=label)

xik = range(1300)
Schi = 1+xik**alpha
ax.plot(xik,Schi,'black')

ax.set(xlabel=f'{(B/A)**(1/alpha):.2f} (\u03C1-{rho_c:.4f})^{-gamma/alpha:.2f} k')
ax.set(ylabel=f'{A:.2f}(\u03C1-{rho_c:.4f})^{gamma:.2f} S')
ax.set(title = f'Renormalized S(k) for different densities near criticality\nIn black: S(k) = {A:.2f}(\u03C1-{rho_c:.4f})^{gamma:.2f} + {B:.2f}|k|^{alpha:.2f}')
plt.legend(fontsize=12)

fig.savefig('collapse.png')
