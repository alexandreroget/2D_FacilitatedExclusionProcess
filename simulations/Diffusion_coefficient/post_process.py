import matplotlib.pyplot as plt

import numpy as np

t = []
N = []

lambdada = np.array([0.1, 0.25, 0.5, 0.75])
epsilon=0.01

plt.figure()
colors = ['#16679c', '#c76005', '#1d7c1d', '#c6393c']

for i in range(len(lambdada)):
    t.clear()
    N.clear()
    
    filename = "data/diffusion_current_lambda0" + str(int(100*lambdada[i])) + ".txt"

    with open(filename, 'r') as file:
        for line in file:
            value_1, value_2 = map(float, line.strip().split(';'))
            t.append(value_1)
            N.append(value_2)
            
    label = "\u03BB = " + str(lambdada[i])
    # plt.plot(t, N, label=label, color=colors[i])
    
    t2 = []
    N2 = []
    div = 5000
    n_points = int(len(t)/div)
    for j in range(n_points):
    	t2.append(t[j*div])
    	N2.append(N[j*div])

    plt.plot(t2, N2, label=label, linewidth='2', color=colors[i])
    
    
plt.plot(t,np.array(t)*epsilon,label="N=\u03B5t",linewidth='2',color='black')

plt.title('Diffusion current evolution')
plt.xlabel('t')
plt.ylabel('N')
plt.grid(True)
plt.legend(fontsize=14)

plt.savefig('diffusion_current.png')
