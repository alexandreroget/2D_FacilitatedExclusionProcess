import numpy as np
import matplotlib.pyplot as plt

L = [5, 10, 20, 40, 60, 80, 100]

for i in range(2):
    abs_time = []
    for j in range(np.size(L)):
        if i == 1: 
            filename = f"data/L{L[j]}_subcritical.txt"
        else:
            filename = f"data/L{L[j]}_supercritical.txt"
        
        file = open(filename,'r')
        
        Lines = file.read().splitlines()
        rho = float(Lines[0])
        
        values = []
        for k in range(np.size(Lines)-1):
            value = float(Lines[k].strip())
            values.append(value)

        median = np.median(values)
        abs_time.append(median)

    plt.figure(i)

    if i == 1: 
        plt.plot(np.log10(L), abs_time, marker='o')
        title = 'Absorption time evolution in subcritical regime (\u03C1 = ' + str(rho) + ')'
        plt.title(title)
        plt.xlabel('log(L)')
        plt.ylabel('t')
        plt.grid(True)

        plt.savefig('abs_time_subcritical.png')
    else:
        plt.plot(np.log10(L), np.log10(abs_time), marker='o')
        title = 'Absorption time evolution in supercritical regime (\u03C1 = ' + str(rho) + ')'
        plt.title(title)
        plt.xlabel('log(L)')
        plt.ylabel('log(t)')
        plt.grid(True)

        plt.savefig('abs_time_supercritical.png')
