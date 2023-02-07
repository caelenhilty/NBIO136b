import random as r
import numpy as np
import matplotlib.pyplot as plt

#declaring parameters
leakPot = -70*10**-3 #V
totMemRes = 5*10**6 #Ohms
totMemCapac = 2*10**-9 #F

spikeThreshold = -50*10**-3 #V
resetPot = -65*10**-3 #V

#LIR equation, solved for Vm-dot
def f(x,i):
    return ((leakPot-x)/totMemRes + Iapps[i])/totMemCapac

#generate vectors
dt = 0.1 * 10**-3
duration = 2
time = np.arange(0,duration+dt, dt)
Iapps = np.arange(3,6,0.1) * 10**-9
firingRates = np.zeros(len(Iapps))

#forward Euler method
for i in range(len(Iapps)):
    V = np.ones(len(time)) * leakPot
    spikes = 0
    for j in range(1, len(time)):
        x = V[j-1]
        V[j] = x + dt*f(x,i)
        if V[j] > spikeThreshold:
            V[j] = resetPot
            spikes += 1
    firingRates[i] = spikes/duration

#plot it
plt.figure(layout = 'constrained')
plt.title("2.1.1(c) f-I curve")
plt.plot(Iapps*10**9, firingRates)
plt.xlabel("$I_{app}$ (pA)")
plt.ylabel("Average Firing Rate (Hz)")
plt.show()
plt.close()