import random as r
import numpy as np
import matplotlib.pyplot as plt
from math import log, e, sqrt

#declaring parameters
leakPot = -70*10**-3 #V
totMemRes = 5*10**6 #Ohms
totMemCapac = 2*10**-9 #F

spikeThreshold = -50*10**-3 #V
resetPot = -65*10**-3 #V

Tm = totMemCapac * totMemRes #timescale
sigma = 0.01

#LIR equation, solved for Vm-dot
def f(x,i):
    return ((leakPot-x)/totMemRes + Iapps[i])/totMemCapac

#generate vectors
dt = 0.1 * 10**-3
duration = 2
time = np.arange(0,duration+dt, dt)
Iapps = np.arange(3,6,0.1) * 10**-9



### first trial ###
firingRates = np.zeros(len(Iapps))
theoreticalRates = np.zeros(len(Iapps))

#forward Euler method
for i in range(len(Iapps)):
    V = np.ones(len(time)) * leakPot
    spikes = 0
    for j in range(1, len(time)):
        x = V[j-1]
        V[j] = x + dt*f(x,i) + r.gauss(0,1) * sigma * sqrt(dt)
        if V[j] > spikeThreshold:
            V[j] = resetPot
            spikes += 1
    firingRates[i] = spikes/duration

    #theoretical curve
    if Iapps[i]*totMemRes+leakPot-resetPot > 0 and Iapps[i]*totMemRes+leakPot-spikeThreshold > 0:
        theoreticalRates[i] = 1/(Tm*log(Iapps[i]*totMemRes+leakPot-resetPot, e)-Tm*log(Iapps[i]*totMemRes+leakPot-spikeThreshold, e))

plt.figure(layout = 'constrained')
plt.plot(Iapps*10**9, firingRates, label = "$\sigma$ = 0.01")



### second trial ###
Iapps = np.arange(3,6,0.1) * 10**-9
firingRates = np.zeros(len(Iapps))
sigma = 0.1

#forward Euler method
for i in range(len(Iapps)):
    V = np.ones(len(time)) * leakPot
    spikes = 0
    for j in range(1, len(time)):
        x = V[j-1]
        V[j] = x + dt*f(x,i) + r.gauss(0,1) * sigma * sqrt(dt)
        if V[j] > spikeThreshold:
            V[j] = resetPot
            spikes += 1
    firingRates[i] = spikes/duration
plt.plot(Iapps*10**9, firingRates, label = "$\sigma$ = 0.1")



#plot it
plt.plot(Iapps*10**9, theoreticalRates, "r--", label = "analytical solution")
plt.title("2.1.2(b) f-I curve with noise")
plt.xlabel("$I_{app}$ (pA)")
plt.ylabel("Average Firing Rate (Hz)")
plt.legend()
plt.show()
plt.close()