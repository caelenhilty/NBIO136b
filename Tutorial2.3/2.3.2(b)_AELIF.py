import random as r
import numpy as np
import matplotlib.pyplot as plt
import math as m
def AELIF(Iapp):
    ### voltage dot ###
    def Vdot(i):
        return (leakCond*(leakPot - V[i-1] + deltaTh*m.e**((V[i-1]-Vth)/deltaTh))-Isra[i-1]+Iapps[i])/membraneCapacitance

    ### dynamic spike rate adaptation current ###
    def IsraDot(i):
        return ((a*(V[i-1]-leakPot)-Isra[i-1])/taoSRA)

    ### declaring parameters ###
    leakPot = -75e-3 #V
    Vth = -50e-3 #V
    Vmax = 100e-3 #V
    Vreset = -80e-3 #V
    deltaTh = 2e-3 #V
    leakCond = 10e-9 #S
    membraneCapacitance = 100e-12 #F
    a = 2e-9 #S
    b = 0.02e-9 #A
    taoSRA = 200e-3 #s

    ### vectors ###
    dt = 0.1e-3
    duration = 5 #s
    time = np.arange(0,duration+dt, dt)
    Iapps = np.ones(len(time)) * Iapp
    spikes = []

    V = np.ones(len(time)) * leakPot
    Isra = np.zeros(len(time))

    ### forward Euler ###
    for i in range(1, len(time)):
        V[i] = V[i-1] + dt*Vdot(i)
        Isra[i] = Isra[i-1] + dt*IsraDot(i)
        if V[i] > Vmax:
            V[i-1] = Vmax
            V[i] = Vreset
            Isra[i] += b
            spikes.append(i)

    ### calculate first ISI ###
    if len(spikes) > 1:
        firstISI = time[spikes[1]] - time[spikes[0]]
    else:
        firstISI = 0

    ### calculate average ISI ###
    if len(spikes) > 1:
        ISIs = []
        prev = spikes[0]
        for spike in spikes[1:]:
            ISIs.append(time[spike] - time[prev])
            prev = spike
        sum = 0
        for ISI in ISIs:
            sum += ISI
        ISIaverage = sum / len(ISIs)
    else:
        ISIaverage = 0
    
    return firstISI, ISIaverage

Iapps = np.arange(0, 600, 20) * 10**-12
firstISIs = []
averageISIs = []
for Iapp in Iapps:
    firstISI, ISIaverage = AELIF(Iapp)
    if firstISI > 0:
        firstISIs.append(1/firstISI)
    else:
        firstISIs.append(0)
    if ISIaverage > 0:
        averageISIs.append(1/ISIaverage)
    else:
        averageISIs.append(0)

### plot it ###
plt.figure(layout = 'constrained')
plt.plot(Iapps*10**12, firstISIs,"ro" , label = "1/ISI(1)")
plt.plot(Iapps*10**12, averageISIs, label = "final rate")
plt.xlabel("$I_{app} (pA)$")
plt.ylabel("Spike Rate (Hz)")
plt.title("2.3.2(b) f-I curve")
plt.legend()
plt.show()
plt.close()