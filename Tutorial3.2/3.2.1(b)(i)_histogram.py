import random as r
import numpy as np
import matplotlib.pyplot as plt
import math as m
import numba as nb

def AELIF(sigma = 50e-12, dt = 0.01e-3, duration = 100, b = 0):
    ### ODEs ###
    ### voltage dot ###
    def Vdot(i):
        return (leakCond*(leakPot - V[i-1] + deltaTh*m.e**((V[i-1]-Vth)/deltaTh))-Isra[i-1]+Iapps[i])/membraneCapacitance

    ### dynamic spike rate adaptation current ###
    def IsraDot(i):
        return ((a*(V[i-1]-leakPot)-Isra[i-1])/taoSRA)

    ### declaring parameters ###
    leakPot = -70e-3 #V
    Vth = -50e-3 #V
    Vmax = 100e-3 #V
    Vreset = -80e-3 #V
    deltaTh = 2e-3 #V
    leakCond = 10e-9 #S
    membraneCapacitance = 100e-12 #F
    a = 2e-9 #S
    taoSRA = 150e-3 #s

    ### vectors ###
    time = np.arange(0,duration, dt)

    V = np.ones(len(time)) * leakPot
    Isra = np.zeros(len(time))

    Iapps = np.random.normal(0, sigma/m.sqrt(dt), len(time))
    spikes = []

    ### forward Euler ###
    for i in range(1, len(time)):
        V[i] = V[i-1] + dt*Vdot(i)
        Isra[i] = Isra[i-1] + dt*IsraDot(i)
        if V[i] > Vmax:
            V[i-1] = Vmax
            V[i] = Vreset
            Isra[i] += b
            spikes.append(i)

    ### calculate ISIs ###
    if len(spikes) > 0:
        ISIs = np.zeros(len(spikes))
        prev = spikes[0]
        index=0
        for spike in spikes:
            ISIs[index] = (time[spike] - time[prev])
            prev = spike
            index+=1

    return time, spikes, ISIs, Iapps, V

time, spikes, ISIs, Iapps, V = AELIF(50e-12, 0.01e-3, 100, b=1e-9)
plt.figure(layout = "constrained")
plt.hist(ISIs, bins = 25)
plt.ylabel('frequency')
plt.xlabel("ISI (s)")
plt.title(r'ISI Distribution, $\sigma = 50 pA\cdot s^{0.5}$')
plt.show()