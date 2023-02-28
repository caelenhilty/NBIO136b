import random as r
import numpy as np
import matplotlib.pyplot as plt
import math as m

def AELIF(Iapps):
    ### ODEs ###
    ### voltage dot ###
    def Vdot(i):
        return (leakCond*(leakPot - V[i-1] + deltaTh*m.e**((V[i-1]-Vth)/deltaTh))-Isra[i-1]+Iapps[i])/membraneCapacitance

    ### dynamic spike rate adaptation current ###
    def IsraDot(i):
        return ((a*(V[i-1]-leakPot)-Isra[i-1])/taoSRA)

    ### declaring parameters ###
    leakPot = -60e-3 #V
    Vth = -50e-3 #V
    Vmax = 100e-3 #V
    Vreset = -80e-3 #V
    deltaTh = 2e-3 #V
    leakCond = 8e-9 #S
    membraneCapacitance = 100e-12 #F
    a = 10e-9 #S
    b = 0.5e-9 #A
    taoSRA = 50e-3 #s

    ### vectors ###
    dt = 0.02e-3
    duration = 200 #s
    time = np.arange(0,duration, dt)
    spikes = np.zeros(len(time))

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
            spikes[i] = 1
    
    return time, spikes

values = np.random.uniform(-0.5, 0.5, 40000) * 10e-9
time = np.linspace(start = 0, stop = 40e3*5, num = 10000000)
Iapps = np.zeros(len(time))
n = 0
for value in values:
    for i in range(250):
        Iapps[n+i] = value
    n += 249

times, spikes = AELIF(Iapps)
plt.figure(layout = 'constrained')
plt.plot(times, spikes)
plt.title("Spikes")
plt.xlabel("time (s)")
plt.show()