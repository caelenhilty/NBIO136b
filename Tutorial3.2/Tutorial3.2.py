import random as r
import numpy as np
import matplotlib.pyplot as plt
import math as m
import numba as nb

def AELIF(sigma = 50e-12, dt = 0.01e-3, duration = 100, b = 0, Iconstant = 0):
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

    Iapps = np.random.normal(0, sigma/m.sqrt(dt), len(time)) + Iconstant
    spikes = []

    ### forward Euler ###
    for i in range(1, len(time)):
        V[i] = V[i-1] + dt*Vdot(i)
        Isra[i] = Isra[i-1] + dt*IsraDot(i)
        if V[i] > Vmax:
            V[i-1] = Vmax
            V[i] = Vreset
            Isra[i] += b
            spikes.append(time[i])

    return time, spikes, Iapps, V

def fanoFactor(time, spikes, window,dt):
    bins = np.zeros(int(len(time)*dt//window)+1)
    for spike in spikes:
        bins[int(spike//window)] += 1
    variance = np.var(bins)
    mean = np.average(bins)
    return variance/mean

plt.figure(layout = "constrained")

windows = np.linspace(0.01, 1, num = 1000)
Iconstants = [0, 0.1e-9, 0.2e-9]
a=1
for I in Iconstants:
    fanoFactors = np.zeros(len(windows))
    time, spikes, Iapps, V = AELIF(sigma = 20e-12, b = 0e-9, Iconstant=I)
    for i in range(len(windows)):
        fanoFactors[i] = fanoFactor(time, spikes, windows[i], dt = 0.01e-3)
    plt.subplot(3,1,a)
    plt.plot(windows, fanoFactors)
    plt.title('Fano Factor vs Window Size, b = 0nA, ' + r"$\sigma = 20 pA\cdot s^{0.5}$, " + r'$I_{constant}$ =  ' + f"{I} A")
    plt.ylabel('Fano Factor')
    plt.xlabel("Window Size (s)")
    a+=1
plt.show()