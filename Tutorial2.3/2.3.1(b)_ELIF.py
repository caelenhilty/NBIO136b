import random as r
import numpy as np
import matplotlib.pyplot as plt
def ELIF(Iapp):
    ### voltage-dot ###
    def Vdot(V, G,i):
        return ((leakPot - V)/membraneRes+G*(potassiumPot-V) + Iapps[i])/membraneCond
    ### conductance-dot ###
    def Gdot(G):
        return -G/taoSRA

    ### declaring parameters ###
    leakPot = -75e-3 #V
    Vth = -50e-3 #V
    Vreset = -80e-3 #V
    membraneRes = 100e6 #Ohms
    membraneCond = 100e-12 #F
    potassiumPot = -80e-3 #V
    deltaGSRA = 1e-9 #S
    taoSRA = 200e-3 #s

    ### vectors ###
    dt = 0.1e-3
    duration = 5 #s
    time = np.arange(0,duration+dt, dt)
    Iapps = np.ones(len(time))*Iapp

    V = np.ones(len(time)) * leakPot
    G = np.zeros(len(time))
    spikes=[]

    ### forward Euler method ###
    for i in range(1, len(time)):
        G[i] = G[i-1] + dt * Gdot(G[i-1])
        V[i] = V[i-1] + dt * Vdot(V[i-1], G[i-1],i)
        if V[i] > Vth:
            V[i] = Vreset
            G[i] += deltaGSRA
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
    firstISI, ISIaverage = ELIF(Iapp)
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
plt.title("2.3.1(b) f-I curve")
plt.legend()
plt.show()
plt.close()