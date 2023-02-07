import random as r
import numpy as np
import matplotlib.pyplot as plt

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
duration = 1.5 #s
time = np.arange(0,duration+dt, dt)
Iapps = np.zeros(len(time))
for i in range(5000, 10000):
    Iapps[i] = 500e-12
firingRates = np.zeros(len(Iapps))

V = np.ones(len(time)) * leakPot
G = np.zeros(len(time))

### voltage-dot ###
def Vdot(V, G,i):
    return ((leakPot - V)/membraneRes+G*(potassiumPot-V) + Iapps[i])/membraneCond
### conductance-dot ###
def Gdot(G):
    return -G/taoSRA

### forward Euler method ###
for i in range(1, len(time)):
    G[i] = G[i-1] + dt * Gdot(G[i-1])
    V[i] = V[i-1] + dt * Vdot(V[i-1], G[i],i)
    if V[i] > Vth:
        V[i] = Vreset
        G[i] += deltaGSRA


### plot it ###
plt.figure(layout = 'constrained')
plt.subplot(311)
plt.title("2.3.1(a) ELIF Current Pulse")
plt.plot(time, Iapps*10**12)
plt.ylabel("$I_{app} (pA)$")

plt.subplot(312)
plt.plot(time, V)
plt.ylabel("Membrane Voltage (V)")

plt.subplot(313)
plt.plot(time, G)
plt.ylabel("Spike Rate Adaptation Conductance (F)")
plt.xlabel("time (s)")

plt.show()
plt.close()