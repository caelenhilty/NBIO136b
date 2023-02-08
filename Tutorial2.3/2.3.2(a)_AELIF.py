import random as r
import numpy as np
import matplotlib.pyplot as plt
import math as m

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
duration = 1.5 #s
time = np.arange(0,duration+dt, dt)
Iapps = np.zeros(len(time))
for i in range(5000, 10000):
    Iapps[i] = 280e-12

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

plt.figure(layout = 'constrained')
plt.subplot(211)
plt.title("2.3.2(a) AELIF Current Pulse")
plt.plot(time, Iapps*10**12)
plt.ylabel("$I_{app} (pA)$")

plt.subplot(212)
plt.plot(time, V)
plt.ylabel("Membrane Voltage (V)")
plt.xlabel("time (s)")

plt.show()
plt.close()