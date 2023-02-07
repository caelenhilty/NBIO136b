import random as r
import numpy as np
import matplotlib.pyplot as plt

#declaring parameters
leakPot = -70e-3 #V
totMemRes = 5e6 #Ohms
totMemCapac = 2e-9 #F

resetPot = -65e-3 #V
KPot = -80e-3

#LIR equation, solved for Vm-dot
def f(x,i):
    return ((leakPot-x)/totMemRes + Grefs[i]*(KPot - x) + Iapp[i])/totMemCapac

#dynamic threshold equation
def Vth(V):
    return (-50e-3-V)/1e-3

#dynamic refractory conductance
def Gref(G):
    return (-G/0.2e-3)

#generate vectors
dt = 0.1 * 10**-3
time = np.arange(0,0.2+dt, dt)

V = np.ones(len(time)) * leakPot
Vths = np.ones(len(time)) * -50e-3
Grefs = np.zeros(len(time))
dG = 0.2e-6 #S

I0 = 4.1e-9 #A
Iapp = np.ones(len(time)) * I0

spikes = []

#forward Euler method
for i in range(1, len(time)):
    Grefs[i] = Grefs[i-1] + dt*Gref(Grefs[i-1])
    
    x = V[i-1]
    V[i] = x + dt*f(x,i)

    Vths[i] = Vths[i-1] + dt * Vth(Vths[i-1])
    spikeThreshold = Vths[i]
    if V[i] > spikeThreshold:
        Vths[i] = 200e-3
        Grefs[i] += dG
        spikes.append(i)

for i in spikes:
    V[i] = 50e-3

#plot it
plt.figure(layout = 'constrained')
plt.title("Membrane Potential with Dynamic Theshold")
plt.plot(time, V, label = "$I_{app}$ = 4.01 pA")
plt.xlabel("Time (s)")
plt.ylabel("Membrane Potential (V)")
plt.legend()
plt.show()
plt.close()