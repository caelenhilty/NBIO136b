import random as r
import numpy as np
import matplotlib.pyplot as plt

#declaring parameters
leakPot = -70e-3 #V
totMemRes = 5e6 #Ohms
totMemCapac = 2e-9 #F

spikeThreshold = -50e-3 #V
resetPot = -65e-3 #V

#LIR equation, solved for Vm-dot
def f(x,i):
    return ((leakPot-x)/totMemRes + Iapp[i])/totMemCapac

#generate vectors
dt = 0.1 * 10**-3
time = np.arange(0,0.2+dt, dt)

V = np.ones(len(time)) * leakPot

I0 = 4.01e-9 #A
Iapp = np.ones(len(time)) * I0

tref = 2.6
spikes = []

#forward Euler method
for i in range(1, len(time)):
    tref+=dt
    if tref > 2.5e-3:
        x = V[i-1]
        V[i] = x + dt*f(x,i)
    else:
        V[i] = resetPot
    if V[i] > spikeThreshold:
        V[i] = resetPot
        tref = 0
        spikes.append(i)

for i in spikes:
    V[i] = 50e-3

#plot it
plt.figure(layout = 'constrained')
plt.title("Membrane Potential with Votage Clamp")
plt.plot(time, V, label = "$I_{app}$ = 4.01 pA")

plt.xlabel("Time (s)")
plt.ylabel("Membrane Potential (V)")
plt.legend()
plt.show()
plt.close()