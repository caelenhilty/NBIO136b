import random as r
import numpy as np
import matplotlib.pyplot as plt

#declaring parameters
leakPot = -70*10**-3 #V
totMemRes = 5*10**6 #Ohms
totMemCapac = 2*10**-9 #F

spikeThreshold = -50*10**-3 #V
resetPot = -65*10**-3 #V

#LIR equation, solved for Vm-dot
def f(x,i):
    return ((leakPot-x)/totMemRes + Iapp[i])/totMemCapac

#generate vectors
dt = 0.1 * 10**-3
time = np.arange(0,0.2+dt, dt)

V = np.ones(len(time)) * leakPot

I0 = 3.999*10**-9 #A
Iapp = np.ones(len(time)) * I0

#forward Euler method
for i in range(1, len(time)):
    x = V[i-1]
    V[i] = x + dt*f(x,i)
    if V[i] > spikeThreshold:
        V[i] = resetPot

#plot it
plt.figure(layout = 'constrained')
plt.title("2.1.1(b) Minimum applied current")
plt.plot(time, V, label = "$I_{app}$ = 3.999 pA")

###   do it again   ###
V = np.ones(len(time)) * leakPot
I0 = 4.001*10**-9 #A
Iapp = np.ones(len(time)) * I0

#forward Euler method
for i in range(1, len(time)):
    x = V[i-1]
    V[i] = x + dt*f(x,i)
    if V[i] > spikeThreshold:
        V[i] = resetPot

plt.plot(time, V, label = "$I_{app}$ = 4.001 pA")

plt.xlabel("Time (s)")
plt.ylabel("Membrane Potential (V)")
plt.legend()
plt.show()
plt.close()