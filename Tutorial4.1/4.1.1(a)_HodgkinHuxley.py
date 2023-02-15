import random as r
import numpy as np
import matplotlib.pyplot as plt
import math as math

### parameters ###
Gleak = 30e-9 #S
maxGNa = 12e-6 #S
maxGK = 3.6e-6 #S
NaPot = 45e-3 #V
KPot = -82e-3 #V
leakPot = -60e-3 #V
Cmembrane = 100e-12 #F

### vectors ###
dt = 0.1e-3
duration = 0.35 #s
time = np.arange(0,duration+dt, dt)
Iapps = np.zeros(len(time))
V = np.ones(len(time)) * leakPot
m = np.zeros(len(time))
h = np.zeros(len(time))
n = np.zeros(len(time))

### gating variable rate constants ###
def alphaM(V):
    am = (1e5*(-V-0.045))/(math.exp(100*(-V-0.045))-1)
    # print(f"alphaM: {am:.4f}")
    return am
def betaM(V):
    bm = 4e3*math.exp((-V-0.07)/0.018)
    # print(f"betaM: {bm:.4f}")
    return bm

def alphaH(V):
    ah = 70*math.exp(50*(-V-0.07))
    # print(f"alphaH: {ah:.4f}")
    return ah
def betaH(V):
    bh = 1e3 / (1+math.exp(100*(-V-0.04)))
    # print(f"betaH: {bh:.4f}")
    return bh

def alphaN(V):
    if V == -0.06:
        return 0
    an = (1e4*(-V-0.06))/(math.exp(100*(-V-0.06))-1)
    # print(f"alphaN: {an:.4f}")
    return an
def betaN(V):
    bn = 125*math.exp((-V-0.07)/0.08)
    # print(f"betaN: {bn:.4f}")
    return bn

### system of ODEs ###
def mdot(i):
    return alphaM(V[i])*(1-m[i]) - betaM(V[i])*m[i]
def hdot(i):
    return alphaH(V[i])*(1-h[i]) - betaH(V[i])*h[i]
def ndot(i):
    return alphaN(V[i])*(1-n[i]) - betaN(V[i])*n[i]
def Vdot(i):
    return (Gleak*(leakPot-V[i]) + maxGNa*m[i]**3*h[i]*(NaPot-V[i]) + maxGK*n[i]**4*(KPot-V[i]) + Iapps[i])/Cmembrane

### Forward Euler ###
for i in range(1, len(time)):
    V[i] = V[i-1] + dt*Vdot(i-1)
    m[i] = m[i-1] + dt*mdot(i-1)
    h[i] = h[i-1] + dt*hdot(i-1)
    n[i] = n[i-1] + dt*ndot(i-1)
    # print(f"V: {V[i]:.5f} m: {m[i]:.5f} h: {h[i]:.5f} n: {n[i]:.5f}")

### plot it ###
plt.figure()
plt.plot(time, V)
plt.xlabel("time (s)")
plt.ylabel("Membrane Potential (V)")
plt.title("4.1.1(a) Hodgkin-Huxley, $I_{app}$ = 0")
plt.show()