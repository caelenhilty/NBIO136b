import random as r
import numpy as np
import matplotlib.pyplot as plt
import math as math

def PIR(baseline, step):

    ### parameters ###
    Gleak = 10e-9 #S
    maxGNa = 3.6e-6 #S
    maxGK = 1.6e-6 #S
    maxGT = 0.22e-6 #S
    NaPot = 55e-3 #V
    KPot = -90e-3 #V
    CaPot = 120e-3 #V
    leakPot = -70e-3 #V
    membraneCapacitance = 100e-12 #F

    ### vectors ###
    dt = 0.1e-5 #s
    duration = 0.75 #s
    time = np.arange(0,duration+dt, dt)

    V = np.ones(len(time)) * -0.065
    h = np.ones(len(time)) * 0
    n = np.ones(len(time)) * 0
    hT = np.ones(len(time)) * 0

    ### applied current calculations ###
    Iapps = np.ones(len(time)) * baseline
    for i in range(int(0.25/dt), int((0.4+dt)/dt)):
        Iapps[i] = step

    ### gating variable rate constants ###
    def alphaM(V):
        return (10**5*(V+0.035))/(1-math.exp(-100*(V+0.035)))
    def betaM(V):
        return 4000*math.exp(-(V+0.06)/0.018)
    def mInf(V):
        aM = alphaM(V)
        return aM/(aM + betaM(V))

    def alphaH(V):
        return 350*math.exp(-50*(V+0.058))
    def betaH(V):
        return 5000/(1+math.exp(-100*(V+0.028)))

    def alphaN(V):
        return (5e4*(V+0.034))/(1-math.exp(-100*(V+0.034)))
    def betaN(V):
        return 625*math.exp(-12.5*(V+0.044))

    def mTinf(V):
        return 1/(1+math.exp(-(V+0.052)/0.0074))

    def hTinf(V):
        return 1/(1+math.exp(500*(V+0.076)))
    def taoHT(V):
        if V < -0.080:
            return 0.001*math.exp(15*(V+0.467))
        else:
            return 0.028 + 0.001*math.exp(-(V+0.022)/0.0105)

    ### system of ODEs ###
    def hdot(i):
        return alphaH(V[i])*(1-h[i]) - betaH(V[i])*h[i]
    def ndot(i):
        return alphaN(V[i])*(1-n[i]) - betaN(V[i])*n[i]
    def hTdot(i):
        return (hTinf(V[i]) - hT[i])/taoHT(V[i])
    def Vdot(i):
        v = V[i]
        return (Gleak*(leakPot-v) + maxGNa*mInf(v)**3*h[i]*(NaPot-v) + maxGK*n[i]**4*(KPot-v)+maxGT*mTinf(v)**2*hT[i]*(CaPot-v) + Iapps[i])/membraneCapacitance

    ### forward Euler method ###
    for i in range(1, len(time)):
        V[i] = V[i-1] + dt*Vdot(i-1)
        h[i] = h[i-1] + dt*hdot(i-1)
        n[i] = n[i-1] + dt*ndot(i-1)
        hT[i] = hT[i-1] + dt*hTdot(i-1)

    return time, V

time, V = PIR(0, 50e-12)

plt.plot(time,V)
plt.show()
plt.close()