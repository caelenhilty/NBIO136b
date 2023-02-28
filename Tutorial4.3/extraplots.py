import random as r
import numpy as np
import matplotlib.pyplot as plt
import math as math

def PR(Glink = 50e-9):
    ### parameters ###
    somaArea = 1/3
    somaticLeakCond = somaArea * 5e-9 #S
    maxNaCond = somaArea * 3e-6 #S
    maxDelayedRectifierCond = somaArea * 2e-6 #S
    somaCapacitance = somaArea * 100e-12 #F
    somaIapp = 0

    dendriteArea = 1 - somaArea
    dendriticLeakCond = dendriteArea * 5e-9 #S
    maxCaCond = dendriteArea * 2e-6 #S
    maxKCaCond = dendriteArea * 2.5e-6 #S
    maxKAHPCond = dendriteArea * 40e-9 #S
    dendriteCapacitance = dendriteArea * 100e-12 #F
    dendriteIapp = 0

    linkCond = Glink #S

    NaPot = 60e-3 #V
    CaPot = 80e-3 #V
    taoCa = 50e-3 #S
    KPot = -75e-3 #V
    leakPot = -60e-3 #V

    kConvert = 5e6/dendriteArea #MC^-1

    ### vectors ###
    dt = 0.2e-5 #s
    duration = 2 #s
    time = np.arange(0,duration+dt, dt)

    Vs = np.ones(len(time)) * leakPot
    Vd = np.ones(len(time)) * leakPot
    m = np.ones(len(time)) * 0
    h = np.ones(len(time)) * 0
    n = np.ones(len(time)) * 0
    mCa = np.ones(len(time)) * 0
    mKCa =np.ones(len(time)) * 0
    mKAHP = np.ones(len(time)) * 0
    Ca = np.ones(len(time)) * 0 # concentration

    ### rate constants for gating variables ###
    def alphaM(Vs):
        return (320*10**3*(Vs+0.0469)) / (1-math.exp(-250*(Vs+0.0469)))
    def betaM(Vs):
        return (280*10**3*(Vs+0.0199)) / (math.exp(200*(Vs+0.0199))-1)

    def alphaH(Vs):
        return 128*math.exp(-(Vs+0.043)/0.018)
    def betaH(Vs):
        return 4000/(1+math.exp(-200*(Vs+0.020)))

    def alphaN(Vs):
        return (16e3*(Vs+0.0249))/(1-math.exp(-200*(Vs + 0.0249)))
    def betaN(Vs):
        return 250*math.exp(-25*(Vs+0.040))

    def alphaMCa(Vd):
        return 1600/(1+math.exp(-72*(Vd-0.005)))
    def betaMCa(Vd):
        return (2e4*(Vd+0.0089)) / (math.exp(200*(Vd+0.0089))-1)

    def alphaMKCa(Vd):
        if Vd > -0.010:
            return 2000*math.exp(-(Vd+0.0535)/0.027)
        else:
            return math.exp((Vd+0.05)/0.011-(Vd+0.0535)/0.027) / 0.018975
    def betaMKCa(Vd):
        if Vd > -0.010:
            return 0
        else:
            return 2000*math.exp(-(Vd+0.0535)/0.027) - alphaMKCa(Vd)

    def X(i):
        return min([4000*Ca[i], 1])

    def alphaMKAHP(i):
        return min([20, 20000*Ca[i]])
    betaMKAHP = 4

    ### ODEs ###
    def VsDot(i):
        return (somaticLeakCond*(leakPot-Vs[i]) + maxNaCond*m[i]**2*h[i]*(NaPot-Vs[i]) + maxDelayedRectifierCond*n[i]**2*(KPot-Vs[i]) + linkCond*(Vd[i]-Vs[i]) + somaIapp)/somaCapacitance
    def mDot(i):
        return alphaM(Vs[i])*(1-m[i]) - betaM(Vs[i])*m[i]
    def hDot(i):
        return alphaH(Vs[i])*(1-h[i]) - betaH(Vs[i])*h[i]
    def nDot(i):
        return alphaN(Vs[i])*(1-n[i]) - betaN(Vs[i])*n[i]

    def VdDot(i):
        return (dendriticLeakCond*(leakPot-Vd[i]) + maxCaCond*mCa[i]**2*(CaPot-Vd[i]) + maxKCaCond*mKCa[i]*X(i)*(KPot-Vd[i]) + maxKAHPCond*mKAHP[i]*(KPot-Vd[i]) - linkCond*(Vd[i]-Vs[i]) + dendriteIapp)/dendriteCapacitance
    def mCaDot(i):
        return alphaMCa(Vd[i])*(1-mCa[i]) - betaMCa(Vd[i])*mCa[i]
    def mKCaDot(i):
        return alphaMKCa(Vd[i])*(1-mKCa[i]) - betaMKCa(Vd[i])*mKCa[i]
    def mKAHPDot(i):
        return alphaMKAHP(i)*(1-mKAHP[i]) - betaMKAHP*mKAHP[i]

    def caDot(i):
        return -Ca[i]/taoCa+kConvert*maxCaCond*mCa[i]**2*(CaPot-Vd[i])

    ### Forward Euler ###
    canSpike = True
    for i in range(1, len(time)):
        Vs[i] = Vs[i-1] + dt*VsDot(i-1)
        m[i] = m[i-1] + dt*mDot(i-1)
        h[i] = h[i-1] + dt*hDot(i-1)
        n[i] = n[i-1] + dt*nDot(i-1)
        Vd[i] = Vd[i-1] + dt*VdDot(i-1)
        mCa[i] = mCa[i-1] + dt*mCaDot(i-1)
        mKCa[i] = mKCa[i-1] + dt*mKCaDot(i-1)
        mKAHP[i] = mKAHP[i-1] + dt*mKAHPDot(i-1)
        Ca[i] = Ca[i-1] + dt*caDot(i-1)
        if Vs[i] < -30e-3:
            canSpike = True
        if Vs[i] > -10e-3 and canSpike:
            print("spiked", time[i])
            canSpike = False
    
    return time, Vs, Vd, mCa, mKCa, mKAHP

plt.figure(layout = 'constrained')

plt.subplot(211)
time, Vs, Vd, mCa, mKCa, mKAHP = PR(Glink = 100e-9)
plt.plot(time, mKCa, label = "mKCa")
plt.plot(time, mCa, label = 'mCa')
plt.plot(time, mKAHP, label = 'mKAHP')
plt.ylabel('Gating Variable')
plt.xlabel('Time (s)')
plt.title("Pinksy Rinzel, $G_{link}$ = 100 nS")
plt.legend(loc = 'upper left')

plt.subplot(212)
plt.plot(time, Vs, label = "Somatic")
plt.plot(time, Vd, label = "Dendritic")
plt.ylabel('Membrane Potential (V)')
plt.xlabel('Time (s)')
plt.legend(loc = 'upper left')

plt.show()