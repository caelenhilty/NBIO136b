import random as r
import numpy as np
import matplotlib.pyplot as plt
import math as math


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

def alphaMKAHP(i):
    return min([20, 20000*Ca[i]])
betaMKAHP = 4

V = np.arange(-85, 50, 5) * 10**-3
Ca = np.linspace(0, 2, len(V)) * 10**-3


alpham = np.zeros(len(V))
betam = np.zeros(len(V))
alphah = np.zeros(len(V))
betah = np.zeros(len(V))
alphan = np.zeros(len(V))
betan = np.zeros(len(V))
alphamca = np.zeros(len(V))
betamca = np.zeros(len(V))
alphamkca = np.zeros(len(V))
betamkca = np.zeros(len(V))
alphamkahp = np.zeros(len(V))
betamkahp = np.zeros(len(V))

for i in range(len(V)):
    v = V[i]
    alpham[i] = alphaM(v)
    betam[i] = betaM(v)
    alphah[i] = alphaH(v)
    betah[i] = betaH(v)
    alphan[i] = alphaN(v)
    betan[i] = betaN(v)
    alphamca[i] = alphaMCa(v)
    betamca[i] = betaMCa(v)
    alphamkca[i] = alphaMKCa(v)
    betamkca[i] = betaMKCa(v)
    alphamkahp[i] = alphaMKAHP(i)
    betamkahp[i] = betaMKAHP

plt.figure(layout = 'constrained')
plt.subplot(121)
plt.title("Voltage Dependence")
plt.xlabel("Membrane voltage (V)")
plt.ylabel("Gating Variable Rate Constant")
plt.plot(V, alpham, label = r'$\alpha_{m}$')
plt.plot(V, betam, label = r'$\beta_{m}$')
plt.plot(V, alphah, label = r'$\alpha_{h}$')
plt.plot(V, betah, label = r'$\beta_{h}$')
plt.plot(V, alphan, label = r'$\alpha_{n}$')
plt.plot(V, betan, label = r'$\beta_{n}$')
plt.plot(V, alphamca, label = r'$\alpha_{mCa}$')
plt.plot(V, betamca, label = r'$\beta_{mCa}$')
plt.plot(V, alphamkca, label = r'$\alpha_{mKCa}$')
plt.plot(V, betamkca, label = r'$\beta_{mKCa}$')
plt.legend()

plt.subplot(122)
plt.title("[Ca] Dependence")
plt.xlabel("[Ca] (M)")
plt.ylabel("Gating Variable Rate Constant")
plt.plot(V, alphamkahp, label = r'$\alpha_{mKAHP}$')
plt.plot(V, betamkahp, label = r'$\beta_{mKAHP}$')
plt.legend()
plt.show()

