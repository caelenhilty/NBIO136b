import random as r
import numpy as np
import matplotlib.pyplot as plt

### parameters ###
Sinitials = np.array([0.5,0.25,0.75])
alpha = 1
beta = 2

### ODE, f(s) = s-dot ###
def f(s):
    return alpha*(1-s) - beta*s

### vectors ###
dt = 10**-3
simulation_duration = 3 #s
t = np.arange(0, simulation_duration + dt, dt)
num_steps = len(t)

### forward Euler
for Sinitial in Sinitials:
    S = np.ones(num_steps) * Sinitial
    for i in range(1, num_steps):
        S[i] = S[i-1] + f(S[i-1])*dt
    plt.plot(t,S, label = f"S = {Sinitial}")

plt.xlabel("t (s)")
plt.ylabel("S = $N_O / N_T$")
plt.legend()
plt.show()
plt.close()