import random as r
import numpy as np
import matplotlib.pyplot as plt

w = [1,2]
x0 = 1
v0 = 0

def X(v):
    return v

def V(x, w, v):
    return -w**2*x-0.5*v

dt = 10**-3 #s
simulation_duration = 10 #s
t = np.arange(0, simulation_duration + dt, dt)
num_steps = len(t)
x = np.zeros(num_steps)
x[0] = x0
v = np.zeros(num_steps)
v[0] = v0

for W in w:
    for i in range(1, num_steps):
        x[i] = x[i-1] + X(v[i-1]) * dt
        v[i] = v[i-1] + V(x[i-1], W, v[i-1]) * dt
    plt.plot(x,v, label = f"w = {W}")

plt.xlabel("x")
plt.ylabel("v")
plt.legend()
plt.show()
plt.close()