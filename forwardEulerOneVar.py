import random as r
import numpy as np
import matplotlib.pyplot as plt

K = 1
A = 1
B = 1
TExt = 25

def f(x):
    return (-A+B*(TExt-x))/K

x0 = 10
dt = 10**-3 #s
simulation_duration = 60 #s

t = np.arange(0, simulation_duration + dt, dt)
num_steps = len(t)
x = np.zeros(num_steps)
x[0] = x0

for i in range(1, num_steps):
    x[i] = x[i-1] + f(x[i-1]) * dt

plt.plot(t,x)
plt.show()
plt.close()