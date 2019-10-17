import numpy as np
import matplotlib.pyplot as plt

N = 1000
alpha = 2

x1 = np.linspace(-10, 10, N)
x2 = x1.copy()
y1 = x1.copy()
y2 = x1.copy()
z1 = x1.copy()
z2 = x1.copy()

r1 = np.sqrt(x1**2+y1**2+z1**2)
r2 = np.sqrt(x2**2+y2**2+z2**2)
r12 = 1/np.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
f = np.exp(-alpha*(r1))

plt.figure()
plt.plot(x1, f)
plt.grid()


plt.show()
