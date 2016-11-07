import numpy as np
import matplotlib.pyplot as plt
data = np.loadtxt("transprob", delimiter=",")
E = data[:,0]
T = data[:,1]
plt.plot(E,T,'r-.o')
plt.xlabel("$E/V_0$")
plt.ylabel("T")
plt.title("Transmission Probability by FFT")
plt.savefig("transprob1.png")
plt.show()
