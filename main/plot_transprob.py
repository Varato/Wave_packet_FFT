import numpy as np
import matplotlib.pyplot as plt
data = np.loadtxt("transprob", delimiter=",")
E = data[:,0]
T = data[:,1]
plt.plot(E,T,'bo')
plt.show()