import math as mt
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from scipy.interpolate import interp1d


wde1=np.genfromtxt('printwde.dat')
z1=wde1[:,0]
w1=wde1[:,1]

l=len(z1)


wde2=np.genfromtxt('printwde_sm.dat')
z2=wde2[:,0]
w2=wde2[:,1]

u=np.zeros(l)
for i in range (0,l):
	u[i]=-1

plt.plot(z1,w1)
#plt.plot(z2,w2)
#plt.plot(z1,u)
plt.show()

