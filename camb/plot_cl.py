import math as mt
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from scipy.interpolate import interp1d


cls=np.genfromtxt('test_lensedCls.dat')
l=cls[:,0]
cl=cls[:,1]

plt.plot(l,cl)
plt.xscale('log')
plt.title('TT')
plt.xlabel('l')
plt.ylabel('l(l+1)CTT/2pi')
plt.show()

