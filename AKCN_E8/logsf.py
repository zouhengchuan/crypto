from scipy.stats import chi2
from math import log
import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(100,500,400)
y = np.zeros(shape=(1,400))
y = chi2.logsf(x, 8) / log(2)
y1 = -0.71*x + 15
z = chi2.sf(x,8) / log(2)

plt.plot(x,y,label='logsf')
plt.plot(x,y1,label='-0.71x + 15')
# plt.plot(x,z,label='sf')
plt.title('Log Survival Function of Chi-Squared Distribution')
plt.xlabel('(dis/s)^2')
plt.ylabel('logsf(x,8)/log(2)')
plt.legend()
plt.grid(True)
plt.show()