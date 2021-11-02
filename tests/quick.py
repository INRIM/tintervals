

import numpy as np
from matplotlib import pyplot
ion()
close('all')



t = np.arange(100) + np.random.normal(0, 0.05, 100)

t[t>50] += 50

t[30] +=1


deltat = median(diff(t))

tp = np.roll(t,1) + deltat
tm = np.roll(t,-1) - deltat


tnorm = median(column_stack((tp, t, tm)), axis=-1)

figure()
plot(tnorm, t-tnorm)

figure()
plot(diff(tnorm))