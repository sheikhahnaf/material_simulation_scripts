import numpy as np
import scipy
import scipy.integrate

from matplotlib import animation
from IPython.display import HTML
from matplotlib import pyplot as plt
from pydmd import DMD
from pydmd import FbDMD



A=np.loadtxt('dmd.txt',delimiter=",")
B=A.T
C=B[:,0:5000]

print(C.shape)
dmd = DMD(svd_rank=1, tlsq_rank=2, exact=True, opt=True)
dmd.fit(C)
fbdmd = FbDMD(exact=True)
fbdmd.fit(C)
print(fbdmd.reconstructed_data.shape)
fbdmd.dmd_time['dt'] *= .5
fbdmd.dmd_time['tend'] += 10

plt.plot(fbdmd.dmd_timesteps, fbdmd.dynamics.T.real)

plt.show()


print("Shape before manipulation: {}".format(dmd.reconstructed_data.shape))
dmd.dmd_time['dt'] *= .25
dmd.dmd_time['tend'] *= 3
print("Shape after manipulation: {}".format(dmd.reconstructed_data.shape))
