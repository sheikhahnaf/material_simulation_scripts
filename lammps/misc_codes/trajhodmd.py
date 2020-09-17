import matplotlib.pyplot as plt
import numpy as np
import time

from pydmd import HODMD



A=np.loadtxt('dmd.txt',delimiter=",")
B=A.T
C=B[:,0:3000]
original=C

print(C.shape)
hodmd = HODMD(svd_rank=0, exact=True, opt=True, d=30).fit(C)
print(hodmd.reconstructed_data.shape)
hodmd.plot_eigs()
hodmd.original_time['dt'] = hodmd.dmd_time['dt'] 
hodmd.original_time['t0'] = hodmd.dmd_time['t0'] 
hodmd.original_time['tend'] = hodmd.dmd_time['tend'] 

plt.plot(hodmd.original_timesteps, C[0,:], '.', label='snapshots')
plt.plot(hodmd.original_timesteps, original[0,:], '-', label='original function')
plt.plot(hodmd.dmd_timesteps, hodmd.reconstructed_data[0].real, '--', label='DMD output')
plt.legend()
plt.show()
hodmd.dmd_time['tend'] = 4000

fig = plt.figure(figsize=(15, 5))
plt.plot(hodmd.original_timesteps, C[0,:], '.', label='snapshots')
plt.plot(hodmd.dmd_timesteps, B[0,0:4001], '-', label='original function')
plt.plot(hodmd.dmd_timesteps, hodmd.reconstructed_data[0].real, '--', label='DMD output')
plt.legend()
plt.show()


print("Shape after manipulation: {}".format(hodmd.reconstructed_data.shape))
