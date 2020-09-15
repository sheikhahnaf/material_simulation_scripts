# VDOS - Vibrational density of states
# This script makes use of LAMMPS MD trajectories which
# include velocity information. The dump file format should
# be of the following:
# ID TYPE x y z vx vy vz.
# The use of numpy.correlate and numpy.fft are implemented
# for efficiency. However, one should be cautious of memory
# usage for large number of atoms or long times.

# Author: Stefan Bringuier
# Affiliation: University of Arizona - Dept. MSE
# Date: Feb. 12 2014

# SEE if __name__ == "__main__":
# for options to run

#TODO - Implement C++ version
#TODO - Input script driven
#TODO - Clean up routines if where possible

#!/usr/bin/python

import numpy as np
from scipy import fftpack
from scipy import signal
import matplotlib.pyplot as plt
import glob

def read_dump(FileName,NumAtoms,Nfield, SnapShots):
    """read lammps dumpfile with header (not saved)"""

    File = open(FileName,'r')
    data = np.ndarray((NumAtoms,Nfield,SnapShots),dtype=float)
    
    t = 0
    while (t < SnapShots):
        #read header
        h1 = File.readline()
        time = File.readline()
        h2 = File.readline()
        numatoms = File.readline()
        h3 = File.readline()
        xlen = File.readline()
        ylen = File.readline()
        zlen = File.readline()
        h4 = File.readline()

        for a in range(NumAtoms):
            #Read string -> strip '\n' char -> split into new list
            line = File.readline().strip('\n').split()
            data[a,:,t] = line
            
        t += 1

    File.close()
    return data


def autocorr(X):
    """ the convolution is actually being done here
    meaning from -inf to inf so we only want half the
    array"""

    result = np.correlate(X, X, mode='full')
    return result[result.size/2:]

    
def fft_autocorr(AutoCorr,dt):
    """FFT of autocorrelation function"""
    #fft_arry = fftpack.dct(AutoCorr)*dt
    fft_arry = np.fft.rfft(AutoCorr)*dt
    return fft_arry


def  process(corlen,num_atoms,nfield,total):
    """ Function to get the ensemble average of the velocity auto-correlation
    for a given run, averaged over multiple runs. File name has to be dump.vel*.
    Input
    corlen - correlation length (last 3 must be vx,vy,vz)
    num_atoms - number of atoms in dump file
    nfield - number of columns in dump file
    total - total number of ouput configurations in dump file
    Output:
    TVACF_v - Velocity auto-correlation averaged """

    #Velocity AutoCorr sums
    TVACF_v = np.zeros(corlen,dtype=float)
    VACF_vx = np.zeros(corlen,dtype=float)
    VACF_vy = np.zeros(corlen,dtype=float)
    VACF_vz = np.zeros(corlen,dtype=float)

    #Loop over different runs
    flist = glob.glob('d2v.txt')
    for ff in flist:
        f = ff # 'dump.vel'
        data = read_dump(f,num_atoms,nfield,total)
        print("Processing file: ", f)
        VACF_vx.fill(0.000)
        VACF_vy.fill(0.000)
        VACF_vz.fill(0.000)
        
        #Ensemble average of correlation
        blocks = range(corlen,total,corlen)
        for t in blocks:
            for i in range(num_atoms):
                VACF_vx += autocorr(data[i,5,t-corlen:t])
                VACF_vy += autocorr(data[i,6,t-corlen:t])
                VACF_vz += autocorr(data[i,7,t-corlen:t])

        #Average of Vx,Vy,Vz and blocks 
        TVACF_v += ((VACF_vx + VACF_vy + VACF_vz) / 3.00)
        TVACF_v /= len(blocks)
    
    #Avg over flist simulations    
    TVACF_v /= len(flist)
    #Not needed - Normalize with atoms
    #TVACF_v /= num_atoms
    #Normalize v(0)*v(0)
    TVACF_v /= TVACF_v[0]

    x = np.arange(corlen)*dt*(1.0E12) #ps
    #Write VACF to file
    f = open('VACF.dat','w')
    f.write('#Time(ps) Norm. I. (No Units) \n')
    for i in range(x.size):
        f.write('%f %f \n' %(x[i],TVACF_v[i]))
    f.close()

    return x,TVACF_v

if __name__ == "__main__":

    #Driver portion change as needed
    dt = 0.001E-12 * 10
    timesteps = 1.0E6
    corlen = 1000
    total = 10000
    num_atoms = 108
    nfield = 8
    plotvacf = 0

    x,TVACF_v = process(corlen,num_atoms,nfield,total)
    
    if plotvacf == 1:
        plt.plot(x,TVACF_v,'k')
        plt.ylabel(r'$\frac{v(t)\cdot v(t)}{v(0)\cdot v(0)}$',fontsize=20)
        plt.xlabel('Time [ps]')
        plt.show()

    # TODO: Still needs work here
    # FFT - Pad the singal with zeros then smooth with windowing
    # function. This seems to be parameter dependent
    # Other option is to make signal symmetric about 0 prior to FFT
    std = corlen/15.5
    vsize = TVACF_v.size  
    np.lib.pad(TVACF_v,(0,vsize*15),'constant',constant_values=(0))
    TVACF_v *= signal.get_window(('gaussian',std),vsize) 
    fft_v = fft_autocorr(TVACF_v,dt) * 1.0E12 #THz^-1     
    freq = np.fft.rfftfreq(vsize, d=dt) / 1.0E12 #THz
    
    plt.plot(freq,np.abs(fft_v))
    plt.xlim((0,20.0))
    xx,locs = plt.xticks()
    ll = ['%2.0f' % a for a in xx]
    plt.xticks(xx, ll)
    plt.xlabel('THz')
    plt.ylabel(r'VDOS [THz $^{-1}$]')
    plt.show()
