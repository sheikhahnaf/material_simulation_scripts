#!/usr/bin/python
import sys,string
import numpy as np
from funcs import *
import time

start_time = time.time()

# Read input file from twist.inp
a1,a2,a3,angle,alat,t_z,l2_from_file,l2_file,SL_a1,\
SL_a2, plot_lattice = read_input("twist.inp")

print("----------------------TWISTER---------------------\n")
print("If you use this package please cite the following\n\
paper for which this code was developed:\n\
Phys. Rev. Lett. 121, 266401 (2018)")
print("\nInputs read from twist.inp:")
print("Unit cell vectors (normalized):")
print(a1, a2, a3)
print("Twist angle: %f radians"%(angle))
print("Lattice parameters: %f %f %f"%(alat[0],alat[1],alat[2]))

print("Translation of layer 2 in z-direction: %f Angstroms"%(t_z))
print("Layer 2 basis read from file: %s"%(l2_from_file))

if l2_from_file:
  print("Layer 2 basis read from: %s"%(l2_file))

print("Plot the superlattice: %s"%(plot_lattice))

# The superlattice vectors: 
a1_ang = a1*alat[0]
a2_ang = a2*alat[1]
a1_n = SL_a1[0]*a1_ang + SL_a1[1]*a2_ang
a2_n = SL_a2[0]*a1_ang + SL_a2[1]*a2_ang
print("Superlattice vectors:")
print(a1_n,a2_n)

# Estimate the supercell size required:
tmp = a1_n + a2_n
n_uc = int(np.linalg.norm(tmp)/np.amin(alat[:2]))
n_uc = int(n_uc*2.8)

# create a n_uc x n_uc x 1 supercell.
sc = [n_uc,n_uc,1]

# Place axis of rotation to the center of the 
# generated super cell.
axis = int(n_uc/2)*a1_ang + int(n_uc/2)*a2_ang
print("Placing axis at:")
print(axis)

# Generate atom positions in the super cell for layer 1
pos_ang_l1,labels1,nat_l1,alat_sc,A_ang,n_basis = gen_pos("basis_pos_crys",alat,a1,a2,a3,sc)

layer1 = pos_ang_l1

# Generate atom positions in the super cell for layer 2
if l2_from_file:
  pos_ang_l2,labels2,nat_l2,alat_sc,A_ang,n_basis2 = gen_pos(l2_file,alat,a1,a2,a3,sc)
  layer2 = pos_ang_l2
else:
  layer2 = pos_ang_l1
  labels2 = labels1

layer1_t = np.zeros((nat_l1,3))
layer2_t = np.zeros((nat_l1,3))

# Translate layer 1 in-plane so that rotation axis is at (0,0)
layer1_t[:,0] = layer1[:,0] - axis[0]
layer1_t[:,1] = layer1[:,1] - axis[1]
layer2_t[:,0] = layer2[:,0] - axis[0]
layer2_t[:,1] = layer2[:,1] - axis[1]
layer1_t[:,2] = layer1[:,2]

# Translate layer 2 in out-of-plane direction - by interlayer spacing provided
layer2_t[:,2] = layer2[:,2] + t_z

layer1_t = np.array(layer1_t)
layer2_t = np.array(layer2_t)

# Rotate layer 2
norm = [0.,0.,1.0]
layer2_t = layer2_t
layer2_r = Rotate_atoms(layer2_t,norm,angle)

# PLot the lattice and superlattice vectors
if plot_lattice:
  make_plot(layer1_t,layer2_r,a1,a2,a1_n,a2_n)

# Locate the points inside a1_n, a2_n parallelopiped (superlattice area)
nl1, nl2,V, V_in = \
locate_points_acc(a1_n,a2_n,layer1_t,layer2_r,labels1, labels2,nat_l1,alat, "superlattice.dat",alat[0])

# Check the number of atoms in each layer.
print("Number of atoms in layer 1: %d"%(nl1) )
print("Number of atoms in layer 2: %d"%(nl2) )

# Area of superlattice and the unit cell.
Ar_sl = np.linalg.norm(np.cross(a1_n,a2_n))
Ar_uc = np.linalg.norm(np.cross(a1_ang,a2_ang))
print("Areas, uc: %f, sl: %f"%( Ar_uc, Ar_sl))
print("Expected number of atoms in each layer from area ratios: %d"%( round((Ar_sl/Ar_uc)*n_basis)))

print("Time taken: %12.2f s"%(time.time() - start_time))
