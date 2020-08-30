import numpy as np
import seekpath
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms
from phonopy.interface.calculator import read_crystal_structure
unitcell, _ = read_crystal_structure("POSCAR", interface_mode='vasp')
cell=unitcell.get_cell()
positions = unitcell.get_scaled_positions()
numbers = np.unique(unitcell.get_chemical_symbols(), return_inverse=True)[1]
path_data = seekpath.get_path((cell, positions, numbers))
labels = path_data['point_coords']
band_ranges = []
for set in path_data['path']:
                band_ranges.append([labels[set[0]], labels[set[1]]])

dict= {'ranges': band_range , 'labels': path_data['path']}
with open("output.txt","w") as f :
    for i  in range(len(dict['labels'])):
         f.writelines(str(dict['labels'][i])+":"+str(dict['ranges'][i]))
