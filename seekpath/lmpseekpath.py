from phonopy.structure.atoms import PhonopyAtoms
from phonopy import Phonopy 
import numpy as np
import seekpath
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms
from lammps import lammps

lammps_input=input("filename-")
_lammps_commands_list = open(lammps_input).read().split('\n')
print("commandlist")
def mass_to_symbol(mass, tolerance=5e-1):
    from phonopy.structure.atoms import atom_data

    for element in atom_data:
        if element[3] is not None and abs(mass - element[3]) < tolerance:
            return element[1]

    return 'H' 
def get_structure_from_lammps(command_list):
    """
    Get the crystal structure from lammps input
    :param command_list: LAMMPS input commands in list (one item for line)
    :return: numpy array matrix with forces of atoms [Natoms x 3]
    """
    

    
    lmp = lammps()

    lmp.commands_list(command_list)
    lmp.command('run 0')
	
    na = lmp.get_natoms()
    print(na)
    xlo =lmp.extract_global("boxxlo", 1)
    xhi =lmp.extract_global("boxxhi", 1)
    ylo =lmp.extract_global("boxylo", 1)
    yhi =lmp.extract_global("boxyhi", 1)
    zlo =lmp.extract_global("boxzlo", 1)
    zhi =lmp.extract_global("boxzhi", 1)
    xy =lmp.extract_global("xy", 1)
    yz =lmp.extract_global("yz", 1)
    xz =lmp.extract_global("xz", 1)

    unitcell = np.array([[xhi-xlo, xy,  xz],
                           [0,  yhi-ylo,  yz],
                           [0,   0,  zhi-zlo]]).T
    print(unitcell)                       
    positions = lmp.gather_atoms("x", 1, 3)
#    type_mass = lmp.gather_atoms("mass", 1, 1)
    type_mass = lmp.extract_atom("mass", 2)

    type = lmp.gather_atoms("type", 0, 1)

    positions = np.array([positions[i] for i in range(na * 3)]).reshape((na, 3))
    masses = np.array([type_mass[type[i]] for i in range(na)])
    symbols = [mass_to_symbol(masses[i]) for i in range(na)]

    return PhonopyAtoms(positions=positions,
                        masses=masses,
                        symbols=symbols,
                        cell=unitcell)
_structure = get_structure_from_lammps(_lammps_commands_list)
cell = _structure.get_cell()
positions = _structure.get_scaled_positions()
numbers = np.unique(_structure.get_chemical_symbols(), return_inverse=True)[1]

path_data = seekpath.get_path((cell, positions, numbers))

labels = path_data['point_coords']

band_ranges = []
for set in path_data['path']:
     band_ranges.append([labels[set[0]], labels[set[1]]])
dict= {'ranges': band_ranges , 'labels': path_data['path']}
print(dict)
with open("output.txt","w") as f :
    for i  in range(len(dict['labels'])):
         f.writelines(str(dict['labels'][i])+":"+str(dict['ranges'][i]))
