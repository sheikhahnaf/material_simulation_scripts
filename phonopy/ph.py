from phonolammps import Phonolammps
from phonopy import Phonopy
from phonopy import phonon
from phonopy.phonon.band_structure import (get_band_qpoints, get_band_qpoints_by_seekpath)
from phonopy.interface.phonopy_yaml import PhonopyYaml

lammps_filename=input("lammps file name")
phlammps = Phonolammps(lammps_filename,supercell_matrix=[[3, 0, 0],[0, 3, 0],[0, 0, 3]])

unitcell = phlammps.get_unitcell()
force_constants = phlammps.get_force_constants()
supercell_matrix = phlammps.get_supercell_matrix()
npoints=input("npoints")

phonon = Phonopy(unitcell,supercell_matrix)
phonon.set_force_constants(force_constants)
phonon.set_mesh([20, 20, 20])
bands, labels, path_connections = get_band_qpoints_by_seekpath(phonon.primitive, npoints)
phonon.run_band_structure(bands,path_connections=path_connections,labels=labels)
phonon.write_yaml_band_structure(comment=comment)

phonon.set_total_DOS()
phonon.plot_total_DOS().show()

phonon.set_thermal_properties()
phonon.plot_thermal_properties().show()
