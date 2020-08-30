from ase.build import mx2
from ase.io.cif import write_cif
from ase.io import read, write
tmdc=input("TMDC name=\n")
slb=mx2(tmdc)
write(tmdc+"_ase.xyz",slb)
