FC=gfortran

#LAPACK=/NSFCGZ/app/LAPACK/3.5.0/lib
#BLAS=/NSFCGZ/app/BLAS/3.5.0/lib
LAPACK=/usr/lib/x86_64-linux-gnu
BLAS=/usr/lib/x86_64-linux-gnu

exe:
	$(MAKE) PhononUnfolding.out   
clean:
	rm -f *.o *.exe *.mod
clear:
	rm -f *.o *.mod

#-------------------------------------------------------------------------------
PhononUnfoldingModule.o : PhononUnfoldingModule.f90
	$(FC) -c  PhononUnfoldingModule.f90 

PhononUnfolding.o : PhononUnfolding.f90
	$(FC) -c  PhononUnfolding.f90

PhononUnfolding.out:  PhononUnfoldingModule.o PhononUnfolding.o
	$(FC)  PhononUnfoldingModule.o PhononUnfolding.o   $(LAPACK)/liblapack.a $(BLAS)/libblas.a  -o  PhononUnfolding.out

#-------------------------------------------------------------------------------
