!  by Jin-Wu Jiang, jwjiang5918@hotmail.com; jiangjinwu@shu.edu.cn

!  structure for borophene, with 12 atom types.

   program main
   implicit none
   integer, parameter:: i4  = 4
   integer, parameter:: dp  = kind(1.0d0)

   integer(i4), parameter:: natomcell = 12_i4
   integer(i4), parameter:: natype = 12_i4
   integer(i4), parameter:: maxneigh = 6_i4
   integer(i4), parameter:: maxlayer = 100_i4
   real(dp), parameter:: pi = 3.141593_dp
   real(dp), parameter:: bohr = 0.529177_dp
   real(dp), parameter:: bohriv = 1.88973_dp
   real(dp), parameter:: mass = 10.811_dp
   real(dp), parameter:: bigspace = 100.0_dp

   character(len=1):: orientation
   integer(i4):: nlayer
   integer(i4):: ntot
   integer(i4), allocatable:: atnum(:)
   integer(i4), allocatable:: atype(:)
   real(dp):: thicknessiner
   real(dp):: thickness
   real(dp):: lenxyz(3)
   real(dp):: mdbox(3)
   real(dp), allocatable:: xalat(:)
   real(dp), allocatable:: yalat(:)
   real(dp), allocatable:: zalat(:)

   integer(i4):: jxsf
   integer(i4):: jlammps
   integer(i4):: jovito
   integer(i4):: jvibra
   integer(i4):: i
   integer(i4):: j
   integer(i4):: k
   integer(i4):: l
   integer(i4):: il
   integer(i4):: iz
   integer(i4):: ii
   integer(i4):: nx
   integer(i4):: ny
   integer(i4):: nz
   integer(i4):: N
   integer(i4):: Nxy
   logical:: ltube
   real(dp), allocatable:: coord(:,:)
   real(dp):: shift0(3)
   real(dp):: shift(3)
   real(dp):: unitvec(3,3)
   real(dp):: box_x
   real(dp):: box_y
   real(dp):: box_z
   real(dp):: x0
   real(dp):: y0
   real(dp):: z0
   real(dp):: zav
   real(dp):: radius
   real(dp):: radius1
   real(dp):: radius2
   real(dp):: diameter
   real(dp):: theta
   real(dp):: ly

   real(dp):: a1
   real(dp):: a2
   real(dp):: a3
   real(dp):: h

   ! set up parameters for simulation
   ! tube or sheet
   ltube = .false.
   ! orientation, a=armchair, z=zigzag
   orientation = 'a'
   nlayer = 1
   lenxyz(1) = 20.0
   lenxyz(2) = 20.0

   ! for tube, this version only works for single-layer
   if(ltube .and. nlayer.ne.1)then
      write(*,*)'MUST set nlayer = 1 for tube'
      stop
   end if

   ! data from New J. Phys. 18, 073016 (2016)
   a1 = 2.0*2.866_dp
   a2 = 3.0*1.614_dp
   a3 = 20.0_dp
   h = 0.911_dp

   thicknessiner = h
   thickness = 0.5_dp * a3

   ! unit vectors
   unitvec(1,1) = a1*2.0
   unitvec(2,1) = 0.0
   unitvec(3,1) = 0.0
   unitvec(1,2) = 0.0
   unitvec(2,2) = a2*3.0
   unitvec(3,2) = 0.0
   unitvec(1,3) = 0.0
   unitvec(2,3) = 0.0
   unitvec(3,3) = a3

   shift0(:) = 0.0_dp
   allocate(coord(3,natomcell))
   ! atoms in each unit cell, see my note
   ! armchair
   if(orientation.eq.'a')then
      coord(1,1) = 0.0_dp
      coord(2,1) = 0.0_dp
      coord(3,1) = 0.0_dp
      coord(1,2) = 0.25 * a1
      coord(2,2) = (1.0/6.0) * a2
      coord(3,2) = h
      coord(1,3) = 0.5*a1
      coord(2,3) = 0.0_dp
      coord(3,3) = 0.0_dp
      coord(1,4) = 0.75 * a1
      coord(2,4) = (1.0/6.0) * a2
      coord(3,4) = h
      coord(1,5) = 0.0_dp
      coord(2,5) = (1.0/3.0) * a2
      coord(3,5) = 0.0_dp
      coord(1,6) = 0.25 * a1
      coord(2,6) = 0.5 * a2
      coord(3,6) = h
      coord(1,7) = 0.5*a1
      coord(2,7) = (1.0/3.0) * a2
      coord(3,7) = 0.0_dp
      coord(1,8) = 0.75 * a1
      coord(2,8) = 0.5 * a2
      coord(3,8) = h
      coord(1,9) = 0.0
      coord(2,9) = (2.0/3.0) * a2
      coord(3,9) = 0.0
      coord(1,10) = 0.25 * a1
      coord(2,10) = (5.0/6.0) * a2
      coord(3,10) = h
      coord(1,11) = 0.5 * a1
      coord(2,11) = (2.0/3.0) * a2
      coord(3,11) = 0.0
      coord(1,12) = 0.75 * a1
      coord(2,12) = (5.0/6.0) * a2
      coord(3,12) = h
      box_x = a1
      box_y = a2
      box_z = 0.5_dp * a3
      shift0(2) = 0.5_dp * a2
   ! zigzag
   else if(orientation.eq.'z')then
      coord(1,1) = 0.0_dp
      coord(2,1) = 0.0_dp
      coord(3,1) = h
      coord(1,2) = (1.0/6.0) * a2
      coord(2,2) = 0.25 * a1
      coord(3,2) = 0.0
      coord(1,3) = 0.0_dp
      coord(2,3) = 0.5*a1
      coord(3,3) = h
      coord(1,4) = (1.0/6.0) * a2
      coord(2,4) = 0.75 * a1
      coord(3,4) = 0.0
      coord(1,5) = (1.0/3.0) * a2
      coord(2,5) = 0.0_dp
      coord(3,5) = h
      coord(1,6) = 0.5 * a2
      coord(2,6) = 0.25 * a1
      coord(3,6) = 0.0
      coord(1,7) = (1.0/3.0) * a2
      coord(2,7) = 0.5*a1
      coord(3,7) = h
      coord(1,8) = 0.5 * a2
      coord(2,8) = 0.75 * a1
      coord(3,8) = 0.0
      coord(1,9) = (2.0/3.0) * a2
      coord(2,9) = 0.0
      coord(3,9) = h
      coord(1,10) = (5.0/6.0) * a2
      coord(2,10) = 0.25 * a1
      coord(3,10) = 0.0
      coord(1,11) = (2.0/3.0) * a2
      coord(2,11) = 0.5 * a1
      coord(3,11) = h
      coord(1,12) = (5.0/6.0) * a2
      coord(2,12) = 0.75 * a1
      coord(3,12) = 0.0
      box_x = a2
      box_y = a1
      box_z = 0.5_dp * a3
      shift0(1) = 0.5_dp * a2
   else
      write(*,'(''# unknown orientation type'',a2)')orientation
      stop
   end if

   ! reset size
   nx = lenxyz(1)/box_x
   lenxyz(1) = real(nx) * box_x
   ny = lenxyz(2)/box_y
   lenxyz(2) = real(ny) * box_y
   nz = nlayer
   lenxyz(3) = real(nz) * box_z

   ! total atom number
   N = natomcell * nx * ny * nz
   ntot = N
   ! allocate ntot-related arrays
   allocate(xalat(ntot), yalat(ntot), zalat(ntot))
   allocate(atnum(ntot), atype(ntot))

   ! generate structure
   do iz = 1, nz
      z0 = real(iz-1) * box_z
      shift(:) = shift0(:) * real(mod(iz-1,2))
      do i = 1, nx
         do j = 1, ny
            x0 = real(i-1) * box_x
            y0 = real(j-1) * box_y
            k = (iz-1)*nx*ny*natomcell + (i-1)*ny*natomcell + (j-1)*natomcell
            do ii = 1, natomcell
               k = k + 1
               xalat(k) = x0 + coord(1,ii) + shift(1)
               yalat(k) = y0 + coord(2,ii) + shift(2)
               zalat(k) = z0 + coord(3,ii) + shift(3)
               atnum(k) = 5
               atype(k) = ii
            end do !ii
         end do !j
      end do !i
   end do !iz
   if(k.ne.N)then
      write(*,'(''# inconsistent for atom number'', 2i12)')k,N
      stop
   end if

   ! roll up plane into tube, only works for single-layer
   ! inner layer
   zav = 0.0
   do i = 1, ntot
      zav = zav + zalat(i)
   end do !i
   zav = zav/real(ntot)

   diameter = lenxyz(2)/pi
   if(ltube)then
      ly = pi * diameter
      radius1 = 0.5_dp*diameter - 0.5*thicknessiner
      radius2 = 0.5_dp*diameter + 0.5*thicknessiner
      do i = 1, ntot, 1
         ! inner layer
         if(zalat(i).lt.zav)then
            radius = radius1
         ! outer layer
         else
            radius = radius2
         end if
         theta = 2.0_dp*pi*yalat(i)/ly
         xalat(i) = xalat(i)
         yalat(i) = radius * cos(theta)
         zalat(i) = radius * sin(theta)
      end do
   end if

   mdbox(:) = lenxyz(:)
   if(ltube)then
      mdbox(2) = mdbox(2) + bigspace
      mdbox(3) = mdbox(3) + bigspace
   else
      mdbox(3) = mdbox(3) + bigspace
   end if

   ! shift to center of simulation box
   if(ltube)then
      do i = 1, ntot
         yalat(i) = yalat(i) + 0.5_dp * mdbox(2)
         zalat(i) = zalat(i) + 0.5_dp * mdbox(3)
      end do !i
   else
      do i = 1, ntot
         zalat(i) = zalat(i) + 0.5_dp * mdbox(3)
      end do !i
   end if

   ! shift for small value, to avoid possible boundary effects
   xalat(:) = xalat(:)+0.3
   yalat(:) = yalat(:)+0.3
   zalat(:) = zalat(:)+0.3

   ! output structure information
   write(*,'(''# ntot = '', i12)')ntot
   if(ltube)then
      write(*,'(''# system = Boro tube, length, diameter = '',2f16.4)')lenxyz(1),diameter
   else
      write(*,'(''# system = Boro sheet, lenxyz(1:3) = '',3f16.4 )')lenxyz(1:3)
   end if
   write(*,'(''# mdbox(1:3) = '',3f16.4 )')mdbox(1:3)
   write(*,'(''# a1, a2, a3 = '',3f16.4 )')a1,a2,a3
   write(*,'(''# nx, ny, nz = '',3i12 )')nx, ny, nz

   ! out .xsf file
   jxsf = 11
   open(jxsf,file='xyz.xsf',status='unknown',form='formatted')
   write(jxsf,'(''# structure data, .xsf formate, viewed by XCRYSDEN'')')
   write(jxsf,'(''# ntot = '',i12)')ntot
   write(jxsf,'(''ATOMS'')')
   do i = 1, ntot
      write(jxsf,'(i12,3f20.8)')atnum(i), xalat(i), yalat(i), zalat(i)
   end do !i
   close(jxsf)

   ! out lammps.dat file
   jlammps = 12
   open(jlammps,file='lammps.dat',status='unknown',form='formatted')
   write(jlammps,'(''# Lammps sturcture data file'')')
   write(jlammps,'(''# ntot = '',i12)')ntot
   write(jlammps,'(i12,''   atoms'')')ntot
   write(jlammps,'(i12,'' atom types'')')natype
   write(jlammps,*)
   write(jlammps,'(2f20.8,''   xlo xhi'')')0.0_dp,mdbox(1)
   write(jlammps,'(2f20.8,''   ylo yhi'')')0.0_dp,mdbox(2)
   write(jlammps,'(2f20.8,''   zlo zhi'')')0.0_dp,mdbox(3)
   write(jlammps,*)
   write(jlammps,'(''Atoms'')')
   write(jlammps,*)
   do i = 1, ntot
      write(jlammps,'(2i12,3f20.8)')i,atype(i),xalat(i),yalat(i),zalat(i)
   end do !i
   close(jlammps)

   ! out ovito file
   jovito = 13
   open(jovito,file='xyz.ovito',status='unknown',form='formatted')
   write(jovito,'(''ITEM: TIMESTEP'')')
   write(jovito,'(i12)')1
   write(jovito,'(''ITEM: NUMBER OF ATOMS'')')
   write(jovito,'(i12)')ntot
   write(jovito,'(''ITEM: BOX BOUNDS pp pp pp'')')
   write(jovito,'(2f20.8,''   xlo xhi'')')0.0_dp,mdbox(1)
   write(jovito,'(2f20.8,''   ylo yhi'')')0.0_dp,mdbox(2)
   write(jovito,'(2f20.8,''   zlo zhi'')')0.0_dp,mdbox(3)
   write(jovito,'(''ITEM: ATOMS id type x y z c_csym'')')
   do i = 1, ntot
      write(jovito,'(2i12,4f20.8)')i,atype(i),xalat(i),yalat(i),zalat(i),0.0
   end do !i
   close(jovito)

   ! out vibra file
   jvibra = 14
   open(jvibra,file='vibra.fdf',status='unknown',form='formatted')
   write(jvibra,'(''# sturcture data for vibra file, should check common lines'')')
   write(jvibra,'(''SystemName siesta'')')
   write(jvibra,'(''SystemLabel siesta'')')
   write(jvibra,'(''NumberOfAtoms '',i12)')natomcell
   write(jvibra,'(''LatticeConstant 1.0 Bohr'')')
   write(jvibra,'(''%block LatticeVectors'')')
   write(jvibra,'(3f8.3)')bohriv*unitvec(1:3,1)
   write(jvibra,'(3f8.3)')bohriv*unitvec(1:3,2)
   write(jvibra,'(3f8.3)')bohriv*unitvec(1:3,3)
   write(jvibra,'(''%endblock LatticeVectors'')')
   write(jvibra,'(''AtomicCoordinatesFormat  NotScaledCartesianBohr'')')
   write(jvibra,'(''%block AtomicCoordinatesAndAtomicSpecies'')')
   do i = 1, natomcell
      write(jvibra,'(3f16.4,i4,f16.4)')bohriv*coord(1:3,i),1,mass
   end do !i
   write(jvibra,'(''%endblock AtomicCoordinatesAndAtomicSpecies'')')
   write(jvibra,'(''SuperCell_1 0'')')
   write(jvibra,'(''SuperCell_2 0'')')
   write(jvibra,'(''SuperCell_3 0'')')
   write(jvibra,'(''AtomicDispl 0.04 Bohr'')')
   write(jvibra,'(''BandLinesScale ReciprocalLatticeVectors'')')
   write(jvibra,'(''%block BandLines'')')
   write(jvibra,'(''1  0.000  0.000  0.000'')')
   write(jvibra,'(''20  0.500  0.000  0.000'')')
   write(jvibra,'(''20  0.500  0.500  0.000'')')
   write(jvibra,'(''20  0.000  0.500  0.000'')')
   write(jvibra,'(''%endblock BandLines'')')
   write(jvibra,'(''#Eigenvectors .true.'')')
   close(jvibra)

   write(*,'(''Job done, Sir!'')')
   end
