!  by Jin-Wu Jiang, jwjiang5918@hotmail.com; jiangjinwu@shu.edu.cn

!  structure for indiene.

   program main
   implicit none
   integer, parameter:: i4  = 4
   integer, parameter:: dp  = kind(1.0d0)

   integer(i4), parameter:: natomcell = 4_i4
   integer(i4), parameter:: natype = 1_i4
   integer(i4), parameter:: maxneigh = 3_i4
   real(dp), parameter:: bigspace = 100.0_dp
   real(dp), parameter:: pi = 3.141593_dp
   real(dp), parameter:: bohriv = 1.88973_dp

   character(len=1):: orientation
   integer(i4):: nlayer
   integer(i4):: ntot
   integer(i4):: atnumcell(natomcell)
   integer(i4):: atypecell(natomcell)
   integer(i4), allocatable:: atnum(:)
   integer(i4), allocatable:: atype(:)
   real(dp):: bond
   real(dp):: latconst
   real(dp):: thicknessiner
   real(dp):: thickness
   real(dp):: mass
   real(dp):: lenxyz(3)
   real(dp):: mdbox(3)
   real(dp):: shift0(3)
   real(dp), allocatable:: xalat(:)
   real(dp), allocatable:: yalat(:)
   real(dp), allocatable:: zalat(:)
   real(dp), allocatable:: coord(:,:)

   integer(i4):: jxsf
   integer(i4):: jlammps
   integer(i4):: jovito
   integer(i4):: jvibra
   integer(i4):: i
   integer(i4):: j
   integer(i4):: k
   integer(i4):: l
   integer(i4):: ii
   integer(i4):: il
   integer(i4):: is
   integer(i4):: nx
   integer(i4):: ny
   integer(i4):: nz
   integer(i4):: N
   logical:: ltube
   real(dp):: box_x
   real(dp):: box_y
   real(dp):: a
   real(dp):: b
   real(dp):: h
   real(dp):: x0
   real(dp):: y0
   real(dp):: z0
   real(dp):: zav
   real(dp):: radius
   real(dp):: radius0
   real(dp):: radius1
   real(dp):: radius2
   real(dp):: diameter
   real(dp):: theta
   real(dp):: ly

   ! tube or sheet
   ltube = .false.
   ! orientation, a=armchair, z=zigzag
   orientation = 'a'
   nlayer = 1
   lenxyz(1) = 20.0
   lenxyz(2) = 20.0

   ! for tube, this version only works for single-layer
   if(ltube .and. nlayer.ne.1)then
      write(*,*)'MUST set nlayer = 1 for tube, input nlayer = ',nlayer
      stop
   end if

   ! structure parameters from RSC Adv. 6, 8006 (2016)
   bond = 2.89
   latconst = 4.24
   ! thickness value is arbitrarily given, should be modified for few-layer
   ! indiene
   thickness = 5.63

   ! derived structural parameters
   ! inplane 'bond', i.e, projection of bond onto the xy plane
   b = latconst/sqrt(3.0)
   ! height
   h = sqrt(bond*bond-b*b)
   a = latconst

   allocate(coord(3,natomcell))
   ! armchair
   if(orientation.eq.'a')then
      coord(1,1) = 0.0
      coord(2,1) = 0.0
      coord(3,1) = 0.0
      coord(1,2) = 0.5*b
      coord(2,2) = 0.5*sqrt(3.0)*b
      coord(3,2) = h
      coord(1,3) = 1.5*b
      coord(2,3) = 0.5*sqrt(3.0)*b
      coord(3,3) = 0.0
      coord(1,4) = 2.0*b
      coord(2,4) = 0.0
      coord(3,4) = h
      box_x = 3.0*b
      box_y = sqrt(3.0_dp)*b
      shift0(1) = b
      shift0(2) = 0.0
      shift0(3) = 0.0
   ! zigzag
   ! zigzag and armchair are different by the following coordinate
   ! transformation.
   !            0  1  0
   !            1  0  0
   !            0  0 -1
   ! i.e., x --> y, y--> x, and z --> -z
   else if(orientation.eq.'z')then
      coord(2,1) = 0.0
      coord(1,1) = 0.0
      coord(3,1) = -0.0
      coord(2,2) = 0.5*b
      coord(1,2) = 0.5*sqrt(3.0)*b
      coord(3,2) = -h
      coord(2,3) = 1.5*b
      coord(1,3) = 0.5*sqrt(3.0)*b
      coord(3,3) = -0.0
      coord(2,4) = 2.0*b
      coord(1,4) = 0.0
      coord(3,4) = -h
      box_y = 3.0*b
      box_x = sqrt(3.0_dp)*b
      shift0(2) = b
      shift0(1) = 0.0
      shift0(3) = 0.0
   else
      write(*,*)"unknown orientation type ",orientation
      stop
   end if
   atnumcell(:) = 49
   atypecell(:) = 1

   ! reset size
   nx = lenxyz(1)/box_x
   lenxyz(1) = real(nx) * box_x
   ny = lenxyz(2)/box_y
   lenxyz(2) = real(ny) * box_y
   nz = nlayer
   lenxyz(3) = real(nz)*thickness

   ! total atom number
   N = nlayer * natomcell * nx * ny
   ntot = N

   ! allocate ntot-related arrays
   allocate(xalat(ntot), yalat(ntot), zalat(ntot))
   allocate(atnum(ntot), atype(ntot))

   ! generate structure
   do l = 1, nlayer
      z0 = real(l-1)*thickness + shift0(3)
      is = mod(l-1,2)
      do i = 1, nx
         do j = 1, ny
            x0 = real(i-1) * box_x + real(is)*shift0(1)
            y0 = real(j-1) * box_y + real(is)*shift0(2)
            k = (l-1)*nx*ny*natomcell + (i-1)*ny*natomcell + (j-1)*natomcell
            do ii = 1, natomcell
               k = k + 1
               xalat(k) = x0 + coord(1,ii)
               yalat(k) = y0 + coord(2,ii)
               zalat(k) = z0 + coord(3,ii)
               atnum(k) = atnumcell(ii)
               atype(k) = atypecell(ii)
            end do !ii
         end do !j
      end do !i
   end do !l
   if(k.ne.N)then
      write(*,'(''inconsistent for atom number'',2i12)')k,N
      stop
   end if

   ! roll up plane into tube
   if(ltube)then
      zav = 0.0
      do i = 1, ntot
         zav = zav + zalat(i)
      end do !i
      zav = zav/real(ntot)

      diameter = lenxyz(2)/pi
      radius1 = 0.5_dp*diameter - 0.5*h
      radius2 = 0.5_dp*diameter + 0.5*h
      ly = pi * diameter
      do i = 1, ntot, 1
         ! inner layer
         if(zalat(i).lt.(zav-0.1))then
            radius = radius1
         ! outer layer
         else
            radius = radius2
         end if
         theta = 2.0_dp*pi*yalat(i)/ly
         xalat(i) = xalat(i)
         yalat(i) = radius * cos(theta)
         zalat(i) = radius * sin(theta)
      end do !i
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

   ! small shift, to avoid possible boundary effects
   xalat(:) = xalat(:) + 0.1
   yalat(:) = yalat(:) + 0.1
   zalat(:) = zalat(:) + 0.1

   ! output structure information
   write(*,'(''# ntot = '', i12)')ntot
   if(ltube)then
      write(*,'(''# system = MoS2 tube, length, diameter = '',2f16.4)')lenxyz(1),diameter
   else
      write(*,'(''# system = MoS2 sheet, lenxyz(1:3) = '',3f16.4 )')lenxyz(1:3)
   end if
   write(*,'(''# mdbox(1:3) = '',3f16.4 )')mdbox(1:3)
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
   mass = 30.97
   jvibra = 14
   open(jvibra,file='vibra.fdf',status='unknown',form='formatted')
   write(jvibra,'(''# sturcture data for vibra file'')')
   write(jvibra,'(''SystemLabel siesta'')')
   write(jvibra,'(''NumberOfAtoms '',i12)')natomcell
   write(jvibra,'(''LatticeConstant 1.0 Bohr'')')
   write(jvibra,'(''%block LatticeVectors'')')
   write(jvibra,'(3f8.3)')bohriv*lenxyz(1), 0.0, 0.0
   write(jvibra,'(3f8.3)')0.0, bohriv*lenxyz(2), 0.0
   write(jvibra,'(3f8.3)')0.0, 0.0, bohriv*lenxyz(3)
   write(jvibra,'(''%endblock LatticeVectors'')')
   write(jvibra,'(''AtomicCoordinatesFormat  NotScaledCartesianBohr'')')
   write(jvibra,'(''%block AtomicCoordinatesAndAtomicSpecies'')')
   do i = 1, ntot
      write(jvibra,'(3f16.4,i4,f16.4)')bohriv*xalat(i),bohriv*yalat(i),bohriv*zalat(i),1,mass
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
