!  by Jin-Wu Jiang, jwjiang5918@hotmail.com; jiangjinwu@shu.edu.cn

!  structure for puckered arsenene (p-arsenene).

   program main
   implicit none
   integer, parameter:: i4 = 4
   integer, parameter:: dp = kind(1.0d0)

   integer(i4), parameter:: natomcell = 4_i4
   integer(i4), parameter:: natype = 2_i4
   integer(i4), parameter:: maxneigh = 3_i4
   integer(i4), parameter:: maxlayer = 100_i4
   real(dp), parameter:: bigspace = 100.0_dp
   real(dp), parameter:: pi = 3.141593_dp

   character(len=1):: orientation
   integer(i4):: nlayer
   integer(i4):: ntot
   integer(i4):: atypecell(natomcell)
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
   integer(i4):: i
   integer(i4):: j
   integer(i4):: k
   integer(i4):: ii
   integer(i4):: l
   integer(i4):: il
   integer(i4):: iz
   integer(i4):: nx
   integer(i4):: ny
   integer(i4):: nz
   integer(i4):: N
   integer(i4):: Nxy
   logical:: ltube
   real(dp), allocatable:: coord(:,:)
   real(dp):: shift0(3)
   real(dp):: shift(3)
   real(dp):: box_x
   real(dp):: box_y
   real(dp):: box_z
   real(dp):: b
   real(dp):: c
   real(dp):: d
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
   real(dp):: u
   real(dp):: v

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

   ! data from first-principles calculations in arxiv:1604.03422. v is
   ! produced by using a3=11.11 from Applied Physics Express 8, 041801 (2015).
   a1 = 4.77_dp
   a2 = 3.68_dp
   a3 = 11.11_dp
   u = 0.0714_dp
   v = 0.108_dp

   thicknessiner = v * a3 * 2.0_dp
   thickness = 0.5_dp * a3

   allocate(coord(3,natomcell))
   shift0(:) = 0.0_dp
   ! atoms in each unit cell, see my note
   ! armchair
   if(orientation.eq.'a')then
      coord(1,1) = -u * a1
      coord(2,1) = 0.0_dp
      coord(3,1) = -v * a3
      coord(1,2) = u * a1
      coord(2,2) = 0.0_dp
      coord(3,2) = v * a3
      coord(1,3) = (0.5_dp-u) * a1
      coord(2,3) = 0.5_dp * a2
      coord(3,3) = v * a3
      coord(1,4) = (0.5_dp+u) * a1
      coord(2,4) = 0.5_dp * a2
      coord(3,4) = -v * a3
      box_x = a1
      box_y = a2
      box_z = 0.5_dp * a3
      shift0(2) = 0.5_dp * a2
   ! zigzag
   else if(orientation.eq.'z')then
      coord(1,1) = 0.0_dp
      coord(2,1) = u * a1
      coord(3,1) = -v * a3
      coord(1,2) = 0.0_dp
      coord(2,2) = -u * a1
      coord(3,2) = v * a3
      coord(1,3) = 0.5_dp * a2
      coord(2,3) = -(0.5_dp-u) * a1
      coord(3,3) = v * a3
      coord(1,4) = 0.5_dp * a2
      coord(2,4) = -(0.5_dp+u) * a1
      coord(3,4) = -v * a3
      box_x = a2
      box_y = a1
      box_z = 0.5_dp * a3
      shift0(1) = 0.5_dp * a2
   else
      write(*,'(''# unknown orientation type'',a2)')orientation
      stop
   end if
   ! type for atoms in the cell
   atypecell(1) = 1
   atypecell(2) = 2
   atypecell(3) = 2
   atypecell(4) = 1

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
      shift(:) = shift0(:) * real(mod(iz,2))
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
            atnum(k) = 33
            atype(k) = atypecell(ii)
         end do !ii
      end do !j
   end do !i
   end do !iz
   if(k.ne.N)then
      write(*,'(''# inconsistent for atom number'', 2i12)')k,N
      stop
   end if

   ! roll up plane into tube, only works for single-layer BP
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
      do i = 1, ntot
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

   ! small shift, to avoid possible boundary effects
   xalat(:) = xalat(:) + 0.5
   yalat(:) = yalat(:) + 0.5
   zalat(:) = zalat(:) + 0.5

   ! output structure information
   write(*,'(''# ntot = '', i12)')ntot
   if(ltube)then
      write(*,'(''# system = tube, length, diameter = '',2f16.4)')lenxyz(1),diameter
   else
      write(*,'(''# system = sheet, lenxyz(1:3) = '',3f16.4 )')lenxyz(1:3)
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

   write(*,'(''Job done, Sir!'')')
   end
