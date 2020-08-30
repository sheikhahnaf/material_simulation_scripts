!  by Jin-Wu Jiang, jwjiang5918@hotmail.com; jiangjinwu@shu.edu.cn

!  structure for InS       

   program main
   implicit none
   integer, parameter:: i4  = 4
   integer, parameter:: dp  = kind(1.0d0)

   integer(i4), parameter:: natomcell = 8_i4
   integer(i4), parameter:: natype = 4_i4
   integer(i4), parameter:: maxneigh = 4_i4
   real(dp), parameter:: bigspace = 100.0_dp
   real(dp), parameter:: pi = 3.141593_dp
   character(len=1):: orientation
   integer(i4):: nlayer
   integer(i4):: ntot
   integer(i4):: atnumcell(natomcell)
   integer(i4):: atypecell(natomcell)
   integer(i4), allocatable:: atnum(:)
   integer(i4), allocatable:: atype(:)
   real(dp):: bond
   real(dp):: bondMM
   real(dp):: latconst
   real(dp):: thicknessiner
   real(dp):: thickness
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
   real(dp):: bmm
   real(dp):: c
   real(dp):: ci
   real(dp):: d
   real(dp):: h
   real(dp):: x0
   real(dp):: y0
   real(dp):: z0
   real(dp):: zav
   real(dp):: radius
   real(dp):: radius0
   real(dp):: radius1
   real(dp):: radius2
   real(dp):: radius3
   real(dp):: radius4
   real(dp):: diameter
   real(dp):: theta
   real(dp):: xmin
   real(dp):: ymin
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

   latconst = 3.9400
   bond = 2.5600
   bondMM = 2.8200
   ! derived parameters
   ! inplane 'bond', mimic graphene                                                                    
   b = latconst/sqrt(3.0_dp)
   h = sqrt(bond*bond-b*b)
   ! thickness, arbitrary value, should be revised for few-layer or bulk system
   thickness = 6.0

   a = latconst
   bmm = bondMM
   allocate(coord(3,natomcell))
   ! armchair
   if(orientation.eq.'a')then                                                                          
      coord(1,1) = 0.0
      coord(2,1) = 0.0
      coord(3,1) = 0.0
      coord(1,2) = 0.5*b
      coord(2,2) = 0.5*sqrt(3.0)*b
      coord(3,2) = -h
      coord(1,3) = 1.5*b
      coord(2,3) = 0.5*sqrt(3.0)*b
      coord(3,3) = 0.0
      coord(1,4) = 2.0*b
      coord(2,4) = 0.0
      coord(3,4) = -h
      coord(1,5) = 0.0
      coord(2,5) = 0.0
      coord(3,5) = -2.0*h-bmm
      coord(1,6) = 0.5*b
      coord(2,6) = 0.5*sqrt(3.0)*b
      coord(3,6) = -h-bmm
      coord(1,7) = 1.5*b
      coord(2,7) = 0.5*sqrt(3.0)*b
      coord(3,7) = -2.0*h-bmm
      coord(1,8) = 2.0*b
      coord(2,8) = 0.0
      coord(3,8) = -h-bmm
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
      coord(3,2) = h
      coord(2,3) = 1.5*b
      coord(1,3) = 0.5*sqrt(3.0)*b
      coord(3,3) = -0.0
      coord(2,4) = 2.0*b
      coord(1,4) = 0.0
      coord(3,4) = h
      coord(2,5) = 0.0
      coord(1,5) = 0.0
      coord(3,5) = 2.0*h+bmm
      coord(2,6) = 0.5*b
      coord(1,6) = 0.5*sqrt(3.0)*b
      coord(3,6) = h+bmm
      coord(2,7) = 1.5*b
      coord(1,7) = 0.5*sqrt(3.0)*b
      coord(3,7) = 2.0*h+bmm
      coord(2,8) = 2.0*b
      coord(1,8) = 0.0
      coord(3,8) = h+bmm
      box_y = 3.0*b
      box_x = sqrt(3.0_dp)*b
      shift0(2) = b
      shift0(1) = 0.0
      shift0(3) = 0.0
   else
      write(*,*)"unknown orientation type ",orientation                                                   
      stop
   end if

   ! set atom number
   atnumcell(1) =  16
   atnumcell(2) =  49
   atnumcell(3) =  16
   atnumcell(4) =  49
   atnumcell(5) =  16
   atnumcell(6) =  49
   atnumcell(7) =  16
   atnumcell(8) =  49
   ! set atom type
   atypecell(1) = 1
   atypecell(2) = 2
   atypecell(3) = 1
   atypecell(4) = 2
   atypecell(5) = 3
   atypecell(6) = 4
   atypecell(7) = 3
   atypecell(8) = 4

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
      write(*,*)"inconsistent for atom number",k,N
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
      radius1 = 0.5_dp*diameter - h - 0.5*bmm
      radius2 = 0.5_dp*diameter - 0.5*bmm
      radius3 = 0.5_dp*diameter + 0.5*bmm
      radius4 = 0.5_dp*diameter + 0.5*bmm + h
      ly = pi * diameter
      do i = 1, ntot, 1
         if(zalat(i).gt.(zav+0.5*bmm+0.1))then
            radius = radius1
         else if(zalat(i).gt.(zav+0.1))then
            radius = radius2
         else if(zalat(i).gt.(zav-0.5*bmm-0.1))then
            radius = radius3
         else
            radius = radius4
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

   ! small shift, avoid possible boundary effects
   xmin = 1000.0
   ymin = 1000.0
   do i = 1, ntot
      if(xalat(i).lt.xmin) xmin = xalat(i)                                                                
      if(yalat(i).lt.ymin) ymin = yalat(i)                                                                
   end do !i
   do i = 1, ntot
      xalat(i) = xalat(i) - xmin + 0.1
      yalat(i) = yalat(i) - ymin + 0.1
   end do !i

   ! output structure information
   write(*,*)"# orientation, ntot = ",orientation,ntot
   if(ltube)then
      write(*,*)"# system = tube; length, diameter =",lenxyz(1),diameter
   else
      write(*,*)"# system = sheet; lenxyz(1:3) =",lenxyz(1:3)
   end if
   write(*,*)"# mdbox(1:3) = ",mdbox(1:3)
   write(*,*)"# nx, ny, nz =",nx,ny,nz

   ! out .xsf file
   jxsf = 11
   open(jxsf,file='xyz.xsf',status='unknown',form='formatted')                                         
   write(jxsf,*)"# structure data, .xsf formate, viewed byXCRYSDEN"
   write(jxsf,*)"# ntot =",ntot
   write(jxsf,*)"ATOMS"
   do i = 1, ntot
      write(jxsf,'(i12,3f20.8)')atnum(i), xalat(i), yalat(i), zalat(i)                                    
   end do !i
   close(jxsf)

   ! out lammps.dat file
   jlammps = 12
   open(jlammps,file='lammps.dat',status='unknown',form='formatted')                                   
   write(jlammps,*)"# Lammps sturcture data file"
   write(jlammps,*)"# ntot =",ntot
   write(jlammps,*)ntot,"   atoms"
   write(jlammps,*)natype," atom types"
   write(jlammps,*)
   write(jlammps,*)0.0_dp,mdbox(1),"   xlo xhi"
   write(jlammps,*)0.0_dp,mdbox(2),"   ylo yhi"
   write(jlammps,*)0.0_dp,mdbox(3),"   zlo zhi"
   write(jlammps,*)
   write(jlammps,*)"Atoms"
   write(jlammps,*)
   do i = 1, ntot
      write(jlammps,'(2i12,3f12.4)')i,atype(i),xalat(i),yalat(i),zalat(i)                                 
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

   write(*,*)"Job done, Sir!"
   end
