!  by Jin-Wu Jiang, jwjiang5918@hotmail.com; jiangjinwu@shu.edu.cn

!  structure for 1T-WTe2   

   program main
   implicit none
   integer, parameter:: i4  = 4
   integer, parameter:: dp  = kind(1.0d0)

   integer(i4), parameter:: natomcell = 12_i4
   integer(i4), parameter:: natype = 3_i4
   integer(i4), parameter:: maxneigh = 6_i4
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
   real(dp):: latconst
   real(dp):: thicknessiner
   real(dp):: thickness
   real(dp):: lenxyz(3)
   real(dp):: mdbox(3)
   real(dp), allocatable:: xalat(:)
   real(dp), allocatable:: yalat(:)
   real(dp), allocatable:: zalat(:)
   real(dp), allocatable:: coord(:,:,:)

   integer(i4):: jxsf
   integer(i4):: jlammps
   integer(i4):: jovito
   integer(i4):: i
   integer(i4):: j
   integer(i4):: k
   integer(i4):: l
   integer(i4):: ii
   integer(i4):: il
   integer(i4):: nx
   integer(i4):: ny
   integer(i4):: nz
   integer(i4):: N
   logical:: ltube
   real(dp):: box_x
   real(dp):: box_y
   real(dp):: a
   real(dp):: b
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

   latconst = 3.4970
   bond = 2.7202

   ! derived parameters
   ! inplane 'bond', mimic graphene                                                                    
   b = latconst/sqrt(3.0_dp)
   ! heighth of X atom
   h = sqrt(bond*bond-b*b)
   thicknessiner = 2.0*h
   ! thickness, arbitrary value, should be revised forfew-layer or bulk system
   thickness = 6.0

   a = latconst
   c = thickness
   ci = thicknessiner
   allocate(coord(3,natomcell,nlayer))
   ! armchair and zigzag are different by the followingcoordinate transformation.
   !            0  1  0
   !           -1  0  0
   !            0  0  1
   ! i.e., x --> y, y--> -x, and z --> z
   ! armchair
   if(orientation.eq.'a')then                                                                          
      do i = 1, nlayer, 1
         x0 = 0.5*b
         y0 = 0.0
         z0 = real(i-1) * thickness
         coord(1,1,i) = 0.0_dp + x0
         coord(2,1,i) = 0.0_dp + y0
         coord(3,1,i) = 0.0_dp + z0
         coord(1,2,i) = -0.5_dp*b + x0
         coord(2,2,i) = 0.5_dp*sqrt(3.0)*b + y0
         coord(3,2,i) = h + z0
         coord(1,3,i) = 0.5*b + x0
         coord(2,3,i) = 0.5*sqrt(3.0)*b + y0
         coord(3,3,i) = -h + z0
         coord(1,4,i) = 0.0_dp + x0
         coord(2,4,i) = a + y0
         coord(3,4,i) = 0.0_dp + z0
         coord(1,5,i) = -0.5_dp*b + x0
         coord(2,5,i) = a + 0.5_dp*sqrt(3.0)*b + y0
         coord(3,5,i) = h + z0
         coord(1,6,i) = 0.5*b + x0
         coord(2,6,i) = a + 0.5*sqrt(3.0)*b + y0
         coord(3,6,i) = -h + z0
         coord(1,7,i) = 1.5*b + x0
         coord(2,7,i) = a + 0.5*sqrt(3.0)*b + y0
         coord(3,7,i) = 0.0 + z0
         coord(1,8,i) = b + x0
         coord(2,8,i) = a + y0
         coord(3,8,i) = h + z0
         coord(1,9,i) = 2.0*b + x0
         coord(2,9,i) = a + y0
         coord(3,9,i) = -h + z0
         coord(1,10,i) = 1.5*b + x0
         coord(2,10,i) = 0.5*sqrt(3.0)*b + y0
         coord(3,10,i) = 0.0 + z0
         coord(1,11,i) = b + x0
         coord(2,11,i) = 0.0 + y0
         coord(3,11,i) = h + z0
         coord(1,12,i) = 2.0*b + x0
         coord(2,12,i) = 0.0 + y0
         coord(3,12,i) = -h + z0
      end do !i
      box_x = 3.0*b
      box_y = 2.0*a
   ! zigzag
   else if(orientation.eq.'z')then                                                                     
      do i = 1, nlayer, 1
         x0 = 0.5*b
         y0 = 0.0
         z0 = real(i-1) * thickness
         coord(1,1,i) = 0.0_dp + y0
         coord(2,1,i) = 0.0_dp + x0
         coord(3,1,i) = 0.0_dp + z0
         coord(1,2,i) = -0.5_dp*sqrt(3.0)*b - y0
         coord(2,2,i) = -0.5_dp*b + x0
         coord(3,2,i) = h + z0
         coord(1,3,i) = -0.5*sqrt(3.0)*b - y0
         coord(2,3,i) = 0.5*b + x0
         coord(3,3,i) = -h + z0
         coord(1,4,i) = -a - y0
         coord(2,4,i) = 0.0_dp + x0
         coord(3,4,i) = 0.0_dp + z0
         coord(1,5,i) = -a - 0.5_dp*sqrt(3.0)*b - y0
         coord(2,5,i) = -0.5_dp*b + x0
         coord(3,5,i) = h + z0
         coord(1,6,i) = -a - 0.5*sqrt(3.0)*b - y0
         coord(2,6,i) = 0.5*b + x0
         coord(3,6,i) = -h + z0
         coord(1,7,i) = -a - 0.5*sqrt(3.0)*b - y0
         coord(2,7,i) = 1.5*b + x0
         coord(3,7,i) = 0.0 + z0
         coord(1,8,i) = -a - y0
         coord(2,8,i) = b + x0
         coord(3,8,i) = h + z0
         coord(1,9,i) = -a - y0
         coord(2,9,i) = 2.0*b + x0
         coord(3,9,i) = -h + z0
         coord(1,10,i) = -0.5*sqrt(3.0)*b - y0
         coord(2,10,i) = 1.5*b + x0
         coord(3,10,i) = 0.0 + z0
         coord(1,11,i) = -0.0 - y0
         coord(2,11,i) = b + x0
         coord(3,11,i) = h + z0
         coord(1,12,i) = -0.0 - y0
         coord(2,12,i) = 2.0*b + x0
         coord(3,12,i) = -h + z0
      end do !i
      box_x = sqrt(3.0_dp)*b * 2.0
      box_y = 3.0_dp*b
   else
      write(*,*)"unknown orientation type",orientation
      stop
   end if

   ! set atom number
   atnumcell(1) =  74
   atnumcell(2) =  52
   atnumcell(3) =  52
   atnumcell(4) =  74
   atnumcell(5) =  52
   atnumcell(6) =  52
   atnumcell(7) =  74
   atnumcell(8) =  52
   atnumcell(9) =  52
   atnumcell(10) =  74
   atnumcell(11) =  52
   atnumcell(12) =  52

   ! set atom type
   atypecell(1) = 1
   atypecell(2) = 2
   atypecell(3) = 3
   atypecell(4) = 1
   atypecell(5) = 2
   atypecell(6) = 3
   atypecell(7) = 1
   atypecell(8) = 2
   atypecell(9) = 3
   atypecell(10) = 1
   atypecell(11) = 2
   atypecell(12) = 3

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
   do i = 1, nx
      do j = 1, ny
         x0 = real(i-1) * box_x
         y0 = real(j-1) * box_y
         z0 = 0.0_dp
         k = (i-1)*ny*natomcell*nlayer + (j-1)*natomcell*nlayer
         do l = 1, nlayer
            do ii = 1, natomcell
               k = k + 1
               xalat(k) = x0 + coord(1,ii,l)
               yalat(k) = y0 + coord(2,ii,l)
               zalat(k) = z0 + coord(3,ii,l)
               atnum(k) = atnumcell(ii)
               atype(k) = atypecell(ii)
            end do !ii
         end do !l
      end do !j
   end do !i
   if(k.ne.N)then
      write(*,*)"inconsistent for atom number",k,N
      stop
   end if

   ! roll up plane into tube, only works for single-layer
   ! the middle layer is purely bent, the inner layer is compressed, 
   ! the outer layer is tensiled
   if(ltube)then
      zav = 0.0
      do i = 1, ntot
         zav = zav + zalat(i)
      end do !i
      zav = zav/real(ntot)

      diameter = lenxyz(2)/pi
      radius0 = 0.5_dp*diameter
      radius1 = 0.5_dp*diameter - 0.5*thicknessiner
      radius2 = 0.5_dp*diameter + 0.5*thicknessiner
      ly = pi * diameter
      do i = 1, ntot, 1
         ! inner layer
         if(zalat(i).lt.(zav-0.1))then
            radius = radius1
         ! outer layer
         else if(zalat(i).gt.(zav+0.1))then
            radius = radius2
         ! middle layer
         else
            radius = radius0
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
