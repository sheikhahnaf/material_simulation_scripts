!  by Jin-Wu Jiang, jwjiang5918@hotmail.com; jiangjinwu@shu.edu.cn

!  structure for p-SiO     

   program main
   implicit none
   integer, parameter:: i4  = 4
   integer, parameter:: dp  = kind(1.0d0)

   integer, parameter:: pi  = 3.1415926

   character(len=100):: file_tmp
   character(len=1):: orientation
   integer(i4), parameter:: natomcell = 4_i4
   integer(i4), parameter:: natype = 2_i4
   integer(i4), parameter:: maxneigh = 3_i4
   real(dp), parameter:: bigspace = 100.0_dp
   integer(i4):: ntot
   integer(i4):: atnumcell(natomcell)
   integer(i4):: atypecell(natomcell)
   integer(i4), allocatable:: atnum(:)
   integer(i4), allocatable:: atype(:)
   real(dp):: thickness
   real(dp):: lenxyz(3)
   real(dp):: mdbox(3)
   real(dp):: shift0(3)
   real(dp), allocatable:: xalat(:)
   real(dp), allocatable:: yalat(:)
   real(dp), allocatable:: zalat(:)
   real(dp), allocatable:: coord(:,:)
   integer(i4):: jfile
   integer(i4):: jxsf
   integer(i4):: jlammps
   integer(i4):: jovito
   integer(i4):: i
   integer(i4):: j
   integer(i4):: k
   integer(i4):: l
   integer(i4):: ii
   integer(i4):: is
   integer(i4):: nx
   integer(i4):: ny
   integer(i4):: nz
   integer(i4):: N
   integer(i4):: id1
   integer(i4):: id2
   integer(i4):: nlayer
   logical:: ltube
   real(dp):: box_x
   real(dp):: box_y
   real(dp):: x0
   real(dp):: y0
   real(dp):: z0
   real(dp):: zav
   real(dp):: radius
   real(dp):: diameter
   real(dp):: theta
   real(dp):: ly
   real(dp):: xmin
   real(dp):: ymin
   real(dp):: a1
   real(dp):: a2
   real(dp):: a3
   real(dp):: u
   real(dp):: v
   real(dp):: w
   real(dp):: d1
   real(dp):: d2
   real(dp):: t1
   real(dp):: t2
   real(dp):: rtmp(3)
   real(dp):: t123
   real(dp):: r1A
   real(dp):: r34
   real(dp):: t423
   real(dp):: rA4
   real(dp):: t1A4
   real(dp):: xi
   real(dp):: eta
   real(dp):: delta
   real(dp):: alpha
   real(dp):: beta

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

   id1 =   14
   id2 =    8
   a1 =  4.7010 !along x
   a2 =  2.7390 !along y
   ! arbitrary value for a3, should be modified for few-layer or bulk system
   a3 = 10.0
   d1 =  1.8430 !r13
   d2 =  1.8590 !r14
   t1 =  1.6755 !theta134 with apex 1
   t2 =  2.3126 !theta425 with apex 4

   ! derived parameters, see my paper for atom notations
   ! A is at the middle of r23
   ! theta1A4 = ?
   t123 = acos(1.0-0.5*a2*a2/d1/d1)
   r1A = d1*cos(0.5*t123)
   r34 = sqrt(d1*d1+d2*d2-2.0*d1*d2*cos(t1))
   t423 = acos(1.0-0.5*a2*a2/r34/r34)
   rA4 = r34*cos(0.5*t423)
   t1A4 = acos(0.5*(d2*d2+r1A*r1A-rA4*rA4)/d2/r1A)

   ! alpha is the angle between r14 and z-axis, beta is the angle between r1A
   ! and z-axis. They are related by two equations:
   !    alpha + beta = theta1A4   ------------------------(1)
   !    r14*sin(alpha)+r1A*sin(beta)=0.5*a1   ------------(2)
   ! following is the solution of equation (2)
   xi = r1A*sin(t1A4)
   eta = d2-r1A*cos(t1A4)
   delta = 4.0*xi*xi*xi*xi+4.0*eta*eta*xi*xi-xi*xi*a1*a1
   rtmp(1) = 0.5*(eta*a1+sqrt(delta))/(eta*eta+xi*xi)
   rtmp(2) = 0.5*(eta*a1-sqrt(delta))/(eta*eta+xi*xi)
   if(rtmp(1)<rtmp(2))then
      alpha = asin(rtmp(1))
   else
      alpha = asin(rtmp(2))
   end if
   beta = t1A4-alpha

   ! (u, v, w) = ?
   u = 0.5*d2*sin(alpha)/a1
   v = 0.5*d2*cos(alpha)/a3
   w = r1A*cos(beta)/a3

   thickness = 0.5*a3

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
      coord(3,3) = (v+w) * a3
      coord(1,4) = (0.5_dp+u) * a1
      coord(2,4) = 0.5_dp * a2
      coord(3,4) = (-v+w) * a3
      box_x = a1
      box_y = a2
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
      coord(3,3) = (v+w) * a3
      coord(1,4) = 0.5_dp * a2
      coord(2,4) = -(0.5_dp+u) * a1
      coord(3,4) = (-v+w) * a3
      box_x = a2
      box_y = a1
      shift0(1) = 0.5_dp * a2
   else
      write(*,*)"unknown orientation type ",orientation                                                   
      stop
   end if

   ! set atom number
   do i = 1, natomcell, 2
      atnumcell(i) = id1
   end do !i
   do i = 2, natomcell, 2
      atnumcell(i) = id2
   end do !i
   ! set atom type
   do i = 1, natomcell
      atypecell(i) = i
   end do !i

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
      ly = pi * diameter
      do i = 1, ntot, 1
         radius = 0.5*diameter + (zalat(i)-zav)
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
   write(jlammps,*)natomcell," atom types"
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
