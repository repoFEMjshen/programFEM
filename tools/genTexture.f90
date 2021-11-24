      PROGRAM texture
      implicit none
      real(8) :: phi1,theta,phi2
      real(8) :: x,y
      real(8) :: n(3),na1(3),na2(3),na3(3)
      real(8) :: p_tex(3), n_tex(3), tex_intensity, tex_this
      real(8) :: axisRotate(3),AxisCrossN(3),degRotate
      real(8), parameter :: PI=3.14159265359D0
      real(8), allocatable :: grainEuler(:,:),elementGrains(:)
      real(8), allocatable :: grainCAxes(:,:)
      character(len=255) :: strFileNameInitTexture
      integer, allocatable :: grainPhase(:)
      integer :: i, nGrains,nGrains_file, iGrain, nElem, iElem, elemIdx, ierror
      character :: textureForm
      logical :: grainPhasesImported,fileExists,elementGrainListImported
      logical :: doneTransformations

      write(*,*) 'how many grains?'
      read(*,*) nGrains

      write(*,*) 'enter file name to import initial texture.'
      write(*,*) 'enter 0 to start from random texture:'
      read(*,*) strFileNameInitTexture
      
      CALL init_random_seed()
      
      ! read element grains
      elementGrainListImported = .false.
      inquire(file='elementGrains.inp',exist=fileExists)
      if(fileExists) then
         write(*,*) 'reading element grains...'
         open(200,file='elementGrains.inp')
         read(200,*) nElem
         write(*,*) '# of elements:',nElem
         allocate(elementGrains(nElem))
         do elemIdx=1,nElem
            read(200,*,iostat=ierror) iElem, elementGrains(iElem)
         enddo
         close(200)
         if (ierror /= 0) then
            write(*,*) 'file read error: elementGrains.inp'
            stop
         endif
         elementGrainListImported = .true.
      else
         write(*,*) 'WARNING! elementGrains.inp not found.'
         write(*,*) 'Will not create a list of orientations for each element.'
      endif
      
      ! read grain phases
      allocate(grainPhase(nGrains))
      CALL readGrainPhases('grainPhases.inp',nGrains,grainPhase,grainPhasesImported)
      if(grainPhasesImported) write(*,*) 'Grain Phases imported'
      
      ! sample (0001) poles from uniform distribution on S2 sphere
      allocate(grainCAxes(nGrains,3))
      allocate(grainEuler(nGrains,3))
      grainCAxes = 0.d0
      grainEuler = 0.d0
      if (trim(strFileNameInitTexture) /= '0')  then
         inquire(file=trim(strFileNameInitTexture),exist=fileExists)
         if (.not.fileExists) then
            write(*,*) 'file not found:',trim(strFileNameInitTexture)
            stop
         endif
         open(100,file=trim(strFileNameInitTexture))
         read(100,*) nGrains_file
         do i=1,nGrains
            read(100,*) iGrain, grainEuler(iGrain,:)
            CALL getNc_fromEuler(n(1),n(2),n(3),grainEuler(iGrain,1),grainEuler(iGrain,2),grainEuler(iGrain,3))
            grainCAxes(iGrain,:) = n
         enddo
         close(100)
      else
         do iGrain = 1,nGrains
            CALL getRandomDirection(n(1),n(2),n(3),phi1,theta)
            grainCAxes(iGrain,:) = n
            CALL RANDOM_NUMBER(phi2)
            phi2 = phi2*2.D0*PI
            grainEuler(iGrain,3) = phi2
         enddo      
      endif
      
      ! apply transformations: (1) extend/compress (2) rotate
      doneTransformations = .false.
      do while (.not.doneTransformations)
      
         write(*,*) 'extension/compression axis for the transformation of c-axes:'
         read(*,*) p_tex(1:3)
         tex_intensity = DSQRT(DOT_PRODUCT(p_tex,p_tex))
         n_tex = p_tex/tex_intensity   ! texture bundle normal
         
         write(*,*) 'Form texture by (E)xtending or (C)ompressing pole figures?'
         write(*,*) 'extension creates basal textures, compression creates a "ring" of (0001) axes'
         read(*,*) textureForm
         if (textureForm=='e') textureForm = 'E'
         if (textureForm=='c') textureForm = 'C'
         if (textureForm /='C' .and. textureForm /= 'E') then
            write(*,*) 'unrecognized option'
            stop
         endif
         
         write(*,*) 'after extension/compression rotate poles along axis:'
         read(*,*) axisRotate(1:3)
         if (.not.all(axisRotate==0.d0)) &
            axisRotate = axisRotate / DSQRT(DOT_PRODUCT(axisRotate,axisRotate))
         
         write(*,*) 'by degrees:'
         read(*,*) degRotate
         degRotate = degRotate / 180.d0 * PI
         
         do iGrain = 1,nGrains
            n = grainCAxes(iGrain,:)
            ! stretch the direction vectors along p_tex, with a stretch of tex_intensity
            if (tex_intensity.GT.0.000001D0) then ! to avoid 0/0, skip procedure if stretch is close to zero
               if (textureForm=='E') then ! extend along p
                  !n(:) = n(:) + DOT_PRODUCT(p_tex,n)*n_tex(:)
                  n(:) = n(:) - DOT_PRODUCT(n_tex,n)*n_tex(:) + DOT_PRODUCT(n_tex,n)*n_tex(:)*(1.d0 + tex_intensity)
               elseif (textureForm=='C') then ! compress along p
                  n(:) = n(:) - DOT_PRODUCT(n_tex,n)*n_tex(:) + DOT_PRODUCT(n_tex,n)*n_tex(:)*(1.d0/(1.d0 + tex_intensity))
               endif
            endif
            
            ! rotate the axis:
            if (any(axisRotate /= 0.d0) .and. &
                degRotate /= 0.d0) then
               call cross(AxisCrossN,axisRotate,n)
               n = n*cos(degRotate) + AxisCrossN*sin(degRotate) + &
                   axisRotate*dot_product(axisRotate,n)*(1.d0 - cos(degRotate))
            endif

            ! re-normalize the direction
            n(:) = n(:) / DSQRT(DOT_PRODUCT(n,n))
            
            ! save to the list of axes
            grainCAxes(iGrain,:) = n
         enddo
         
         write(*,*) 'finished with transformations? (T/F)'
         read(*,*) doneTransformations
         
      enddo
      
      ! calculate the euler angles, also assign the a-axis (2-1-10) orientations
      do iGrain = 1,nGrains
         n = grainCAxes(iGrain,:)
         ! calculate the two euler angles for the direction
         CALL getPhi1andTheta(n(1),n(2),n(3),phi1,theta)
         ! assign the third one with a uniform random distribution
         !CALL RANDOM_NUMBER(phi2)
         !phi2 = phi2*2.D0*PI

         grainEuler(iGrain,1) = phi1
         grainEuler(iGrain,2) = theta
         !grainEuler(iGrain,3) = phi2
      enddo
      
      ! export the orientations for each element
      if (elementGrainListImported) then
         open(100,file='texture_new.dat')
         do iElem=1,nElem
            iGrain = elementGrains(iElem)
            write(100,*) grainEuler(iGrain,:), grainPhase(iGrain)
         enddo
         close(100)
      endif
      
      ! export the orientations for each grain
      open(100,file='grainTexture_new.inp')
      write(100,*) nGrains
      do iGrain=1,nGrains
         write(100,*) iGrain, grainEuler(iGrain,:)
      enddo
      close(100)
      
      ! export pole figure and spherical plots for HCP
      open(205,file='pole_figures.dat')
      write(205,*) nGrains
      open(101,file='debug_c.dat')
      write(101,'(A)',advance='no') 'ListPointPlot3D[{'
      ! c-axes
      do iGrain=1,nGrains
         CALL getNc_fromEuler(n(1),n(2),n(3),grainEuler(iGrain,1),grainEuler(iGrain,2),grainEuler(iGrain,3))
         
         write(101,'(A,F6.2,A,F6.2,A,F6.2,A)',advance='no') '{',n(1),',',n(2),',',n(3),'},'

         if(n(3).LT.0.0) n(:) = -n(:)
         CALL getStereographicProjection(n,x,y)
         write(205,*) x,y
      enddo
      write(101,'(A,F6.2,A,F6.2,A,F6.2,A)') '{',0.0,',',0.0,',',0.0,'}},BoxRatios->{1,1,1},PlotRange->{{-1,1},{-1,1},{-1,1}}]'
      close(101)

      write(205,*) 3*nGrains
      open(101,file='debug_a1.dat')
      open(102,file='debug_a2.dat')
      open(103,file='debug_a3.dat')
      write(101,'(A)',advance='no') 'ListPointPlot3D[{'
      write(102,'(A)',advance='no') 'ListPointPlot3D[{'
      write(103,'(A)',advance='no') 'ListPointPlot3D[{'
      ! a-axes      
      do iGrain=1,nGrains
         CALL getNas_fromEuler(na1,na2,na3,grainEuler(iGrain,1),grainEuler(iGrain,2),grainEuler(iGrain,3))         
         !if(na1(3).LT.0.0) na1(:) = -na1(:)
         !if(na2(3).LT.0.0) na2(:) = -na2(:)
         !if(na3(3).LT.0.0) na3(:) = -na3(:)
         write(101,'(A,F6.2,A,F6.2,A,F6.2,A)',advance='no') '{',na1(1),',',na1(2),',',na1(3),'},'
         write(102,'(A,F6.2,A,F6.2,A,F6.2,A)',advance='no') '{',na2(1),',',na2(2),',',na2(3),'},'
         write(103,'(A,F6.2,A,F6.2,A,F6.2,A)',advance='no') '{',na3(1),',',na3(2),',',na3(3),'},'

         na1(3) = na1(3) + 1.D0
         na2(3) = na2(3) + 1.D0
         na3(3) = na3(3) + 1.D0
         if(na1(3)>1.0) write(205,*) na1(1)/na1(3), na1(2)/na1(3)
         if(na2(3)>1.0) write(205,*) na2(1)/na2(3), na2(2)/na2(3)
         if(na3(3)>1.0) write(205,*) na3(1)/na3(3), na3(2)/na3(3)
      enddo
      close(205)
      write(101,'(A,F6.2,A,F6.2,A,F6.2,A)') '{',0.0,',',0.0,',',0.0,'}},BoxRatios->{1,1,1},PlotRange->{{-1,1},{-1,1},{-1,1}}]'
      write(102,'(A,F6.2,A,F6.2,A,F6.2,A)') '{',0.0,',',0.0,',',0.0,'}},BoxRatios->{1,1,1},PlotRange->{{-1,1},{-1,1},{-1,1}}]'
      write(103,'(A,F6.2,A,F6.2,A,F6.2,A)') '{',0.0,',',0.0,',',0.0,'}},BoxRatios->{1,1,1},PlotRange->{{-1,1},{-1,1},{-1,1}}]'
      
      close(101)
      close(102)
      close(103)
      
      write(*,*) 'texture created for', nGrains, 'grains.'

      
      if (allocated(elementGrains)) deallocate(elementGrains)
      deallocate(grainEuler)
      deallocate(grainPhase)
      
      END PROGRAM
      
      ! verified to produce a uniform distribution of directions
      SUBROUTINE getRandomDirection(nx,ny,nz,phi1,theta)
      implicit none
      real(8), intent(out) :: nx,ny,nz
      real(8), intent(out) :: phi1,theta
      real(8), parameter :: PI=3.14159265359D0
      
      real(8) :: u,v
      
      CALL RANDOM_NUMBER(u)
      CALL RANDOM_NUMBER(v)
      
      phi1 = 2.D0*PI*u
      theta = acos(2*v-1)
      
      nx = sin(theta)*sin(phi1)
      ny = -sin(theta)*cos(phi1)
      nz = cos(theta)

      END SUBROUTINE
      
      SUBROUTINE getPhi1andTheta(nx,ny,nz,phi1,theta)
      implicit none
      real(8), intent(in) :: nx,ny,nz
      real(8), intent(out) :: phi1,theta
      real(8), parameter :: PI=3.14159265359D0

      if (nx==0.d0 .and. ny==0.d0 .and. nz==0.d0) then
         phi1 = 0.d0
         theta = 0.d0
         return
      endif
      
      theta = acos(nz)
      if (nx==0.d0 .and. ny==0.d0) then   ! singularity here
         phi1=0.d0
      else
         phi1 = acos(-ny/DSQRT(nx*nx+ny*ny))
      endif
      if(nx.LT.0.D0) phi1 = 2.D0*PI-phi1

      END SUBROUTINE

      SUBROUTINE getNc_fromEuler(nx,ny,nz,phi1,theta,phi2)
      implicit none
      real(8), intent(out) :: nx,ny,nz
      real(8), intent(in) :: phi1,theta,phi2
      real(8), parameter :: PI=3.14159265359D0

      nx = sin(theta)*sin(phi1)
      ny = -sin(theta)*cos(phi1)
      nz = cos(theta)

      END SUBROUTINE
      
      SUBROUTINE getStereographicProjection(n,PF_x,PF_y)
      implicit none
      real(8), intent(in) :: n(3)
      real(8), intent(out):: PF_x, PF_y
      real(8) :: temp_n(3)
      
      temp_n = n
      if (temp_n(3) < 0.d0) temp_n = -temp_n
      temp_n(3) = temp_n(3) + 1.D0
      PF_x = temp_n(1) / temp_n(3)
      PF_y = temp_n(2) / temp_n(3)
      
      END SUBROUTINE
      
      SUBROUTINE getNas_fromEuler(na1,na2,na3,phi1,theta,phi2)
      implicit none
      real(8), intent(out) :: na1(3),na2(3),na3(3)
      real(8), intent(in) :: phi1,theta,phi2
      real(8), parameter :: PI=3.14159265359D0
      real(8) :: n0(3), rot(3,3)
      
      !CALL getRot(rot,phi1,theta,phi2)
      CALL euler_slip(phi1,theta,phi2,rot)
      
      n0(1) = 1.D0
      n0(2) = 0.D0
      n0(3) = 0.D0      
      na1(:) = MATMUL(rot(:,:),n0(:))
      n0(1) = cos(2.D0/3.D0*PI)
      n0(2) = sin(2.D0/3.D0*PI)
      n0(3) = 0.D0
      na2(:) = MATMUL(rot(:,:),n0(:))
      n0(1) = cos(4.D0/3.D0*PI)
      n0(2) = sin(4.D0/3.D0*PI)
      n0(3) = 0.D0
      na3(:) = MATMUL(rot(:,:),n0(:))

      END SUBROUTINE
      
      SUBROUTINE getRot(rot,phi1,theta,phi2)
      implicit none
      real(8), intent(out) :: rot(3,3)
      real(8), intent(in) ::  phi1,theta,phi2
      
      rot(1,1) = cos(theta)
      rot(1,2) = -cos(phi2)*sin(theta)
      rot(1,3) = -sin(phi2)*sin(theta)
      rot(2,1) = cos(phi1)*sin(theta)
      rot(2,2) = cos(phi1)*cos(theta)*cos(phi2)-sin(phi1)*sin(phi2)
      rot(2,3) = -cos(phi2)*sin(phi1)-cos(phi1)*cos(theta)*sin(phi2)
      rot(3,1) = sin(phi1)*sin(theta)
      rot(3,2) = cos(phi1)*sin(phi2)+cos(theta)*cos(phi2)*sin(phi1)
      rot(3,3) = cos(phi1)*cos(phi2)-cos(theta)*sin(phi1)*sin(phi2)
      
      END SUBROUTINE
      
      subroutine euler_slip(phi,theta,omega,tlgt)
      implicit double precision (a-h,o-z)  

!      include 'aba_param.inc'

      dimension tlg(3,3),tlgt(3,3)
      real(8) :: phi,theta,omega

      pi=4.d0*datan(1.d0)

  
      sp=dsin(phi)                      
      cp=dcos(phi)                     
      st=dsin(theta)                     
      ct=dcos(theta)                    
      so=dsin(omega)                    
      co=dcos(omega)   
      tlg(1,1)=co*cp-so*sp*ct
      tlg(1,2)=co*sp+so*ct*cp   
      tlg(1,3)=so*st   
      tlg(2,1)=-so*cp-sp*co*ct 
      tlg(2,2)=-so*sp+ct*co*cp
      tlg(2,3)=co*st
      tlg(3,1)=sp*st       
      tlg(3,2)=-st*cp       
      tlg(3,3)=ct

      tlgt = TRANSPOSE(tlg)

      return
      end   
      
      SUBROUTINE readGrainPhases(strGrainFile,nGrains,grainPhase,grainPhasesImported)
      implicit none
      character(len=*), intent(in) :: strGrainFile
      integer, intent(in) :: nGrains
      integer, intent(out):: grainPhase(nGrains)
      logical, intent(out):: grainPhasesImported

      logical :: fileExists
      integer :: iGrain, nEntries, I

      ! assign default grain phase
      grainPhase(:) = 1
      grainPhasesImported = .FALSE.

      inquire(file=strGrainFile,exist=fileExists)
      if (fileExists) then
         open(105,file=strGrainFile)
         read(105,*) nEntries
         do I=1,nEntries
            read(105,*) iGrain,grainPhase(iGrain)
            if(iGrain.GT.nGrains)then
               write(*,*) 'unidentified grain ID ',iGrain, & 
                          ' in',strGrainFile
               write(*,*) '# of grains in model:',nGrains
               STOP
            endif
         enddo
         close(105)
      endif
      
      if(fileExists) grainPhasesImported = .TRUE.
      
      END SUBROUTINE

      SUBROUTINE cross(vc, v1, v2)
      implicit none
      real(8), intent(in) :: v1(3),v2(3)
      real(8), intent(out):: vc(3)
      
      vc(1)=+v1(2)*v2(3)-v1(3)*v2(2)
      vc(2)=-v1(1)*v2(3)+v1(3)*v2(1)
      vc(3)=+v1(1)*v2(2)-v1(2)*v2(1)
      
      END SUBROUTINE
