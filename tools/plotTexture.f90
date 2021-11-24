      PROGRAM plotTexture

      !use meshtools
      use crystal_Titanium
      use options
      use grains
      
      implicit none
      
      include 'PARDIS.H'
      
      integer :: NX, NELX
      real(8) :: G0XYZ(MAXCRD)
      integer :: IJK(MNELX)
      character(len=30) :: meshName
      logical :: solutionExists

      integer, parameter :: nMaxFilesToProcess = 500
      integer :: iFile, nFilesToProcess
      character(len=128) :: listTextureFilesToProcess(nMaxFilesToProcess)
      character(len=128) :: listTecPlotFileNames(nMaxFilesToProcess)
      character(len=128) :: strTextureFileName,strTecPlotFilename,strCommandLine
      character(len=10) :: strStaTextureNumber,strEndTextureNumber
      integer :: staTextureNumber,endTextureNumber
      real(8) :: phi1,theta,phi2
      real(8) :: xy_c(2,MAXN_GRAINS)
      real(8) :: xy_a1(2,MAXN_GRAINS)
      real(8) :: xy_a2(2,MAXN_GRAINS)
      real(8) :: xy_a3(2,MAXN_GRAINS)
      real(8) :: nc(3),na1(3),na2(3),na3(3)
      real(8) :: p_tex(3), tex_magn, tex_this
      real(8), parameter :: PI=3.14159265359D0
      
      ! slip systems/families
      real(8) :: gmt(3,nslip_max),gst(3,nslip_max)    ! gmt, gst - slip planes and slip directions in the global frame

      integer :: islip
      
      integer :: nGrainsInFile
      real(8) :: grainRot(9,MAXN_GRAINS)
      real(8) :: R_cry(3,3), det
      integer :: i, iGrain
      integer :: error, istep, nstep, step_inp
      integer :: lastProgress
      real(8) :: time_inp, temp_inp
      integer :: lenTextureFileName
      logical :: multipleFiles
      
      integer, parameter :: maxNWords = 10
      integer :: wPos(maxNWords,2), iWord, nWords
      logical :: ignoreVariableNames
      
      ! READ IN THE MESH, CONNECTIVITY
      !CALL READGM(G0XYZ,IJK,NX,NELX,meshName)      

      !write(*,*) 'Mesh imported:'
      !write(*,*) '# of elements: ', NELX
      !write(*,*) '# of nodes: ', NX
      
      multipleFiles = .false.
      
      if (iargc().GE.1) then                                   ! User specified texture file
         CALL getarg(1, strCommandLine)
         if (strCommandLine(1:5)=='list:') then
            CALL parseString(trim(strCommandLine),':',wPos(:,:),maxNWords,nWords)
            strStaTextureNumber = strCommandLine(wPos(2,1):wPos(2,2))   ! start texture number
            strEndTextureNumber = strCommandLine(wPos(3,1):wPos(3,2))   ! end texture number
            if (nWords==3) then
               read(strStaTextureNumber,'(I10)') staTextureNumber
               read(strEndTextureNumber,'(I10)') endTextureNumber
               nFilesToProcess = endTextureNumber-staTextureNumber+1
               
               do iFile = 1,nFilesToProcess
                  write(listTextureFilesToProcess(iFile),'(A,I0,A)') 'grainTexture_',iFile,'.inp'
                  write(listTecPlotFileNames(iFile),'(A,I0,A)') 'pole_figures_init_',iFile,'.plt'
               enddo
               
               multipleFiles = .true.
            endif
         endif
      endif      
      
      if (.not.multipleFiles) then
         if (iargc().GE.1) then
            CALL getarg(1, strTextureFileName)
            lenTextureFileName = len(trim(strTextureFileName))
            if (strTextureFileName(lenTextureFileName-3:lenTextureFileName)/='.inp') then
               strTextureFileName = trim(strTextureFileName)//'.inp'
               lenTextureFileName = len(trim(strTextureFileName))
            endif
            strTecPlotFilename = strTextureFileName(1:lenTextureFileName-4)//'.plt'
         else
            strTextureFileName = 'grainTexture.inp'
            strTecPlotFilename = 'pole_figures_init.plt'
         endif
         
         nFilesToProcess = 1
         listTextureFilesToProcess(1) = strTextureFileName
         listTecPlotFileNames(1) = strTecPlotFilename
         
      endif
      
      ! READ TIME STEPPING INFORMATION
      solutionExists = .false.
      
      if (multipleFiles) then
      
         write(*,*) 'Plotting initial pole figures from multiple files...'
         
      else
         
         inquire(file='time.out',exist=solutionExists)
         if (.not.solutionExists) then
            nstep = 0
            write(*,*) 'Cannot find the solution - will plot initial texture only'
         else
            open(102,file='time.out')
            nstep=0
            do
               read(102,*,iostat=error)
               if(error.EQ.-1) EXIT
               nstep=nstep+1
            enddo
            close(102)
         endif

         write(*,*) '# of time steps:', nstep
         
      endif
               
      ! initialize the HCP lattice 
      CALL init_crys_HCP()
      
      do iFile = 1,nFilesToProcess
      
         strTextureFileName = listTextureFilesToProcess(iFile)
         strTecPlotFilename = listTecPlotFileNames(iFile)
         
         !read texture: grainTexture(:)
         CALL readTexture_direct(trim(strTextureFileName),nGrainsInFile,textureImported)
         if(textureImported) write(*,*) 'texture imported'
         if(.NOT.textureImported) then
            write(*,*) 'WARNING ! no texture information ('//trim(strTextureFileName)//') found.'
            stop
         endif
         
         ! plot initial texture on a separate file

         ! open output file
         OPEN(UNIT=900,FILE=strTecPlotFilename)
         
         ! calculate stereographs for all grains and crystallographic planes
         do iGrain = 1,nGrainsInFile
            ! read rotation matrix
            call euler_slip(grainTexture(1,iGrain),grainTexture(2,iGrain),grainTexture(3,iGrain),R_cry)

            ! rotate the lattice
            do islip = 1,nslip_HCP
               gmt(:,islip) = MATMUL(R_cry(:,:),cmt_HCP(:,islip))  ! rotate slip planes
               gst(:,islip) = MATMUL(R_cry(:,:),cst_HCP(:,islip))  ! rotate slip directions
            enddo
            ! calculate the stereographs of the <a> and <c> poles
            nc(:) = gmt(:,1)
            na1(:) = gst(:,1)
            na2(:) = gst(:,2)
            na3(:) = gst(:,3)         
            CALL getStereographicProjection(nc, xy_c(1,iGrain), xy_c(2,iGrain))
            CALL getStereographicProjection(na1,xy_a1(1,iGrain),xy_a1(2,iGrain))
            CALL getStereographicProjection(na2,xy_a2(1,iGrain),xy_a2(2,iGrain))
            CALL getStereographicProjection(na3,xy_a3(1,iGrain),xy_a3(2,iGrain))
         enddo
         
         ! print zone header for c-axis pole figure
         write(900,*) 'VARIABLES = "X", "Y"'
         write(900,*) 'GEOMETRY'
         write(900,*) 'F=POINT, CS=GRID, X=0,Y=0,Z=0, C=BLACK, S=LOCAL, L=SOLID,'
         write(900,*) 'PL=2, LT=1, CLIPPING=CLIPTOVIEWPORT, DRAWORDER=AFTERDATA, MFC="", EP=72, T=CIRCLE 1'
         write(900,'(A,I0,A)') 'ZONE T="{0001}" C=BLACK, I=',nGrainsInFile,', DATAPACKING=POINT'
         do iGrain=1,nGrainsInFile
            write(900,*) xy_c(1:2,iGrain)
         enddo
         ! print zone header for a-axis pole figure
         write(900,'(A,I0,A)') 'ZONE T="(11-20)" C=RED, I=',3*nGrainsInFile,', DATAPACKING=POINT'
         do iGrain=1,nGrainsInFile
            write(900,*) xy_a1(1:2,iGrain)
            write(900,*) xy_a2(1:2,iGrain)
            write(900,*) xy_a3(1:2,iGrain)
         enddo
         
         ! print the geometric center of the poles.
         ! if the symmetry axis of the material aligns with the global frame, these numbers must come out very small:
         !open(unit=901,file='skew_'//strTextureFileName(1:lenTextureFileName-4)//'.log')
         !write(901,*) sum(xy_c(1,1:nGrainsInFile))/nGrainsInFile, sum(xy_c(2,1:nGrainsInFile))/nGrainsInFile, &
         !             sum(xy_c(1,1:nGrainsInFile))/nGrainsInFile/(maxval(xy_c(1,1:nGrainsInFile))-minval(xy_c(1,1:nGrainsInFile))), &
         !             sum(xy_c(2,1:nGrainsInFile))/nGrainsInFile/(maxval(xy_c(2,1:nGrainsInFile))-minval(xy_c(2,1:nGrainsInFile))) 
         !close(901)
         
         write(*,*) 'plotted initial pole figures for:', nGrainsInFile, ' grains in ',trim(strTextureFileName)

         ! close input & output files
         close(900)
         
      enddo
      
      if (.not.solutionExists) goto 700
      
      ! open grain rotation output file
      open(980,FILE='grain_Rot.out')
      
      ! open output file
      OPEN(UNIT=900,FILE='pole_figures_TP.plt')
                  
      write(*,*) 'Processing...'
      lastProgress=0
      do i=1,MIN(60,nstep)
         write(*,'(A)',advance='no') '_'
      enddo
      write(*,*) ''

      do istep=1,nstep
      
         read(980,*,iostat=error) step_inp, time_inp, temp_inp
         if(error.NE.0) then ! end of file
            exit
         endif
         ! READ TEXTURE
         do iGrain = 1,nGrainsInFile
            read(980,*,iostat=error) grainRot(1:9,iGrain)
            if(error.NE.0) then
               exit
            endif
         enddo
         if(error.NE.0) then
            write(*,*) 'unexpected end of file: grain_AVG.out' 
            write(*,*) '# of grains read: ',iGrain,'at time step',step_inp
            exit
         endif
         
         ! calculate stereographs for all grains and crystallographic planes
         do iGrain = 1,nGrainsInFile
            ! read rotation matrix
            R_cry(1,1:3) = grainRot(1:3,iGrain)
            R_cry(2,1:3) = grainRot(4:6,iGrain)
            R_cry(3,1:3) = grainRot(7:9,iGrain)
            call calcDet(R_cry,det)
            R_cry = R_cry / det**(1.d0/3.d0)
            ! rotate the lattice
            do islip = 1,nslip_HCP
               gmt(:,islip) = MATMUL(R_cry(:,:),cmt_HCP(:,islip))  ! rotate slip planes
               gst(:,islip) = MATMUL(R_cry(:,:),cst_HCP(:,islip))  ! rotate slip directions
            enddo
            ! calculate the stereographs of the <a> and <c> poles
            nc(:) = gmt(:,1)
            na1(:) = gst(:,1)
            na2(:) = gst(:,2)
            na3(:) = gst(:,3)
            CALL getStereographicProjection(nc, xy_c(1,iGrain), xy_c(2,iGrain))
            CALL getStereographicProjection(na1,xy_a1(1,iGrain),xy_a1(2,iGrain))
            CALL getStereographicProjection(na2,xy_a2(1,iGrain),xy_a2(2,iGrain))
            CALL getStereographicProjection(na3,xy_a3(1,iGrain),xy_a3(2,iGrain))
         enddo
         
         ! print zone header for c-axis pole figure
         write(900,*) 'VARIABLES = "X", "Y"'
         write(900,*) 'GEOMETRY'
         write(900,*) 'F=POINT, CS=GRID, X=0,Y=0,Z=0, C=BLACK, S=LOCAL, L=SOLID,'
         write(900,*) 'PL=2, LT=1, CLIPPING=CLIPTOVIEWPORT, DRAWORDER=AFTERDATA, MFC="", EP=72, T=CIRCLE 1'
         write(900,'(A,I0,A)') 'ZONE T="{0001}" C=BLACK, I=',nGrainsInFile,', DATAPACKING=POINT, '
         write(900,'(A,F10.3)') 'STRANDID=1, SOLUTIONTIME=',time_inp
         do iGrain=1,nGrainsInFile
            write(900,*) xy_c(1:2,iGrain)
         enddo
         ! print zone header for a-axis pole figure
         write(900,'(A,I0,A)') 'ZONE T="(11-20)" C=RED, I=',3*nGrainsInFile,', DATAPACKING=POINT, '
         write(900,'(A,F10.3)') 'STRANDID=2, SOLUTIONTIME=',time_inp
         do iGrain=1,nGrainsInFile
            write(900,*) xy_a1(1:2,iGrain)
            write(900,*) xy_a2(1:2,iGrain)
            write(900,*) xy_a3(1:2,iGrain)
         enddo
         
      enddo

      write(*,*) 'also exported texture evolution plot: pole_figures_TP.plt'
      
      
700   CONTINUE
      
      if(solutionExists) write(*,*) max(istep,nstep), ' time steps'

      ! close input & output files
      close(980)
      close(900)
            
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
      
      subroutine euler_slip(phi,theta,omega,tlgt)
      implicit double precision (a-h,o-z)  

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

      subroutine calcDet(A,det)
      implicit none
      real(8), intent(in) :: A(3,3)
      real(8), intent(out):: det
      
      det=(A(1,1)*A(2,2)*A(3,3)-A(1,1)*A(2,3)*A(3,2)-A(2,1)*A(1,2)  &
           *A(3,3)+A(2,1)*A(1,3)*A(3,2)+A(3,1)*A(1,2)*A(2,3)-A(3,1)  &
           *A(1,3)*A(2,2))
      end subroutine
      

!*********************************************************
      subroutine READGM(G0XYZ,IJK,NX,NELX,meshName)
      implicit none
      include 'PARDIS.H'
      real(8), intent(out) :: G0XYZ(MAXCRD)
      integer, intent(out) :: IJK(MNELX)
      integer, intent(out) :: NX, NELX
      character*(*), intent(out) :: meshName
      !locals
      integer :: IJKNEL(MAXELN)
      real(8) :: xc(3)
      integer :: MDIM,NDFADD,NDF,NNODE,NGAUSS
      integer :: error
      character(len=50) :: fileName
      character(len=150) :: lineStr
      !dummies
      integer :: icyc,len,i,ii,j,jj,nnum,iElem,iNode,NEL,NEQ
      
      open(201,file='loadtime.inp')
      call READSTR(201,lineStr,error)
      read(201,*) icyc
      call READSTR(201,lineStr,error)
      read(201,*) meshName
      close(201)
      
      len=len_trim(meshName)
      fileName=meshName(1:len)//'.inp'

      open(unit=LR,file=fileName)

      call READSTR(LR,lineStr,error)
      read(LR,*)MDIM,NDFADD
      NDF=MDIM+NDFADD     
      
      call READSTR(LR,lineStr,error)
      read(LR,*) NX
      NEQ=NDF*NX
      
      IF(NDF.GT.NDOFELX)THEN
         WRITE(*,*)'INSUFFICIENT MEM-NDOFELX'
         STOP
      ENDIF

      !error check - deniz
      IF(NX.GT.MAXNODE)THEN
         write(*,*) 'Increase MAXNODE to ', NX
         WRITE(*,*)'Increase MAXNODE to ', NX
         STOP
      ENDIF
      NEQ=NDF*NX
      IF(NDF.GT.MDOFX)THEN
         WRITE(*,*)'INSUFFICIENT MEMORY-MDOFX'
         STOP
      ENDIF
      IF(NX*MDIM.GT.MAXCRD)THEN
         WRITE(*,*)"INSUFFICIENT MEM -MAXCRD"
         STOP
      ENDIF

      IF(MDIM.GT.MAXDIM)THEN
         WRITE(*,*)'INSUFFICIENT MEM-MAXDIM'
         STOP
      ENDIF
      
      ! READ INITIAL POSITIONS
      do iNode=1,NX
         read(LR,*) nnum,(xc(ii),ii=1,MDIM)
         do ii=1,MDIM
            G0XYZ(MDIM*(nnum-1)+ii)=xc(ii)
         enddo
      enddo
      
      
      call READSTR(LR,lineStr,error)
      if(IBRICKORTET.EQ.1)then
         NNODE=8
         NGAUSS=8
      else
         NGAUSS=1
         NNODE=4
      endif

      read(LR,*) NELX

      ! READ IN THE IJK/CONNECTIVITY MATRIX
      jj=0
      do iElem=1,NELX
         read(LR,*) NEL,(IJKNEL(j),j=1,NNODE)
         do j=1,NNODE
            jj=jj+1
            IJK(jj)=IJKNEL(j)
         enddo
      enddo
      
      close(LR)

      return
      end
      
      SUBROUTINE READSTR(fileNum,lineStr,error)
      IMPLICIT REAL*8(A-H,O-Z)
      ! this code skips 1 line of text in the input file.
      integer, intent(in) :: fileNum
      character(150), intent(out) :: lineStr
      integer, intent(out) :: error
      read(fileNum,'(A)',iostat=error) lineStr
      RETURN
      END
      
      SUBROUTINE cross(vc, v1, v2)
      implicit none
      real(8), intent(in) :: v1(3),v2(3)
      real(8), intent(out):: vc(3)
      
      vc(1)=+v1(2)*v2(3)-v1(3)*v2(2)
      vc(2)=-v1(1)*v2(3)+v1(3)*v2(1)
      vc(3)=+v1(1)*v2(2)-v1(2)*v2(1)
      
      END SUBROUTINE

      SUBROUTINE getStereographicProjection(n_in,x,y)
      implicit none
      real(8), intent(in) :: n_in(3)
      real(8), intent(out) :: x,y
      real(8) :: n(3)
      
      n = n_in
      
      if (n(1)==0.d0 .and. n(2)==0.d0 .and. n(3)==0.d0) then
         x = 0.d0
         y = 0.d0
         return
      endif
      
      if (n(3) < 0.d0) n(:) = - n(:)
      n(:) = n(:) / DSQRT(DOT_PRODUCT(n,n))
      n(3) = n(3) + 1.D0
      x = n(1) / n(3)
      y = n(2) / n(3)
      
      END SUBROUTINE