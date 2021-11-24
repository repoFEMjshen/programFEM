      PROGRAM findCriticalHardSoftGrains
      
      use meshtools
      use grains
      
      implicit none
      
      include 'PARDIS.H'
      
      real(8), parameter :: PI=3.14159265359D0
      integer :: NX, NELX, NODE
      real(8) :: G0XYZ(MAXCRD)
      integer :: IJK(MNELX), IJKNEL(MAXELN)
      INTEGER:: IELEMNO(MAXNODELEM*MAXNODE),NINDX(MAXNODE)
      
      
      character :: eulerDirection
      real(8), allocatable :: grainRot(:,:,:)
      real(8) :: cAxis(3)
      
      integer :: iGrain, grainID, iGrainID, jGrainID, nGrains_loc
      integer :: SoftGrainID, HardGrainID
      
      integer :: iNeighbor
      integer, allocatable :: nNeighbors(:)
      real(8) :: listN(10)
      integer :: posN(10)
      integer :: iPos
      integer :: nTotNeighbors
      real(8), allocatable :: measureCriticalityList(:)
      integer, allocatable :: grainNeighborsList(:,:)
      integer, allocatable :: listNeighbors(:,:)
      real(8), allocatable :: nyeGrain(:), cauchyGrain(:)
      real(8) :: times(100000)
      logical :: excludeBoundaryGrains
      integer :: exportNucleationMeasureHardGrain,exportNucleationMeasureSoftGrain, exportNucleationMeasureListPos
      real(8) :: cauchy(6)
      integer :: i,j, iList, error
      character(len=150) :: lineStr
      integer :: loadDir
      integer :: iTime
      real(8) :: tTime
      integer :: iStep, nStep
      
      ! variables for domain geometry
      real(8) :: domainRange(3,2), domainLen(3)
      ! element groupings
      integer, allocatable:: elemFaceNeighbors(:,:)
      real(8), allocatable :: elemCenters(:,:)
      real(8), allocatable :: elemVol(:)

      ! grains
      logical :: fileExists, successReadNeighbors,successReadNodeConnect

      character(len=30) :: meshName, strDataFile, strR_type, solnFileName

      ! READ IN THE MESH, CONNECTIVITY
      CALL READGM(G0XYZ,IJK,NX,NELX,meshName)      

      write(*,*) 'Mesh imported:'
      write(*,*) '# of elements: ', NELX
      write(*,*) '# of nodes: ', NX
      
      NODE = 4
      
      write(*,'(A)',advance='no') ' identifying node-to-element connectivity... '
      CALL readNODE_ELEM_CONNEC('nodeElementConnect.dat',IELEMNO,NINDX,NX,successReadNodeConnect)
      if (.not.successReadNodeConnect) then
         write(*,'(A)',advance='no') 'calculating... '
         CALL NODE_ELEM_CONNEC(IELEMNO,NINDX,IJK,NX,NELX,MNELX,MAXNODELEM)
         write(*,*) 'done.'
      else
         write(*,*) 'imported.'
      endif

      ! ---------------------- IMPORT GRAINS ------------------------!
      ! first try grains.inp, if fails, elementGrains.inp, if fails create a trivial grain.
      ! import grains
      ! first try grains.inp, if fails, elementGrains.inp, if fails create a trivial grain.
      CALL readGrains('grains.inp',NX,NELX)
      if(grainsImported) &
         write(*,'(A,I0,A)') ' ', nGrains,' grains imported from grains.inp'   

      if(.NOT.grainsImported) then
         CALL readGrainsFromElementList('elementGrains.inp',NX,NELX)
         if(grainsImported) &
            write(*,'(A,I0,A)') ' ', nGrains,' grains imported from elementGrains.inp'
      endif
      
      if(.NOT.grainsImported) then
         CALL createSingleGrain(NX,NELX)
         write(*,'(A,I0,A)') 'single grain assumed'
      endif
      
      ! read grain texture - Deniz - Sept 22 - readGrainTexture() called from the grains module
      CALL readTexture('grainTexture.inp')
      if(textureImported) write(*,*) 'texture imported'
      if(.NOT.textureImported) then
         write(*,*) 'WARNING ! no texture information (grainTexture.inp) found.'&
      &'            assuming single crystal, euler:',grainTexture(:,1)
      endif         
      
      ! CALCULATE ELEMENT-ELEMENT CONNECTIVITY
      allocate(elemFaceNeighbors(NELX,4))
      write(*,'(A)',advance='no') ' identifying element-element connectivity... '
      CALL readFaceNeighborsTET('elementNeighbors.dat',elemFaceNeighbors,NELX,successReadNeighbors)
      if (.not.successReadNeighbors) then
         CALL calcFaceNeighborsTET(elemFaceNeighbors(:,:),IJK,IELEMNO,NINDX, & 
                                   NELX,NODE,MAXNODELEM,MAXNODE)
         write(*,*) 'done'
      else
         write(*,*) 'imported'
      endif
      if(.not.successReadNeighbors) then
         CALL writeFaceNeighborsTET('elementNeighbors.dat',elemFaceNeighbors,NELX)
      endif

      ! IDENTIFY THE NODES IN INDIVIDUAL GRAINS
      CALL identifyGrainNodes(NX,NELX,IJK)

      ! IDENTIFY THE NODES AT 
      ! 1) DOMAIN BOUNDARY, DOMAIN CORNERS
      ! 2) GRAIN BOUNDARIES
      ! 3) GRAIN-GRAIN BOUNDARIES
      CALL identifyBoundaryNodes(G0XYZ,NELX,NX,IJK,elemFaceNeighbors,domainRange)
      domainLen(:)=domainRange(:,2)-domainRange(:,1)
      
      write(*,*) 'Domain range detected - X', domainRange(1,:)
      write(*,*) 'Domain range detected - Y', domainRange(2,:)
      write(*,*) 'Domain range detected - Z', domainRange(3,:)
      
      !read grain sizes: grainSize(:)
      CALL readGrainSizes('grainSizes.inp')
      if(grainSizesImported) write(*,*) 'grain sizes imported'
      

      loadDir=0
      write(*,*) 
      do while (loadDir < 1 .or. loadDir > 3)
         write(*,*) 'load direction? enter 1, 2 or 3 for X,Y,Z,'
         read(*,*) loadDir
      enddo
      
      write(*,*) ''
      write(*,*) 'exclude soft and hard grains that touch the boundary?'
      read(*,*) excludeBoundaryGrains
      
      write(*,*) ''
      write(*,*) 'also export nucleation measure for a particular hard-soft pair?'
      write(*,*) '(enter Hard and Soft grain IDs, or 0 0 to skip)'
      read(*,*) exportNucleationMeasureHardGrain, exportNucleationMeasureSoftGrain
      
      eulerDirection='A'
      write(*,*) ''
      do while (eulerDirection.NE.'A'.AND.eulerDirection.NE.'P')
         write(*,*) 'Grain orientations are provided in form of Euler angles'
         write(*,*) 'What is the direction of the active transformation associated with the Euler angles?'
         write(*,*) '(A) Euler angles rotate the specimen coordinate axes onto the crystal coordinate axes '&
                   &'(OIM standard, Kourosh, Deniz)'
         write(*,*) '(P) Euler angles rotate the crystal coordinate axes onto the specimen coordinate axes '&
                   &'(was used by Deka)'
         write(*,*) 'What is the direction of the active transformation associated with the Euler angles?'
         read(*,'(A1)') eulerDirection
         if (eulerDirection.EQ.'a') eulerDirection='A'
         if (eulerDirection.EQ.'p') eulerDirection='P'
      enddo
      write(*,*)
      
      open(201,file='grainNeighborList.dat')
      call READSTR(201,lineStr,error)
      read(201,*) nGrains_loc
      allocate(listNeighbors(nGrains_loc,100))
      allocate(nNeighbors(nGrains_loc))
      allocate(cauchyGrain(nGrains_loc),nyeGrain(nGrains_loc))
      listNeighbors = 0
      call READSTR(201,lineStr,error)
      nTotNeighbors = 0
      do iGrain=1,nGrains_loc
         read(201,*) grainID, nNeighbors(grainID), (listNeighbors(grainID,iNeighbor) ,iNeighbor=1,nNeighbors(grainID))
         nTotNeighbors = nTotNeighbors + nNeighbors(grainID)
      enddo
      close(201)
      allocate(measureCriticalityList(nTotNeighbors))
      allocate(grainNeighborsList(nTotNeighbors,2))
      grainNeighborsList = 0
      iList = 0
      do HardGrainID=1,nGrains_loc                               ! <-- hard grain
         do iNeighbor=1,nNeighbors(HardGrainID)
            SoftGrainID = listNeighbors(HardGrainID,iNeighbor)  ! <-- soft grain
            iList = iList + 1
            grainNeighborsList(iList,1) = HardGrainID
            grainNeighborsList(iList,2) = SoftGrainID
            if (exportNucleationMeasureHardGrain == HardGrainID .and. &
                exportNucleationMeasureSoftGrain == SoftGrainID) then
                
               exportNucleationMeasureListPos = iList
                
            endif
         enddo
      enddo
      deallocate(listNeighbors)
            
      !CALL readGrainTexture(grain_rot,grain_euler,nGrains_loc)
      

      open(100,FILE='time.out')   
      nStep = 0
      iStep = 0
      do
         istep = istep + 1
         read(100,*,iostat=error) times(istep)
         if (error.ne.0) exit 
         nStep=nStep+1
      enddo
      close(100)
      write(*,*) '# time steps:',nStep
      
      open(202,file='cauchy_gra.out')
      open(203,file='nye_norm_gra.out')
      open(301,file='criticalGrainPairs.out')
      write(301,*) 'hard grain, soft grain, product of avg nye and avg cauchy'
      do iStep =  1, nStep  
         read(202,*) iTime,tTime
         read(203,*) iTime,tTime
         do iGrain=1,nGrains_loc
            read(202,*) grainID, cauchy(1:6)
            cauchyGrain(grainID) = cauchy(loadDir)
            read(203,*) grainID, nyeGrain(grainID)
            if (isGrainOnBoundary(grainID) .and. excludeBoundaryGrains) then
               nyeGrain(grainID) = 0.d0
               cauchyGrain(grainID) = 0.d0
            endif
         enddo
         
         ! calculate criticality if each grain pair
         measureCriticalityList = 0.d0
         do iList = 1,nTotNeighbors         
            HardGrainID = grainNeighborsList(iList,1)
            SoftGrainID = grainNeighborsList(iList,2)
            measureCriticalityList(iList) = cauchyGrain(HardGrainID)*nyeGrain(SoftGrainID)
         enddo
         
         CALL findNLargest(measureCriticalityList,nTotNeighbors,10,listN,posN)
                  
         write(301,*) iTime, tTime
         do iPos=1,10
            iList = posN(iPos)
            HardGrainID = grainNeighborsList(iList,1)
            SoftGrainID = grainNeighborsList(iList,2)
            write(301,*) HardGrainID, SoftGrainID, measureCriticalityList(iList), cauchyGrain(HardGrainID), nyeGrain(SoftGrainID)
         enddo
         ! also write nucleation measures for the Soft-Hard pair of interest
         if (exportNucleationMeasureHardGrain /= 0 .and. exportNucleationMeasureSoftGrain == 0) then
            write(301,*) 'nucleation measure for the Hard grain:', exportNucleationMeasureHardGrain
            do iList=1,nTotNeighbors
               if (grainNeighborsList(iList,1) == exportNucleationMeasureHardGrain) then
                  HardGrainID = grainNeighborsList(iList,1)
                  SoftGrainID = grainNeighborsList(iList,2)         
                  write(301,*) HardGrainID, SoftGrainID, &
                               measureCriticalityList(iList),&
                               cauchyGrain(HardGrainID), &
                               nyeGrain(SoftGrainID)
               endif   
            enddo
         elseif (exportNucleationMeasureHardGrain /= 0) then
            HardGrainID = grainNeighborsList(exportNucleationMeasureListPos,1)
            SoftGrainID = grainNeighborsList(exportNucleationMeasureListPos,2)         
            write(301,*) 'nucleation measure for the Hard-Soft Pair:', HardGrainID, SoftGrainID
            write(301,*) HardGrainID, SoftGrainID, &
                         measureCriticalityList(exportNucleationMeasureListPos),&
                         cauchyGrain(HardGrainID), &
                         nyeGrain(SoftGrainID)
         endif
      enddo
      close(301)
      
      write(*,*) 'done.'
      
      END PROGRAM

      

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


!***********************************************************************
!C     READS LINEAR ARRAY OF ELEMENT NUMBERS FOR EACH NODE - deniz
!***********************************************************************
      SUBROUTINE readNODE_ELEM_CONNEC(strConnectFile,IELEMNO,NINDX,NX,success)
      implicit none
      INCLUDE 'PARDIS.H'
      character(len=*), intent(in) :: strConnectFile      
      logical, intent(out) :: success
      integer, intent(out):: IELEMNO(MAXNODELEM*MAXNODE),NINDX(MAXNODE)
      integer, intent(in) :: NX
!***********************************************************************
!C     IELEMNO: STORES ELEMENT NUMBERS
!C     NINDX: STORES START AND END INDEX IN IELEMNO ARRAY FOR
!C            EACH NODE
!***********************************************************************      
      !Local
      logical :: fileExists
      integer :: error, iElem, nElem, iNode, nNodes, I, lastIndex

      success = .false.
      
      inquire(file=trim(strConnectFile),exist=fileExists)
      if(.NOT.fileExists) return
      
      open(716,file=trim(strConnectFile))
      read(716,*,iostat=error) nNodes
      if (error.ne.0 .or. nNodes.NE.NX) then
         close(716)
         return
      endif
      lastIndex = 1
      do iNode=1,nNodes
         read(716,*,iostat=error) nElem
         read(716,*,iostat=error) (IELEMNO(I),I=lastIndex,lastIndex+nElem-1)
         lastIndex = lastIndex + nElem
         NINDX(iNode) = lastIndex
      enddo
      if (error.ne.0 .or. nNodes.NE.NX) then
         close(716)
         return
      endif
      
      success = .true.
      close(716)
      
      END SUBROUTINE

!***********************************************************************
!C     CREATES LINEAR ARRAY OF ELEMENT NUMBERS FOR EACH NODE
!        modified by Kourosh for faster calculation
!***********************************************************************
      SUBROUTINE NODE_ELEM_CONNEC(IELEMNO,NINDX,IJK,NX,NELX,MNELX,MAXNODELEM)
      implicit none
      integer, intent(in) :: IJK(MNELX),NX,NELX,MNELX,MAXNODELEM
      integer, intent(out) :: IELEMNO(MAXNODELEM*NX),NINDX(NX)
!***********************************************************************
!C     IELEMNO: STORES ELEMENT NUMBERS
!C     NINDX: STORES START AND END INDEX IN IELEMNO ARRAY FOR
!C            EACH NODE
!***********************************************************************
      integer, allocatable :: IELEMNO_TMP(:,:),NINDX_TMP(:)
      integer :: nodeID, iNode, elemID, iElem, KOUNT

      ALLOCATE(IELEMNO_TMP(NX,MAXNODELEM),NINDX_TMP(NX))
      NINDX_TMP=0
      DO elemID=1,NELX
         DO iNode=1,4
            nodeID = IJK((elemID-1)*4+iNode)
            NINDX_TMP(nodeID)=NINDX_TMP(nodeID)+1
            IELEMNO_TMP(nodeID,NINDX_TMP(nodeID))=elemID
         ENDDO
      ENDDO
      
      KOUNT=1
      DO nodeID=1,NX
         DO iElem=1,NINDX_TMP(nodeID)
            IELEMNO(KOUNT)=IELEMNO_TMP(nodeID,iElem)
            KOUNT=KOUNT+1
         ENDDO
         NINDX(nodeID)=KOUNT
      ENDDO
      
      DEALLOCATE(IELEMNO_TMP,NINDX_TMP)

      
      END 
      
      
!***********************************************************************
!C     WRITES LINEAR ARRAY OF ELEMENT NUMBERS FOR EACH NODE - deniz
!***********************************************************************
      SUBROUTINE writeNODE_ELEM_CONNEC(strConnectFile,IELEMNO,NINDX, &
                                       NELX,NX,MAXNODELEM,MAXNODE)
                                       
      implicit none
      character(len=*), intent(in) :: strConnectFile      

      integer, intent(in) :: NX,NELX,MAXNODELEM,MAXNODE
      integer, intent(in) :: IELEMNO(MAXNODELEM*MAXNODE),NINDX(MAXNODE)
!***********************************************************************
!C     IELEMNO: STORES ELEMENT NUMBERS
!C     NINDX: STORES START AND END INDEX IN IELEMNO ARRAY FOR
!C            EACH NODE
!***********************************************************************      
      !Local
      logical :: fileExists
      integer :: error, iElem, nElem, iNode, I, lastIndex
      integer :: elemStart, elemEnd

      
      open(716,file=trim(strConnectFile))
      write(716,*,iostat=error) NX
      lastIndex = 1
      do iNode=1,NX
         elemStart=1
         if(iNode.ne.1) elemStart=NINDX(iNode-1)
         elemEnd = NINDX(iNode)-1
         nElem = elemEnd - elemStart + 1
         
         write(716,*,iostat=error) nElem
         write(716,*,iostat=error) (IELEMNO(iElem),iElem=elemStart,elemEnd)
         lastIndex = lastIndex + nElem
      enddo
      
      close(716)
      
      END SUBROUTINE
      
      SUBROUTINE createHistogram(R,nR,maxR,nBins,R_bins)
      implicit none
      real(8), intent(in) :: R(nR)
      integer, intent(in) :: nR
      real(8), intent(in) :: maxR
      integer, intent(in) :: nBins
      integer, intent(out):: R_bins(nBins)
      !locals
      integer :: iBin,iR
      real(8) :: R_step
      real(8) :: maxR_actual
      
      R_bins(:) = 0
      
      maxR_actual = MAXVAL(R)
      if (maxR_actual > maxR) write(*,*) 'maxR prob', maxR_actual,'>',maxR
      
      do iR = 1,nR
         iBin = (R(iR) / maxR_actual * nBins) + 1
         if (iBin == nBins + 1) iBin = nBins
         R_bins(iBin) = R_bins(iBin) + 1
      enddo
      
      END SUBROUTINE
      
      SUBROUTINE readGrainTexture(grain_rot,grain_euler,nGrains)
      implicit none
      
      real(8), intent(out):: grain_rot(3,3,nGrains), grain_euler(3,nGrains)
      integer, intent(in) :: nGrains
      
      integer :: iElem,iGrain,nGrains_in, iGrain_in
      real(8) :: euler(3)

      OPEN(UNIT=201,FILE='grainTexture.inp')
      read(201,*) nGrains_in
      
      if(nGrains_in /= nGrains) then
         write(*,*) 'nGrains in grainTexture is inconsistent with the rest of the input'
         stop
      endif

      DO iGrain=1,nGrains
         read(201,*) iGrain_in,grain_euler(1:3,iGrain)
         call eulerRot(grain_euler(:,iGrain),grain_rot(:,:,iGrain))
      ENDDO
      
      CLOSE(201)
      END SUBROUTINE
      


! output, rot,:
!   an active transformation matrix corresponding to the actual rotation due to euler angles {v}_i = TLGT_ij {v0}_j
!   as a (passive) coordinate transformation matrix, it transforms coordinates from crystal to global spatial frame. {v}_i = TLGT_ij {v_cry}_j
! input - euler angles
!     phi   - euler(1)
!     theta - euler(2)
!     omega - euler(3)
!---------------------------------------------------
      subroutine eulerRot(euler,rot)   ! restored to original
      implicit double precision (a-h,o-z)  

!      include 'aba_param.inc'

      dimension euler(3),tlg(3,3),rot(3,3)

      phi=euler(1)
      theta =euler(2)
      omega  =euler(3)

  
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

      rot = TRANSPOSE(tlg)   ! Kourosh

      return
      end   


      SUBROUTINE findNLargest(x_ref,Size,N,listN,posN)
      implicit none
      real(8), intent(in)                   :: x_ref(1:Size)
      integer, intent(in)                   :: Size
      integer, intent(in)                   :: N
      real(8), intent(out)                  :: listN(N)
      integer, intent(out)                  :: posN(N)
      
      INTEGER                               :: i
      INTEGER                               :: Location
      real(8) :: tempReal
      integer :: tempInteger
      
      logical :: marked(Size)
      integer :: iSearch
      integer :: extremaAt
      real(8) :: extremaThisPass
      
      marked = .false.
      extremaAt = 0
      
      do iSearch = 1, N
         extremaThisPass = tiny(extremaThisPass)
         do i = 1, Size
            if (x_ref(i) >= extremaThisPass .and. .not.marked(i)) then
               extremaThisPass = x_ref(i)
               extremaAt = i
            endif
         enddo
         marked(extremaAt) = .true.
         listN(iSearch) = extremaThisPass
         posN(iSearch) = extremaAt
      enddo
      
      END SUBROUTINE