

      program nodalRecoveryProg
      use meshtools
      use grains
      use nodalRecovery
      !use options
      
      implicit none
      include 'PARDIS.H'
      
      logical, parameter :: DEBUG_nodalValues = .true.
      integer :: p_deg
      integer :: recoveryMethod
      integer :: localNonlocal
      
      integer :: NX, NELX, NX_dupl, NODE
      integer :: nstepTime
      integer :: nstep_max, lastProgress
      integer :: nFailedPatches
      
      integer :: ENTRIESPERLINE, MAX_NX_DUPL      
      parameter(MAX_NX_DUPL=10*MAXNODE)

      real(8) :: nodalPositions_t      ! time of the imported positions
      integer :: nodalPositions_N      ! step number of the imported positions
      real(8), allocatable :: nodalPositions(:,:)    ! stores the node positions (read from disp.out)
      real(8), allocatable :: gaussValues(:,:)       ! stores the gauss values of the solution
      real(8), allocatable :: nodalValues(:,:)       ! stores the nodal values of the solution
      real(8), allocatable :: elemVol(:)
      real(8) :: times(100000)
      real(8) :: values_t              ! time of the imported solution
      integer :: values_N              ! step number of the imported solution
      logical :: found
      real(8) :: xNodes(3,4), elemJacMat(3,3), elemJac

      real(8) :: G0XYZ(MAXCRD)
      integer :: IJK(MNELX), IJKNEL(MAXELN)
      INTEGER:: IELEMNO(MAXNODELEM*MAXNODE),NINDX(MAXNODE)
      real(8) :: xc(3)
      integer :: info
      
   ! variables for domain geometry
      real(8) :: domainRange(3,2), domainLen(3)
      
      ! element groupings
      integer, allocatable:: elemFaceNeighbors(:,:)

      ! grains
      character(len=50) :: grainFileName
      character(len=150) :: lineStr
      logical :: fileExists, successReadNeighbors,successReadNodeConnect
      logical :: foundDispBlock, plotThisTime
      integer :: iGrain, iPatchNode
      integer, allocatable :: grainNodesAll(:)
      integer, allocatable :: grainNodeIndices(:,:)
      integer :: lastIndex, iNodeIndex, iGloNode

      real(8) :: xset(4),yset(4),zset(4), xpos, ypos, zpos
      real(8) :: distance, weight, totWeight
      integer :: elemNodes(4)
      integer :: val_ind
      character(len=30) :: meshName, fileName, solnFileName
      character :: commentsInData, commentsInDISP, timeawareData, hasGrainInfo
      integer :: iElem,iNode,iNodeBlock,iElemBlock
      integer :: iEuler,iDOF,nQuantDOF,nOutQuantDOF
      logical :: errorL
      integer :: len,i,j,k,kk,ITEMP,error
      integer :: step_N, this_N
      
!------------------------------------------------------------------------------      

      ! READ IN THE MESH, CONNECTIVITY
      CALL READGM(G0XYZ,IJK,NX,NELX,meshName)

      write(*,*) 'Mesh imported:'
      write(*,*) '# of elements: ', NELX
      write(*,*) '# of nodes: ', NX
            
      allocate(nodalPositions(3,NX))
      do iNode=1,NX
         nodalPositions(1:3,iNode)=G0XYZ((iNode-1)*3+1:iNode*3)
      enddo

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
      
      ! READ TIME STEPPING INFORMATION
      open(102,file='time.out')
      nstepTime=0
      times(:) = 0.d0
      do
         read(102,*,iostat=error) times(nstepTime+1)
         if(error.EQ.-1) EXIT
         nstepTime=nstepTime+1
      enddo
      close(102)
      if (nstepTime.EQ.0) then
         write(*,*) 'Cannot find file: time.out'
         STOP
      endif
      
      write(*,*) '# of time steps:', nstepTime

   !------------------------------------------------------------------------------

      write(*,*) 'Specify the file that contains the FE solution:'
      read(*,*) solnFileName
      write(*,*)

      write(*,*) 'Enter # of components of the field variable.'
      write(*,*) '(Temperature: 1, Cauchy stress: 6, Fp: 9, etc.)'
      read(*,*) nQuantDOF
      write(*,*)
      
      commentsInData=' '
      do while (commentsInData.NE.'Y'.AND.commentsInData.NE.'N')
         write(*,*) 'Are there comment/info lines in at the beginning of each time step in data files? (Y/N)'
         read(*,'(A1)') commentsInData
         if (commentsInData.EQ.'y') commentsInData='Y'
         if (commentsInData.EQ.'n') commentsInData='N'
      enddo
      write(*,*)

      commentsInDISP=' '
      do while (commentsInDISP.NE.'Y'.AND.commentsInDISP.NE.'N')
         write(*,*) 'Are there comment/info lines at the beginning of each time step in DISP.OUT? (Y/N)'
         read(*,'(A1)') commentsInDISP
         if (commentsInDISP.EQ.'y') commentsInDISP='Y'
         if (commentsInDISP.EQ.'n') commentsInDISP='N'
      enddo
      write(*,*)
      
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
      
      allocate(elemVol(NELX))
      do iElem=1,NELX
         CALL getElemXNodes(iElem,IJK,G0XYZ,NODE,xNodes)
         CALL calcJacTET4(xNodes, elemJacMat, elemJac, elemVol(iElem))   
      enddo
      
      CALL initialize_Patches('listPatches.dat',nPatches, &
                                  NX,NELX,IJK,G0XYZ)
      if (patchesInitialized) then
         write(*,*) 'patches imported:',nPatches
      else
         write(*,*) 'error importing listPatches.dat'
         stop
      endif
            
      allocate(gaussValues(nQuantDOF,NELX))
      allocate(nodalValues(nQuantDOF,nPatches))
      gaussValues(:,:) = 0.D0
      nodalValues(:,:) = 0.D0

      
      ! open input/output fules
      open(103,file=trim(solnFileName)//'.out')
      open(603,file=trim(solnFileName)//'_N.out')
      
      write(603,*) nPatches
      
      write(*,*) 'p_deg ?'
      read(*,*) p_deg
      
      write(*,*) 'recoveryMethod?'
      read(*,*) recoveryMethod
      
      write(*,*) 'localNonlocal?'
      read(*,*) localNonlocal
      
      
      ! now plot the field variables  !
      ! ------------------------------!
      write(*,*) 'Processing...'
      lastProgress=0
      do i=1,MIN(60,nstepTime)
         write(*,'(A)',advance='no') '_'
      enddo
      write(*,*) ''
         
      nFailedPatches = 0
      
      val_ind = 0

      do while (val_ind < nstepTime)
      
         ! -------- read the solution at -------!
         ! --------- gauss points --------!

         ! ---------- IMPORT THE SOLUTION ----------!
         
         val_ind = val_ind + 1

         error = 0
      
         if(commentsInData.EQ.'Y') then
            read(103,*,iostat=error) values_N, values_t
            if(error.NE.0) exit
         else
            values_N = val_ind
            values_t = val_ind
         endif
         do iElem = 1, NELX      
            read(103,*,iostat=error) (gaussValues(iDOF,iElem),iDOF=1,nQuantDOF)
            if(error.NE.0) exit
         enddo

         if(error.NE.0) then
            write(*,*) 'Solution file cannot be opened, or is corrupt.'
            stop
         end if
         
      
         ! recover nodal values  !
         ! ----------------------!   
         
         ! extract gauss values of Fp
         if(recoveryMethod.EQ.1)then
            ! nodal averaging
            
            nodalValues = 0.d0   ! good to initialize before nodal recovery
            do iGrain=1,nGrains
               ! recover nodal values for each individual grain

               CALL nodalAveraging(nodalValues,gaussValues, &
                                   nQuantDOF,DEBUG_nodalValues)
            enddo
         elseif(recoveryMethod.EQ.2)then
            ! SPR method
                           
            !p_deg, nNeighDeg, methodSPR imported from options.inp
          
            ! method used in case of insufficient elems per patch. 
            ! 0, do not apply patch, rely on contributions from surr. patches
            ! 1, extend the patch to more layers of neighbors
            ! 2, also increase p_deg while extending the patch
            ! recover nodal values for each individual grain
            CALL recoverySPR(nodalValues,gaussValues, &
        nQuantDOF,p_deg,localNonlocal,nFailedPatches, &
                              info,DEBUG_nodalValues)

         endif
! -------------- SUCCESSFULLY RECOVERED NODAL VALUES ------- !
         write(603,*) values_N, values_t
         do iPatchNode=1,nPatches
            write(603,*) nodalValues(:,iPatchNode)
         enddo


         ! --------- DONE ------- !
         !progress bar fancyness
         if((60*val_ind)/nstepTime.GT.lastProgress) then
            write(*,'(A)',advance='no') '*'
            call flush(6)
            lastProgress=(60*val_ind)/nstepTime
         endif
         
      enddo
      
      if (nFailedPatches /= 0) then
         write(*,*) '# of failed patches:',nFailedPatches
      endif

      ! close files
      close(101)
      close(103)
      close(603)
      write(*,*)
      
      CALL grains_Destruct()
      ! element groupings
      deallocate(elemFaceNeighbors)
      
      deallocate(elemVol)
      
      deallocate(nodalPositions)
      deallocate(nodalValues)
      deallocate(gaussValues)
      if (allocated(grainNodesAll)) deallocate(grainNodesAll)
      if (allocated(grainNodeIndices)) deallocate(grainNodeIndices)
     
      END PROGRAM
      
      ! calculates the TET interpolation of a field Phi(s) at the natural coordinates s,
      ! given the nodal values
      SUBROUTINE interpolateTETs(sPoint,interpValue,nodalValuesElem, & 
                                 nQuantDOF)
      implicit none
      real(8), intent(in) :: sPoint(3)
      real(8), intent(out):: interpValue(nQuantDOF)
      real(8), intent(in) :: nodalValuesElem(nQuantDOF,4)
      integer, intent(in) :: nQuantDOF
      
      real(8) :: N(4)
      
      N(1) = sPoint(1)
      N(2) = sPoint(2)
      N(3) = sPoint(3)
      N(4) = 1-sPoint(1)-sPoint(2)-sPoint(3)
      
      interpValue = MATMUL(nodalValuesElem,N)
      
      END SUBROUTINE
      
      ! calculates the TET interpolation of a field Phi(s) at the spatial coordinates x,
      ! given the nodal values Phi_n and element nodal positions x_n
      SUBROUTINE interpolateTETx(sPoint,interpValue,nodalValuesElem, & 
                                 nQuantDOF)
      implicit none
      real(8), intent(in) :: sPoint(3)
      real(8), intent(out):: interpValue(nQuantDOF)
      real(8), intent(in) :: nodalValuesElem(nQuantDOF,4)
      integer, intent(in) :: nQuantDOF
      
      ! NOT YET CODED
      interpValue = 0.d0
      
      END SUBROUTINE

      SUBROUTINE solveAXb(A,N,b,err)
      integer, intent(in):: N
      real(8), intent(in):: A(N,N)
      real(8), intent(inout):: b(N)
      integer, intent(out):: err
         CALL DPOSV('U',N,1,A,N,b,N,err)
      END SUBROUTINE
      
      function factorial(n) result(n_fact)
         integer, intent(in) :: n
         integer:: n_fact
         n_fact = 1
         do i=1,n
            n_fact=n_fact*i
         enddo
      end function factorial

      SUBROUTINE printIntMatrix(A,N,M,matName,printTo)
      implicit none
      integer, intent(in):: N,M
      integer, intent(in):: A(N,M)
      character*(*),intent(in):: matName
      integer, intent(in):: printTo
      integer:: i,j
      
      if(LEN_TRIM(matName).GT.0)then
      
         if(printTo.GT.0) then
            write(printTo,*) matName,':'
         else
            write(*,*) matName,':'
         endif
      endif
      do i=1,N
         if(printTo.GT.0) then
            write(printTo,*) (A(i,j), j=1,M)
         else
            write(*,*) (A(i,j), j=1,M)
         endif
      enddo
      END SUBROUTINE
      
      SUBROUTINE printVector(V,N,nEachRow,vecName,printTo)
      implicit none
      integer, intent(in):: N,nEachRow
      real(8), intent(in):: V(N)
      character*(*),intent(in):: vecName
      integer, intent(in):: printTo
      integer:: i,j
      if(LEN_TRIM(vecName).GT.0)then
         if(printTo.GT.0) then
            write(printTo,*) vecName,':'
         else
            write(*,*) vecName,':'
         endif
      endif
      do i=1,N/nEachRow
         if(printTo.GT.0) then
            write(printTo,*) (V((i-1)*nEachRow+j), j=1,nEachRow)
         else
            write(*,*) (V((i-1)*nEachRow+j), j=1,nEachRow)
         endif
      enddo
      ! print the remaining coefficients, if any
      if(N-(i-1)*nEachRow.NE.0)then
         if(printTo.GT.0) then
            write(printTo,*) (V((i-1)*nEachRow+j), j=1,N-(i-1)*nEachRow)
         else
            write(*,*) (V((i-1)*nEachRow+j), j=1,N-(i-1)*nEachRow)
         endif
      endif
      END SUBROUTINE
      
      SUBROUTINE printIntVector(V,N,nEachRow,vecName,printTo)
      implicit none
      integer, intent(in):: N,nEachRow
      integer, intent(in):: V(N)
      character*(*),intent(in):: vecName
      integer, intent(in):: printTo
      integer:: i,j
      if(LEN_TRIM(vecName).GT.0)then
         if(printTo.GT.0) then
            write(printTo,*) vecName,':'
         else
            write(*,*) vecName,':'
         endif
      endif
      do i=1,N/nEachRow
         if(printTo.GT.0) then
            write(printTo,*) (V((i-1)*nEachRow+j), j=1,nEachRow)
         else
            write(*,*) (V((i-1)*nEachRow+j), j=1,nEachRow)
         endif
      enddo
      ! print the remaining coefficients, if any
      if(N-(i-1)*nEachRow.NE.0)then
         if(printTo.GT.0) then
            write(printTo,*) (V((i-1)*nEachRow+j), j=1,N-(i-1)*nEachRow)
         else
            write(*,*) (V((i-1)*nEachRow+j), j=1,N-(i-1)*nEachRow)
         endif
      endif
      END SUBROUTINE

      SUBROUTINE printMatrix(A,N,M,matName,printTo)
      implicit none
      integer, intent(in):: N,M
      real(8), intent(in):: A(N,M)
      character*(*),intent(in):: matName
      integer, intent(in):: printTo
      integer:: i,j
      if(LEN_TRIM(matName).GT.0)then
         if(printTo.GT.0) then
            write(printTo,*) matName,':'
         else
            write(*,*) matName,':'
         endif
      endif
      do i=1,N
         if(printTo.GT.0) then
            write(printTo,*) (A(i,j), j=1,M)
         else
            write(*,*) (A(i,j), j=1,M)
         endif
      enddo
      END SUBROUTINE
      

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
      logical :: fileExists, successReadNodeConnect
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

      
      END SUBROUTINE