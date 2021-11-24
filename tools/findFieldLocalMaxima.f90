! this program find the maximum

      program PostProcessing
      use meshtools
      use grains
      implicit none
      include 'PARDIS.H'
      
      integer :: NX, NELX, NODE, elemID
      real(8) :: pointX(3)
      real(8) :: xElemNodes(3,4)
      logical :: inside
      integer :: nstepTime, nstepInp,nstep_max, lastProgress
      integer, allocatable :: plotTimeSteps(:)
      real(8) :: plotLastSeconds
      integer :: ENTRIESPERLINE      
      parameter(ENTRIESPERLINE = 500)
      
      integer :: valuesAtNodes
      integer, parameter :: valuesUnique = 1    ! enum for valuesAtNodes
      integer, parameter :: valuesPerGrain = 2

      real(8), allocatable :: fieldValue(:)       ! stores the gauss values of the solution
      real(8) :: maxFieldValue
      integer :: maxFieldValueAtElemID
      real(8) :: maxFieldValueAtTime
      real(8) :: times(100000), plotBeginTime
      real(8) :: values_t              ! time of the imported solution
      integer :: values_N              ! step number of the imported solution
      integer :: connZone
      logical :: connExported, found
      real(8), allocatable :: AveStress(:,:)
      real(8) :: xNodes(3,4), elemJacMat(3,3), elemJac
      real(8) :: G0XYZ(MAXCRD)
      integer :: IJK(MNELX), IJKNEL(MAXELN)
      INTEGER:: IELEMNO(MAXNODELEM*MAXNODE),NINDX(MAXNODE)
      real(8) :: xc(3)
      
   ! variables for domain geometry
      real(8) :: domainRange(3,2), domainLen(3)
      
      ! RVE-averaged quantities
      real(8), allocatable :: elemVol(:)
      real(8), allocatable :: elemCenters(:,:)
      real(8) :: totVol
      
      ! element groupings
      integer, allocatable:: elemFaceNeighbors(:,:)

      ! grains
      character(len=50) :: grainFileName
      character(len=150) :: lineStr
      logical :: fileExists, successReadNeighbors,successReadNodeConnect
      integer :: iGrain
      integer :: lastIndex, iNodeIndex, iGloNode
      integer :: iGrainNode, nNodes_file

      real(8) :: xset(4),yset(4),zset(4), xpos, ypos, zpos
      real(8) :: distance, weight, totWeight
      integer :: elemNodes(4)
      integer :: val_ind, outputTimeStepSkip
      character(len=30) :: meshName, fileName, solnFileName, strOutputFileName
      character :: plotConfig, commentsInData, timeawareData, hasGrainInfo
      integer :: iElem,iNode,iNodeBlock,iElemBlock
      integer :: iEuler,iDOF,nQuantDOF,searchQuantDOF
      logical :: errorL
      integer :: len,i,j,k,kk,ITEMP,error
      integer :: step_N, this_N, disp_ind
      
!------------------------------------------------------------------------------      

      ! READ IN THE MESH, CONNECTIVITY
      CALL READGM(G0XYZ,IJK,NX,NELX,meshName)      

      write(*,*) 'Mesh imported:'
      write(*,*) '# of elements: ', NELX
      write(*,*) '# of nodes: ', NX
      
      NODE = 4
      
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
      write(*,*)

   !------------------------------------------------------------------------------

      write(*,*) 'Specify the file that contains the field data:'
      read(*,*) solnFileName
      write(*,*)

      write(*,*) 'Enter # of components of the field variable.'
      write(*,*) '(Temperature: 1, Cauchy stress: 6, Fp: 9, etc.)'
      read(*,*) nQuantDOF
      write(*,*)
      
      write(*,*) 'Enter the component index to search for the extremum'
      read(*,*) searchQuantDOF
      write(*,*)
      
      commentsInData=' '
      do while (commentsInData.NE.'Y'.AND.commentsInData.NE.'N')
         write(*,*) 'Are there comment/info lines in at the beginning of each time step in data files? (Y/N)'
         read(*,'(A1)') commentsInData
         if (commentsInData.EQ.'y') commentsInData='Y'
         if (commentsInData.EQ.'n') commentsInData='N'
      enddo
      write(*,*)

      allocate(elemVol(NELX))
      allocate(elemCenters(3,NELX))
      do iElem=1,NELX
         CALL getElemXNodes(iElem,IJK,G0XYZ,NODE,xNodes)
         CALL getElemXCenter(iElem,IJK,G0XYZ,NODE,elemCenters(:,iElem))
      enddo
      

      
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
      
      open(103,file=trim(solnFileName)//'.out')
      
      allocate(fieldValue(nQuantDOF))
      
      ! OUTPUT AS TEXT FILE
      ! ------------------------------!
      
      ! now plot the field variables  !
      ! ------------------------------!
      write(*,*) 'Processing...'
      lastProgress=0
      do i=1,MIN(60,nstepTime)
         write(*,'(A)',advance='no') '_'
      enddo
      write(*,*) ''
         
      nstepInp=0
      val_ind = 0
      disp_ind = 0
      connExported = .false.
      
      maxFieldValue = 0.d0
      maxFieldValueAtElemID = 0
      maxFieldValueAtTime = 0.d0
      do while (disp_ind < nstepTime)
      
         

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
            read(103,*,iostat=error) fieldValue(1:nQuantDOF)
            if(error.NE.0) exit
            if (dabs(fieldValue(searchQuantDOF)) > maxFieldValue .and. &
                elemCenters(1,iElem) > domainRange(1,1)+domainLen(1)*0.1d0 .and. &
                elemCenters(1,iElem) < domainRange(1,2)-domainLen(1)*0.1d0 .and. &
                elemCenters(2,iElem) > domainRange(2,1)+domainLen(2)*0.1d0 .and. &
                elemCenters(2,iElem) < domainRange(2,2)-domainLen(2)*0.1d0 .and. &
                elemCenters(3,iElem) > domainRange(3,1)+domainLen(3)*0.1d0 .and. &
                elemCenters(3,iElem) < domainRange(3,2)-domainLen(3)*0.1d0 ) then
               maxFieldValue = dabs(fieldValue(searchQuantDOF))
               maxFieldValueAtElemID = iElem
               maxFieldValueAtTime = values_t
            endif
         enddo

         if(error.NE.0) then
            write(*,*) 'Solution file cannot be opened, or is corrupt.'
            stop
         end if
         
         !progress bar fancyness
         if((60*val_ind)/nstepTime.GT.lastProgress) then
            write(*,'(A)',advance='no') '*'
            call flush(6)
            lastProgress=(60*val_ind)/nstepTime
         endif
         
      enddo
      write(*,*)
      
      write(*,*) 'max local field value:',maxFieldValue
      write(*,*) 'Element at which maximum field value occurs: ',maxFieldValueAtElemID
      write(*,*) 'centered at: ',elemCenters(:,maxFieldValueAtElemID)
      write(*,*) 'Time at which maximum field value occurs: ',maxFieldValueAtTime

      ! close files
      close(103)
      
      CALL grains_Destruct()
      ! element groupings
      deallocate(elemFaceNeighbors)
      
      deallocate(elemVol)
      deallocate(elemCenters)
      
      if (allocated(fieldValue)) deallocate(fieldValue)
      
      end
        
      
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

      
      END SUBROUTINE