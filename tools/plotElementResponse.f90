! development log
!----------------
!
! June 2015 - grain IDs, phases and Euler angles are plotted for inspecting consistency between D3D file and the input files (needs debug)

      program PostProcessing
      use meshtools
      use grains
      implicit none
      include 'PARDIS.H'
      
      integer :: NX, NELX, NODE, elemID
      real(8) :: pointX(3)
      real(8) :: xElemNodes(3,4)
      logical :: inside
      integer :: nstepTime, nstepInp, nstepTP,nstep_max, lastProgress, beginStepTP,nstepsSpecific
      integer, allocatable :: plotTimeSteps(:)
      real(8) :: plotLastSeconds
      integer :: ENTRIESPERLINE      
      parameter(ENTRIESPERLINE = 500)
      
      integer :: valuesAtNodes
      integer, parameter :: valuesUnique = 1    ! enum for valuesAtNodes
      integer, parameter :: valuesPerGrain = 2

      real(8), allocatable :: fieldValue(:)       ! stores the gauss values of the solution
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
      integer :: iEuler,iDOF,nQuantDOF,nOutQuantDOF
      integer, allocatable :: outQuantDOF(:)
      logical, allocatable :: export_QuantDOF(:)
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
      write(*,*)

   !------------------------------------------------------------------------------

      write(*,*) 'Specify the file that contains the field data:'
      read(*,*) solnFileName
      write(*,*)

      write(*,*) 'Enter # of components of the field variable.'
      write(*,*) '(Temperature: 1, Cauchy stress: 6, Fp: 9, etc.)'
      read(*,*) nQuantDOF
      write(*,*)
      
      allocate(outQuantDOF(nQuantDOF))
      allocate(export_QuantDOF(nQuantDOF))
      
      write(*,*) 'specify the components of the field variable you''d like to export:'
      write(*,*) 'example: F F T F T F --> 3rd and 5th components. T: export, F: do not export'
      read(*,*) export_QuantDOF(1:nQuantDOF)
      write(*,*)
      nOutQuantDOF=0
      do iDOF=1,nQuantDOF
         if(export_QuantDOF(iDOF))then
            nOutQuantDOF=nOutQuantDOF+1
            outQuantDOF(nOutQuantDOF) = iDOF
         endif
      enddo
      
      write(*,*) 'Enter element ID to extract the evolution of field values:'
      write(*,*) 'or enter 0 to search for the element with point coordinates'
      read(*,*) elemID
      write(*,*)
      
      if (elemID == 0) then 
         write(*,*) 'Enter point coordinates (X,Y,Z) to pick an element:'
         read(*,*) pointX(1:3)
         do iElem=1,NELX
            ! check if point is inside element iElem
            CALL getElemXNodes(iElem,IJK,G0XYZ,NODE,xElemNodes)
            CALL checkInside(xElemNodes,pointX,inside)
            if (inside) then
               elemID = iElem
               write(*,*) 'found element: ',iElem
               write(*,*) 
               exit
            endif
         enddo
      endif
      
      commentsInData=' '
      do while (commentsInData.NE.'Y'.AND.commentsInData.NE.'N')
         write(*,*) 'Are there comment/info lines in at the beginning of each time step in data files? (Y/N)'
         read(*,'(A1)') commentsInData
         if (commentsInData.EQ.'y') commentsInData='Y'
         if (commentsInData.EQ.'n') commentsInData='N'
      enddo
      write(*,*)

      timeawareData='Y'
      do while (timeawareData.NE.'Y'.AND.timeawareData.NE.'N')
         write(*,*) 'include time information in techPlot output?'
         write(*,*) 'Requirements:'
         write(*,*) '  1) solution field has time information (SOLUTIONTIME field in ZONE header)'
         write(*,*) '  2) disp.out contains the displacements at the corresponding times'
         read(*,'(A1)') timeawareData
         if (timeawareData.EQ.'y') timeawareData='Y'
         if (timeawareData.EQ.'n') timeawareData='N'
      enddo
      write(*,*)
      
      beginStepTP = 1
      nstepTP = nstepTime ! # of time steps to export to TP file.

      
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
      
      !read grain sizes: grainSize(:)
      CALL readGrainSizes('grainSizes.inp')
      if(grainSizesImported) write(*,*) 'grain sizes imported'
            
      !read grain phases: grainPhase(:)
      CALL readGrainPhases('grainPhases.inp')
      if(grainPhasesImported) write(*,*) 'grain phases imported'
      
      !read texture: grainTexture(:)
      CALL readTexture('grainTexture.inp')
      if(textureImported) write(*,*) 'texture imported'
      
      allocate(elemVol(NELX))
      do iElem=1,NELX
         CALL getElemXNodes(iElem,IJK,G0XYZ,NODE,xNodes)
         CALL calcJacTET4(xNodes, elemJacMat, elemJac, elemVol(iElem))   
      enddo
      
      ! calculate the nodal values, if not provided
      ! removed in this version.
      ! if nodal values are given, field is plotted through the nodal points
      ! if gauss values are given, field is plotted through the gauss points (TECPLOT: CELL-centered variables, Block data packing)
   
   
      ! allocate buffers, open files for in/out
      
      open(103,file=trim(solnFileName)//'.out')
      
      allocate(fieldValue(nQuantDOF))
      write(strOutputFileName,'(A,A,I0,A)') trim(solnFileName),'_elem',elemID,'.out'
      open(603, file=trim(strOutputFileName))

      
      ! OUTPUT AS TEXT FILE
      ! ------------------------------!
      
      ! now plot the field variables  !
      ! ------------------------------!
      write(*,*) 'Processing...'
      lastProgress=0
      do i=1,MIN(60,nstepTP)
         write(*,'(A)',advance='no') '_'
      enddo
      write(*,*) ''
         
      nstepInp=0
      val_ind = 0
      disp_ind = 0
      connExported = .false.
      
      do while (disp_ind < MIN(nstepTime,nstepTP))
      
         

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
            if (iElem==elemID) then
               read(103,*,iostat=error) fieldValue(1:nQuantDOF)
            else
               call READSTR(103,lineStr,error)
            endif
            if(error.NE.0) exit
         enddo

         if(error.NE.0) then
            write(*,*) 'Solution file cannot be opened, or is corrupt.'
            stop
         end if


         ! now plot the field variables  !
         ! ------------------------------!

      
         ! EXPORT GAUSS PT VALUES / CELL-CENTERED VALUES
         ! -------------------------------------------------- !
         do j=1,nOutQuantDOF
            write(603,*) values_t, fieldValue(outQuantDOF(j))
         enddo  

         !progress bar fancyness
         if((60*val_ind)/nstepTP.GT.lastProgress) then
            write(*,'(A)',advance='no') '*'
            call flush(6)
            lastProgress=(60*val_ind)/nstepTP
         endif
         
      enddo
      write(*,*)
      
      write(*,*) 'Exported output file: ', strOutputFileName

      ! close files
      close(603)
      close(103)
      
      
      CALL grains_Destruct()
      ! element groupings
      deallocate(elemFaceNeighbors)
      
      deallocate(elemVol)
      
      deallocate(outQuantDOF,export_QuantDOF)
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