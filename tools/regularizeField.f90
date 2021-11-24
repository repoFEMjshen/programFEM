      ! this program applies regularization/averaging filters on FEM fields: gaussian averaging, grain-averaging
      ! also calculates and exports time-evolving rve-averages of the fields
      program regularizeField
      
      use meshtools
      use grains
      use regularization
      
      implicit none
      include 'PARDIS.H'
      
      integer :: regularizationMethod
      real(8) :: characteristicLength
      
      integer :: NX, NELX, NODE
      integer :: nstepTime
      integer :: nstep_max, lastProgress
      
      real(8), allocatable :: fieldLocal(:,:)       ! stores the gauss values of the solution
      real(8), allocatable :: fieldGrains(:,:)      ! stores the grain-averaged values of the solution, per grain
      real(8), allocatable :: fieldRegularized(:,:)       ! stores the gauss values of the solution (regularized)
      real(8), allocatable :: rveAverage(:)        ! stores the rve-averaged value of the field
      real(8) :: times(100000)
      real(8) :: values_t              ! time of the imported solution
      integer :: values_N              ! step number of the imported solution
      logical :: found

      real(8) :: G0XYZ(MAXCRD)
      integer :: IJK(MNELX), IJKNEL(MAXELN)
      INTEGER:: IELEMNO(MAXNODELEM*MAXNODE),NINDX(MAXNODE)
      real(8) :: xc(3)
      integer :: info
      
   ! variables for domain geometry
      real(8) :: domainRange(3,2), domainLen(3)
      
      ! element groupings
      integer, allocatable:: elemFaceNeighbors(:,:)

      ! export time-series data for specific grains/components
      character :: exportSpecificComponents_Grains
      logical, allocatable :: export_QuantDOF(:)
      integer, allocatable :: outQuantDOF(:)
      integer :: nOutGrains
      integer, allocatable :: export_GrainIDs(:)

      ! grains
      character(len=50) :: grainFileName
      character(len=150) :: lineStr
      logical :: fileExists, success,successReadNeighbors,successReadNodeConnect
      logical :: foundDispBlock, plotThisTime
      integer :: iGrain, iPatchNode
      integer :: lastIndex, iNodeIndex, iGloNode

      real(8) :: xset(4),yset(4),zset(4), xpos, ypos, zpos
      real(8) :: distance, weight, totWeight
      integer :: elemNodes(4)
      integer :: val_ind
      character(len=30) :: meshName, fileName, solnFileName
      character :: commentsInData, exportRegularizedField, timeawareData, hasGrainInfo
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
      
      !read grain sizes: grainSize(:) - Deniz - Sept 22
      CALL readGrainSizes('grainSizes.inp')
      if(grainSizesImported) write(*,*) 'grain sizes imported'
      if(.NOT.grainSizesImported) then
         write(*,*) 'WARNING ! no grain size information (grainSizes.inp) found. assuming default grain size of:',grainSize(1)
      endif
      
      write(*,*) 'regularizationMethod? 1: grain-AVG, 2: gaussian averaging, 3: rve-averaging'
      read(*,*) regularizationMethod
      
      if (regularizationMethod==reg_RVEAVG .or. &
          regularizationMethod==reg_GRAINAVG) then
            
         exportRegularizedField=' '
         do while (exportRegularizedField.NE.'Y'.AND.exportRegularizedField.NE.'N')
            write(*,*) 'export the local, regularized field (large file), besides the domain averages (small file)? (Y/N)'
            read(*,'(A1)') exportRegularizedField
            if (exportRegularizedField.EQ.'y') exportRegularizedField='Y'
            if (exportRegularizedField.EQ.'n') exportRegularizedField='N'
         enddo
         write(*,*)
      else
         exportRegularizedField = 'Y'
      endif
      
      exportSpecificComponents_Grains = ' '
      if (regularizationMethod==reg_GRAINAVG) then
      
         do while(exportSpecificComponents_Grains /= 'Y' .and. &
                  exportSpecificComponents_Grains /= 'N' .and. &
                  exportSpecificComponents_Grains /= '0')
                     
            write(*,*) 'export specific components of the field at specific grains as time series?'
            write(*,*) 'Y: choose grains, components. 0: export all. N: do not export'
            read(*,*) exportSpecificComponents_Grains
            write(*,*) ''
            if (exportSpecificComponents_Grains=='y') exportSpecificComponents_Grains = 'Y'
            if (exportSpecificComponents_Grains=='n') exportSpecificComponents_Grains = 'N'
         enddo
         
         if (exportSpecificComponents_Grains == 'Y') then
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
            
            write(*,*) 'Enter the total # grains you would like to export:'
            read(*,*) nOutGrains
            write(*,*)
            allocate(export_GrainIDs(nOutGrains))
            write(*,*) 'Enter grain IDs that you would like to export a time series for:'
            read(*,*) export_GrainIDs(1:nOutGrains)
            write(*,*)
                        
         elseif (exportSpecificComponents_Grains == '0') then
            allocate(outQuantDOF(nQuantDOF))
            nOutQuantDOF = nQuantDOF
            do iDOF=1,nQuantDOF
               outQuantDOF(iDOF) = iDOF
            enddo
            nOutGrains = nGrains
            allocate(export_GrainIDs(nOutGrains))
            do iGrain=1,nGrains
               export_GrainIDs(iGrain) = iGrain
            enddo
         endif
      endif
     
      if (regularizationMethod==reg_GAUSSIAN) then
         write(*,*) 'characteristicLength?'
         read(*,*) characteristicLength
         
         CALL generateElementNonLocalList(characteristicLength,grainSize(:), &
                                          NELX,NX,IJK,G0XYZ,&
                                          grainIDelem,grainElementIndx,grainElements,nGrains,&
                                          1,nGrains,0,success)      
                                          
         if (success) then
            write(*,*) 'generated non-local element list'
         else
            write(*,*) 'error generating non-local element list'
            stop
         endif
      endif

      ! calculate element centers, volumes
            
      allocate(fieldLocal(nQuantDOF,NELX))
      allocate(fieldGrains(nQuantDOF,nGrains))
      allocate(fieldRegularized(nQuantDOF,NELX))
      allocate(rveAverage(nQuantDOF))
      fieldLocal(:,:) = 0.D0
      fieldGrains(:,:) = 0.d0
      fieldRegularized(:,:) = 0.D0
      rveAverage(:) = 0.d0
      
      ! open input/output fules
      open(103,file=trim(solnFileName)//'.out')
      if (exportRegularizedField=='Y') open(603,file=trim(solnFileName)//'_reg.out')
      if (regularizationMethod==reg_GRAINAVG) open(604,file=trim(solnFileName)//'_gra.out')
      if (regularizationMethod==reg_GRAINAVG .and. &
          exportSpecificComponents_Grains /= 'N') then
         open(606,file=trim(solnFileName)//'_gra_t.out')
         write(606,*) 'selected grains:', export_GrainIDs(:)
         write(606,*) 'selected components:', outQuantDOF(1:nOutQuantDOF)
         write(606,*) 'time selected-grain-1-components.. selected-grain-2-components..'
      endif
      if (regularizationMethod==reg_RVEAVG) open(605,file=trim(solnFileName)//'_rve.out')
      
      ! now plot the field variables  !
      ! ------------------------------!
      write(*,*) 'Processing...'
      lastProgress=0
      do i=1,MIN(60,nstepTime)
         write(*,'(A)',advance='no') '_'
      enddo
      write(*,*) ''
      
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
            if (exportRegularizedField=='Y') write(603,*) values_N, values_t
         else
            values_N = val_ind
            values_t = val_ind
         endif
         do iElem = 1, NELX      
            read(103,*,iostat=error) (fieldLocal(iDOF,iElem),iDOF=1,nQuantDOF)
            if(error.NE.0) then
               write(*,*) 'Solution file cannot be opened, or is corrupt.'
               stop
            end if
         enddo
      
         ! regularize field/calculated averages  !
         ! ------------------!   
         if(regularizationMethod==reg_GRAINAVG) then ! grain-averaging
         
            CALL regularizeFieldGrainAVG(fieldLocal, nQuantDOF, fieldRegularized, &
                                         fieldGrains, &
                                         NELX,NX,IJK,G0XYZ,&
                                         grainIDelem,grainElementIndx,grainElements,nGrains,&
                                         1,nGrains, &
                                         success)
         
         elseif(regularizationMethod==reg_GAUSSIAN) then ! gaussian regularization
         
            CALL regularizeFieldGaussian(fieldLocal, nQuantDOF, fieldRegularized, success)
         elseif(regularizationMethod==reg_RVEAVG)  then  ! rve-average

            CALL calculateRveAverage(fieldLocal, nQuantDOF, fieldRegularized, &
                                     rveAverage, &
                                     NELX,NX,IJK,G0XYZ,&
                                     success)         
         endif
         
         if (.not.success) then
            write(*,*) 'regularization failed'
            exit
         endif
         
         if (exportRegularizedField=='Y') then
            do iElem=1,NELX
               write(603,*) fieldRegularized(:,iElem)
            enddo
         endif
         
         ! print grain-averaged values (per grain)
         if (regularizationMethod==reg_GRAINAVG) then
            write(604,*) values_N, values_t
            do iGrain=1,nGrains
               write(604,*) iGrain, fieldGrains(:,iGrain)
            enddo
            if (exportSpecificComponents_Grains /= 'N') then
               write(606,*) values_t, fieldGrains(outQuantDOF(1:nOutQuantDOF), &
                                                  export_GrainIDs(1:nOutGrains))
            endif
         elseif(regularizationMethod==reg_RVEAVG) then
            ! print rve-averaged field
            write(605,*) values_N, values_t, rveAverage(:)
         endif
         
         ! --------- DONE ------- !
         !progress bar fancyness
         if((60*val_ind)/nstepTime.GT.lastProgress) then
            write(*,'(A)',advance='no') '*'
            call flush(6)
            lastProgress=(60*val_ind)/nstepTime
         endif
         
      enddo

      ! close files
      close(103)
      if (exportRegularizedField=='Y') close(603)
      if (regularizationMethod==reg_GRAINAVG) close(604)
      if (regularizationMethod==reg_GRAINAVG .and. &
          exportSpecificComponents_Grains /= 'N') close(606)
      if (regularizationMethod==reg_RVEAVG) close(605)

      write(*,*)
      
      CALL reg_destruct()
      CALL grains_Destruct()
      
      ! element groupings
      deallocate(elemFaceNeighbors)
      
      deallocate(rveAverage)
      
      deallocate(fieldLocal)
     
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