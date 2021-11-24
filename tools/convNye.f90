      PROGRAM convNye
      
      use meshtools
      !use optionsCI
      
      use grains
      use nodalRecovery

      IMPLICIT NONE
      
      include 'mpif.h'    ! KOUROSH_MPI
      INCLUDE 'PARDIS.H'
      
      ! instead of optionsCI
      logical :: commentsInDataFile = .true.
      integer :: nStepsToAnalyse

      real(8), parameter :: PI=3.14159265359D0

      REAL(8):: G0XYZ(MAXCRD)
      
      INTEGER::IERROR,IPROC_MPI,NPROCS_MPI

      INTEGER :: IJK(MNELX)
      INTEGER :: IELEMNO(MAXNODELEM*MAXNODE),NINDX(MAXNODE)
      INTEGER :: NX,NELX,NEQ,NODE,MDIM,NGAUSS,NDF
      real(8) :: times(100000)
      integer :: staGrainIdx,endGrainIdx,nGrainIdx
      
      integer :: NFACES, NEDGES, NELEMPERBLOCK
   ! initialize arrays for nodal stress recovery
      real(8), allocatable :: gauss_values_Fp(:,:)
      real(8), allocatable :: nodal_values_Fp(:,:)
      real(8) :: elemNodalFp_T(3,3,4)
      real(8) :: elemNodalNyeArr(9,4)
      real(8) :: elemNodalNye(3,3,4)
      real(8) :: gradFp_T(3,3,3), gradCurlFp(3,3,3)
      real(8), allocatable:: elemNye(:,:)
      integer :: grainIdx
      
      real(8) :: Fp_t              ! time of the imported solution
      integer :: Fp_N              ! step number of the imported solution
      integer, parameter :: nDOF_Fp = 9
      
   ! crack nucleation, Stress, Nye
      real(8) :: Nye(3,3), Fp(3,3)
      real(8) :: dummy(3)
      real(8), allocatable:: nodalPositions(:,:)
      
      integer, allocatable:: elemFaceNeighbors(:,:)
      
   ! some dummy variables
      integer :: iElem,nElem,iNode,iGloNode,iLocNode,nodeIdx
      integer :: istep,nstep,nstepMax,nstepInp
      integer :: I,II,III,J,JJ,JJJ,K,iDim
      integer :: iStr,iDOF,iBlock,iGrain,jGrain,iGrGrBnd
      integer :: iGrainIdx,jGrainIdx
      logical :: fileExists, successReadNodeConnect, successReadNeighbors
      logical :: disp_found
      integer :: error, info, lastProgress
      
      logical :: errorL
      integer :: len,kk,ITEMP
      integer :: step_N, this_N, disp_ind, disp_Ind_last
      real(8) :: step_Time, this_Time
      
   ! other
      character(10) :: strNSTEP
      character(10) :: strPOSQ
      real(8) :: edgeLength
      real(8) :: xNodes(3,4), xNode(3)
      real(8), allocatable :: xElemCenters(:,:) !(NELX,3)      
      real(8) :: elemJacMat(3,3)
      real(8) :: elemJacInv(3,3)
      real(8) :: elemJac, elemVol
      
      character(len=30) :: meshName, fileName, solnFileName
      
      real(8) :: domainRange(3,2), domainLen(3)
      logical :: nearBoundary
      
      CALL MPI_INIT(IERROR)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,IPROC_MPI,IERROR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NPROCS_MPI,IERROR)
      
!------------------------------------------------------------------------------      
      ! READ IN THE MESH, CONNECTIVITY
      CALL READGM(G0XYZ,IJK,NX,NELX,meshName)      
      ! set element defaults
      MDIM=3
      NDF=3
      NGAUSS=1
      NODE=4
      NFACES=4
      NEDGES=6
      NEQ=NX*NDF
      
      if(IPROC_MPI==0) then
         write(*,*) 'Mesh imported:'
         write(*,*) '# of elements: ', NELX
         write(*,*) '# of nodes: ', NX
      endif
      
      allocate(nodalPositions(NX,3))
      do iNode=1,NX
         nodalPositions(iNode,1:3)=G0XYZ((iNode-1)*3+1:iNode*3)
      enddo
      
      if(IPROC_MPI==0) write(*,'(A)',advance='no') ' identifying node-to-element connectivity... '
      CALL readNODE_ELEM_CONNEC('nodeElementConnect.dat',IELEMNO,NINDX,NX,successReadNodeConnect)
      if (.not.successReadNodeConnect) then
         if(IPROC_MPI==0) write(*,'(A)',advance='no') 'calculating... '
         CALL NODE_ELEM_CONNEC(IELEMNO,NINDX,IJK,NX,NELX,MNELX,MAXNODELEM)
         if(IPROC_MPI==0) write(*,*) 'done.'
      else
         if(IPROC_MPI==0) write(*,*) 'imported.'
      endif

      ! READ TIME STEPPING INFORMATION
      open(102,file='time.out')
      nstep=0
      times(:) = 0.d0
      do
         read(102,*,iostat=error) times(nstep+1)
         if(error.EQ.-1) EXIT
         nstep=nstep+1
      enddo
      close(102)
      if (nstep.EQ.0) then
         if(IPROC_MPI==0) write(*,*) 'Cannot find file: time.out'
         STOP
      endif
      
      if(IPROC_MPI==0) &
         write(*,*) '# of time steps in time.out:', nstep
      
      ! import CI program options
      !call readOptionsCI('optionsCI.inp')
      !if (IPROC_MPI==0) CALL summarizeOptionsCI()
      ! assume these values:
      commentsInDataFile = .true.
      nStepsToAnalyse = 0
      
      ! determine the # of time steps to analyse
      nstepMax = nStepsToAnalyse
      if (nStepsToAnalyse==0) nstepMax = nstep
      if (nStepsToAnalyse > nstep) nstepMax = nstep

      CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
      
   !------------------------------------------------------------------------------

      ! ---------- IMPORT THE SOLUTION ----------!
      ! -------- read node displacements --------!
      allocate(elemNye(9,NELX))
      elemNye(:,:)=0.D0
      
      ! ----------- OPEN SOLN FILE: FP ---------------!
      INQUIRE(file='fp.out',exist=fileExists)
      if(.NOT.fileExists) then
         if(IPROC_MPI==0) write(*,*) 'Error reading the Gauss pt values of FP (fp.out)'
         stop
      endif
      open(103,file='fp.out')
      
      if(IPROC_MPI==0) then
         open(104,file='conv_nye.out')   ! export
         open(105,file='conv_nye_norm.out')   ! export
      endif
      ! ---------------------- IMPORT GRAINS ------------------------!
      ! first try grains.inp, if fails, elementGrains.inp, if fails create a trivial grain.
      CALL readGrains('grains.inp',NX,NELX)
      if(.NOT.grainsImported) then
           CALL readGrainsFromElementList('elementGrains.inp',NX,NELX)
      endif     
      if(.NOT.grainsImported) CALL createSingleGrain(NX,NELX)

      if (grainsImported.AND.IPROC_MPI.EQ.0) write(*,'(I0,A)') nGrains,' grains imported.'   ! KOUROSH_MPI
      
      ! CALCULATE ELEMENT-ELEMENT CONNECTIVITY
      allocate(elemFaceNeighbors(NELX,NFACES))
      if (IPROC_MPI==0) write(*,'(A)',advance='no') ' identifying element-element connectivity... '
      CALL readFaceNeighborsTET('elementNeighbors.dat',elemFaceNeighbors,NELX,successReadNeighbors)
      if (.not.successReadNeighbors) then
         CALL calcFaceNeighborsTET(elemFaceNeighbors(:,:),IJK,IELEMNO,NINDX, & 
                                   NELX,NODE,MAXNODELEM,MAXNODE)
         if (IPROC_MPI==0)write(*,*) 'done'
      else
         if (IPROC_MPI==0)write(*,*) 'imported'
      endif

      ! IDENTIFY THE NODES IN INDIVIDUAL GRAINS
      CALL identifyGrainNodes(NX,NELX,IJK)
      
      ! IDENTIFY THE NODES AT 
      ! 1) DOMAIN BOUNDARY, DOMAIN CORNERS
      ! 2) GRAIN BOUNDARIES
      ! 3) GRAIN-GRAIN BOUNDARIES
      CALL identifyBoundaryNodes(G0XYZ,NELX,NX,IJK,elemFaceNeighbors,domainRange)
      domainLen(:)=domainRange(:,2)-domainRange(:,1)
      
      if (IPROC_MPI==0) write(*,*) 'Domain range detected - X', domainRange(1,:)
      if (IPROC_MPI==0) write(*,*) 'Domain range detected - Y', domainRange(2,:)
      if (IPROC_MPI==0) write(*,*) 'Domain range detected - Z', domainRange(3,:)
            
      !read grain sizes: grainSize(:) - Deniz - Sept 22
      CALL readGrainSizes('grainSizes.inp')
      if(grainSizesImported.and.IPROC_MPI.EQ.0) write(*,*) 'grain sizes imported'
      if(.NOT.grainSizesImported.and.IPROC_MPI==0) then
         write(*,*) 'WARNING ! no grain size information (grainSizes.inp) found. assuming default grain size of:',grainSize(1)
      endif
            
      ! grain phases: grainPhase(:) - Deniz - Sept 22
      CALL readGrainPhases('grainPhases.inp')
      if(grainPhasesImported.and.IPROC_MPI==0) write(*,*) 'grain phases imported'
      if(.NOT.grainPhasesImported.and.IPROC_MPI==0) then
         write(*,*) 'WARNING ! no grain phase information (grainPhases.inp) found. assuming default phase of:',grainPhase(1)
      endif
      
      ! read grain texture - Deniz - Sept 22 - readGrainTexture() called from the grains module
      CALL readTexture('grainTexture.inp')
      if(textureImported.and.IPROC_MPI.EQ.0) write(*,*) 'texture imported'
      if(.NOT.textureImported.and.IPROC_MPI==0) then
         write(*,*) 'WARNING ! no texture information (grainTexture.inp) found.'&
      &'            assuming single crystal, euler:',grainTexture(:,1)
      endif

      CALL initialize_GrainNodes(NX,NELX,NODE,nGrains,nTotGrainNodes)
      
      ! initialize parameters/arrays for the analysis      
      ! -------- allocate solution vectors -------!
      allocate(gauss_values_Fp(nDOF_Fp,NELX))
      allocate(nodal_values_Fp(nDOF_Fp,nTotGrainNodes))
      nodal_values_Fp(:,:) = 0.D0
      gauss_values_Fp(:,:) = 0.d0
            
      if (IPROC_MPI==0) write(*,*) 'Processing...'
      lastProgress=0
      do i=1,MIN(60,nstep)
         if (IPROC_MPI==0) write(*,'(A)',advance='no') '_'
      enddo
      if (IPROC_MPI==0) write(*,*) ''

      nstepInp = 0
      do istep=1,nstep
 
         ! import Fp for this time step
         error = 0
         if (commentsInDataFile) then
            read(103,*,iostat=error) Fp_N, Fp_t
         else
            Fp_N = istep
            Fp_t = times(istep)
         endif
         if(error.NE.0) then   !end of input file      
            write(*,*) 'End of input file.'
            exit 
         endif 
         do iElem = 1, NELX      
            read(103,*,iostat=error) (gauss_values_Fp(iDOF,iElem),iDOF=1,nDOF_Fp)
            if(error.NE.0) exit
         enddo
         if(error.NE.0) then      
            write(*,*) 'Corrupt input file: fp.out.'
            exit
         endif

         this_Time = Fp_t       ! current time of the solution
         this_N = Fp_N          ! current time step of the solution
      
         write(strNSTEP,'(1I10)') this_N

         ! ----------- recover nodal values of Nye ----------- !
         ! -----------         at each grain       ----------- !
         nodal_values_Fp = 0.d0      
      
         CALL PARTITION_CALCULATOR(staGrainIdx,endGrainIdx,nGrainIdx, &
                                   nGrains,NPROCS_MPI,IPROC_MPI) 
      
         CALL recoveryLSQ(nodal_values_Fp, &
                          gauss_values_Fp, &
                          nDOF_Fp, &
                          G0XYZ,IJK, &
                          grainElementIndx,grainElements,IJK_grainNodes, &
                          nGrainNodes,grainNodesIndices, &
                          info,.false.,&
                          staGrainIdx,&
                          endGrainIdx)
                          
         if (info /= 0) then
            write(*,*) 'nodal recovery failed'
            CALL MPI_ABORT(IERROR)
         endif

         CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)      
         call MPI_ALLREDUCE(MPI_IN_PLACE,nodal_values_Fp,nDOF_Fp*nTotGrainNodes,&
            MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR)
         
         if (IPROC_MPI==0) then
            ! calculate the Nye tensor from nodal FP values
         
            ! we already have the nodal_values_Fp(:,:)
            do iElem=1,NELX
               !get the grain ID of the element
               iGrain=grainIDelem(iElem)
               CALL getElemXNodes(iElem,IJK,G0XYZ,NODE,xNodes)
               CALL calcJacTET4(xNodes, elemJacMat, elemJac, elemVol)   
               CALL MATINV3_UEL(elemJacMat,elemJacInv,elemVol)
               ! store the nodal values of transpose of Fp
               CALL getElemNodalValues(iElem,elemNodalNyeArr, &
                                     nodal_values_Fp,nDOF_Fp, &
                                         IJK_grainNodes,info)
                                       
               do iNode=1,4
                  elemNodalFp_T(1:3,1,iNode) = elemNodalNyeArr(1:3,iNode)
                  elemNodalFp_T(1:3,2,iNode) = elemNodalNyeArr(4:6,iNode) 
                  elemNodalFp_T(1:3,3,iNode) = elemNodalNyeArr(7:9,iNode) 
               enddo
               CALL calcGradTensorTET(elemJacInv,elemNodalFp_T, & 
                                                      gradFp_T)
               ! Nye = - (curl x (Fp)^T)^T
               CALL calcCurlTensor(-1.0*gradFp_T(:,:,:),Nye)
               Nye(:,:) = TRANSPOSE(Nye(:,:))
               
               ! save to array for further processing (GNDs, error indicators etc)
               elemNye(1:3,iElem)=Nye(1,1:3) ! ROW-MAJOR
               elemNye(4:6,iElem)=Nye(2,1:3) ! ROW-MAJOR
               elemNye(7:9,iElem)=Nye(3,1:3) ! ROW-MAJOR
            enddo
   ! ------------ NYE TENSOR CALCULATED: elemNye ------------ !
            
            ! EXPORT NYE TENSOR FIELD
            if (commentsInDataFile) &
               write(104,*,iostat=error) Fp_N, Fp_t
            do iElem = 1, NELX      
               write(104,'(9ES12.4)',iostat=error) (elemNye(iDOF,iElem),iDOF=1,9)
            enddo
            
            ! EXPORT the NORM
            if (commentsInDataFile) &
               write(105,*,iostat=error) Fp_N, Fp_t
            do iElem = 1, NELX         
               write(105,'(ES12.4)',iostat=error) DSQRT(DOT_PRODUCT(elemNye(:,iElem),elemNye(:,iElem)))
            enddo
            
         endif
         
         !progress bar fancyness
         if(IPROC_MPI==0) then
            if((60*istep)/nstep.GT.lastProgress) then
               write(*,'(A)',advance='no') '*'
               lastProgress=(60*istep)/nstep
            endif
         endif
      enddo
      
      if(IPROC_MPI==0) write(*,*) 'Total # of steps processed:', nstep
       
      write(*,*) ' '
      
      close(103)      ! close solution file for Fp
      if (IPROC_MPI==0) then
         close(104)      ! close output file for Nye
         close(105)      ! close output file for Nye Norm
      endif
            
      ! deallocate memory
      
      CALL grains_Destruct()
      
      if(allocated(gauss_values_Fp)) deallocate(gauss_values_Fp)
      if(allocated(nodal_values_Fp)) deallocate(nodal_values_Fp)
      !deallocate SPR arrays
      if(allocated(elemFaceNeighbors)) deallocate(elemFaceNeighbors)
      
      CALL MPI_FINALIZE(IERROR)
                  
      END



      SUBROUTINE READSTR(LR)
      IMPLICIT REAL*8(A-H,O-Z)
      !old version:
      !DIMENSION FLAG(20)
      !READ(LR,505) (FLAG(I),I=1,20)
!505   FORMAT(20A4)

! **** modified by deniz
      ! this code skips 1 line of text in the input file.
      character(150)::dummyString
      read(LR,'(A)') dummyString

      RETURN
      END
!*************************************************

      SUBROUTINE SHAPE_TET4(XEL,IPT,SHP,DVOL,XC)
      IMPLICIT NONE

      REAL(8):: SHP1(3,4),SHP(3,4),XJAC(3,3),XJACINV(3,3)
      REAL(8):: XEL(3,4),XC(3)
      REAL(8):: X(4),Y(4),Z(4),XX(3),YY(3),ZZ(3)
      REAL(8):: DVOL
      INTEGER:: I,J, INODE, IPT

      DO I=1,3
      DO J=1,4
      SHP1(I,J)=0.D0
      SHP(I,J)=0.D0
      ENDDO
      ENDDO

      SHP1(1,1)=1.D0

      SHP1(2,2)=1.D0

      SHP1(3,3)=1.D0


      SHP1(1,4)=-1.D0

      SHP1(2,4)=-1.D0

      SHP1(3,4)=-1.D0


   !		WRITE(*,*) 'XEL   IS:  ', XEL

      DO I=1,4
         X(I)=XEL(1,I)
         Y(I)=XEL(2,I)
         Z(I)=XEL(3,I)
      ENDDO

      DO I=1,3
         XX(I)=XEL(1,I)-XEL(1,4)
         YY(I)=XEL(2,I)-XEL(2,4)
         ZZ(I)=XEL(3,I)-XEL(3,4)
      ENDDO

      XJAC(1:3,1)=XX(1:3)
      XJAC(1:3,2)=YY(1:3)
      XJAC(1:3,3)=ZZ(1:3)



      CALL MATINV3_UEL(XJAC,XJACINV,DVOL)
      IF(DVOL.LT.1.D-16)THEN
         WRITE(*,*)'NEGATIVE VOLUME'
         STOP
      ENDIF


      DVOL=1.D0/6.D0*DABS(DVOL)

      DO INODE=1,4
         SHP(:,INODE)=MATMUL(XJACINV,SHP1(:,INODE))
      ENDDO



      RETURN
      END


!------------------------------------------------
      SUBROUTINE MATINV3_UEL(A,AI,DET)
      IMPLICIT NONE
      REAL(8), intent(in) :: A(3,3)
      REAL(8), intent(out):: AI(3,3)
      REAL(8), intent(out):: DET

      DET=(A(1,1)*A(2,2)*A(3,3)-A(1,1)*A(2,3)*A(3,2)-A(2,1)*A(1,2)  &
            *A(3,3)+A(2,1)*A(1,3)*A(3,2)+A(3,1)*A(1,2)*A(2,3)-A(3,1)  &
           *A(1,3)*A(2,2))

      AI(1,1) =  ( A(2,2)*A(3,3)-A(2,3)*A(3,2))/DET
      AI(1,2) = -( A(1,2)*A(3,3)-A(1,3)*A(3,2))/DET
      AI(1,3) = -(-A(1,2)*A(2,3)+A(1,3)*A(2,2))/DET
      AI(2,1) = -( A(2,1)*A(3,3)-A(2,3)*A(3,1))/DET
      AI(2,2) =  ( A(1,1)*A(3,3)-A(1,3)*A(3,1))/DET
      AI(2,3) = -( A(1,1)*A(2,3)-A(1,3)*A(2,1))/DET
      AI(3,1) =  ( A(2,1)*A(3,2)-A(2,2)*A(3,1))/DET
      AI(3,2) = -( A(1,1)*A(3,2)-A(1,2)*A(3,1))/DET
      AI(3,3) =  ( A(1,1)*A(2,2)-A(1,2)*A(2,1))/DET
      RETURN
      END

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
      character(len=50) :: fileName
      !dummies
      integer :: icyc,len,i,ii,j,jj,nnum,iElem,iNode,NEL,NEQ
      
      open(201,file='loadtime.inp')
      call READSTR(201)
      read(201,*) icyc
      call READSTR(201)
      read(201,*) meshName
      close(201)
      
      len=len_trim(meshName)
      fileName=meshName(1:len)//'.inp'

      open(unit=LR,file=fileName)

      call READSTR(LR)
      read(LR,*)MDIM,NDFADD
      NDF=MDIM+NDFADD     
      
      call READSTR(LR)
      read(LR,*) NX
      NEQ=NDF*NX
      
      IF(NDF.GT.NDOFELX)THEN
         WRITE(*,*)'INSUFFICIENT MEM-NDOFELX'
         STOP
      ENDIF

      !error check - deniz
      IF(NX.GT.MAXNODE)THEN
         write(*,*) 'Increase MAXNODE to ', NX
         WRITE(*,*)' MAXNODE is',MAXNODE
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
      
      call READSTR(LR)
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
      if(LEN(TRIM(matName)).GT.0)then
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
      if(LEN(TRIM(vecName)).GT.0)then
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
      if(LEN(TRIM(vecName)).GT.0)then
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
      if(LEN(TRIM(matName)).GT.0)then
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
      
       SUBROUTINE PARTITION_CALCULATOR(NELST,NELEND,NELNOX,N_TO_BE_DISTRIBUTED,NPROCS,IPROC)
       IMPLICIT NONE
       INTEGER::N_TO_BE_DISTRIBUTED,NPROCS,IPROC,I1,I2,IA,NELST,NELNOX,NELEND



       I1=N_TO_BE_DISTRIBUTED/NPROCS+1
       I2=N_TO_BE_DISTRIBUTED/NPROCS

       IA=N_TO_BE_DISTRIBUTED-I2*NPROCS


       IF((IPROC+1).LE.IA) THEN
              NELST=IPROC*I1+1
              NELNOX=I1
       ELSE
              NELST=IA*I1+(IPROC-IA)*I2+1
              NELNOX=I2
       ENDIF

       NELEND=NELST+NELNOX-1

       RETURN
       END      
