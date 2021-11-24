
      ! FBar Patch Module
      !----------------------
      !
      !   this module applies FBar-patch stabilization to TET4 elements.
      !     1) import/generate patches
      !     2) calculate FBar stabilizer multipliers and non-local gradient operators
      !
      !  Reference: J. Cheng, A. Shahba, S. Ghosh, Computational Mechanics 2016
      !
      !-----------------------
      ! Author: Deniz Ozturk
      !-----------------------
      module fBarPatch
      
      implicit none
            
      integer, private :: idMyProc, nProcs
      private :: partitionElements

      integer :: filePatch
      
      ! patch partitioning
      integer, allocatable :: staPatch_PROC(:)
      integer, allocatable :: endPatch_PROC(:)
      integer, allocatable :: nPatches_PROC(:)
      
      ! information about element partitioning
      integer, private, allocatable :: staElem_PROC(:)
      integer, private, allocatable :: endElem_PROC(:)
      integer, private, allocatable :: nElems_PROC(:)
      ! basic info on the mesh
      integer, private :: nElements, nNodes, nDOFs, nNodesPerElement
      
      ! Patch data structures
      integer :: nPatches,MaxElemInPatch,MinElemInPatch
      integer, allocatable :: patchIDfromElemID(:)
      integer, allocatable :: NElementsInPatch(:)
      integer, allocatable :: ElementsInPatch(:,:)
      
      ! FBar Stabilization: 'repositories' for hydrostatic multipliers and gradient matrices
      real(8), allocatable :: ksi_local(:,:)    !(2,nelnox)
      real(8), allocatable :: G_local(:,:,:,:)  !(3,4,MaxElemInPatch,nelnox)
      real(8), allocatable :: VolTau_local(:,:) !(MaxElemInPatch,nelnox)
      integer, parameter   :: timeTau = 1, timeT = 2, time0 = 3

      
      logical :: initializedFBarPatches

      contains

      subroutine initializePatches(nElements_l,nNodes_l,nDOFs_l,nNodesPerElement_l, &
                                   idMyProc_l,nProcs_l, &
                                   staElem_PROC_l,endElem_PROC_l,nElems_PROC_l, &
                                   error)
      implicit none
      
      integer, intent(in) :: nElements_l,nNodes_l,nDOFs_l,nNodesPerElement_l
      integer, intent(in) :: idMyProc_l, nProcs_l
      integer, intent(in) :: staElem_PROC_l(0:nProcs_l-1)
      integer, intent(in) :: endElem_PROC_l(0:nProcs_l-1)
      integer, intent(in) :: nElems_PROC_l(0:nProcs_l-1)
      integer, intent(out):: error
      
      integer :: idProc, elemID, elemIdx, I, J
      integer :: patchID,patchElemID,patchElemIdx
      logical :: fileExists
      
      character(len=50) :: strFileName
      
      if (allocated(NElementsInPatch)) deallocate(NElementsInPatch)
      if (allocated(ElementsInPatch)) deallocate(ElementsInPatch)
      if (allocated(patchIDfromElemID)) deallocate(patchIDfromElemID)
      initializedFBarPatches = .false.
      error = 0
      
      if (idMyProc_l > nProcs_l .or. idMyProc < 0 .or. nProcs_l == 0) then
         error = -100
         return
      endif
      
      nElements = nElements_l
      nNodes = nNodes_l
      nDOFs = nDOFs_l
      nNodesPerElement = nNodesPerElement_l
      
      ! read in patches
      strFileName = 'PATCH.inp'
      inquire(file=trim(strFileName),exist=fileExists)
      if(.NOT.fileExists) then
         error = -200
         return
      endif
      !open(newUnit = filePatch, file = trim(strFileName))
      filePatch = 7487
      open(filePatch,FILE='PATCH.inp')
      read(filePatch,*,iostat=error) nPatches,MaxElemInPatch,MinElemInPatch
      if (error /= 0) return
      
      allocate(NElementsInPatch(nPatches))
      allocate(ElementsInPatch(nPatches,MaxElemInPatch))
      
      NElementsInPatch=0
      ElementsInPatch=0
      
      do I=1,nPatches
         read(filePatch,*,iostat=error) NElementsInPatch(I), &
                     (ElementsInPatch(I,J),J=1,NElementsInPatch(I))      
         if (error /= 0) return
      enddo
      
      if (sum(NElementsInPatch(:)) /= nElements) then
         write(*,*) 'total number of elements in patches is not equal to the number of elements in the mesh'
         stop
      endif
      
      close(filePatch)
      ! generate map: patchID <- elementID
      allocate(patchIDfromElemID(nElements))
      do patchID=1,nPatches
         do elemIdx=1,NElementsInPatch(patchID)
            elemID=ElementsInPatch(patchID,elemIdx)
            patchIDfromElemID(elemID) = patchID
         enddo
      enddo
      
      idMyProc = idMyProc_l
      nProcs = nProcs_l
      
      allocate(staPatch_PROC(0:nProcs-1))
      allocate(endPatch_PROC(0:nProcs-1))
      allocate(nPatches_PROC(0:nProcs-1))
      do idProc=0,nProcs-1
         CALL partitionElements(nPatches,nProcs,idProc, & 
                                staPatch_PROC(idProc), & 
                                endPatch_PROC(idProc), & 
                                nPatches_PROC(idProc))
      enddo   
      
      allocate(staElem_PROC(0:nProcs-1))
      allocate(endElem_PROC(0:nProcs-1))
      allocate(nElems_PROC(0:nProcs-1))
      staElem_PROC = staElem_PROC_l
      endElem_PROC = endElem_PROC_l
      nElems_PROC = nElems_PROC_l
      initializedFBarPatches = .true.
      error = 0

      end subroutine
      
      SUBROUTINE fBar_INIT_SIZE_SLU(IJK,IPROC,N_UPDATE,IN_ST,NNZ_LOC, &
                                    IELEMNO,NINDX)
      IMPLICIT none
      integer, intent(in) :: IJK(nNodesPerElement*nElements),IPROC,N_UPDATE,IN_ST
      integer, intent(out):: NNZ_LOC
      integer, intent(in) :: IELEMNO(:), NINDX(:)
      
      ! locals
      integer :: INCOMP(nNodes),INDX,I,J,K,L,LL,kNode
      integer :: NNO,ICOMP,ISTART,IEND,IELNO,IMATCH,NN
      integer :: patchID,patchElemIdx,patchElemID
      
      if (.not.initializedFBarPatches) then
         write(*,*) 'Called fBar_INIT_SIZE_SLU before fBar patches were initialized'
         stop
      endif

      INDX=1

      DO I=1,N_UPDATE

         NNO=(IN_ST+I-1)/nDOFs+1
         ICOMP=0

         IF(NNO.EQ.1)THEN
            ISTART=1
         ELSE
            ISTART=NINDX(NNO-1)
         ENDIF

         IEND=NINDX(NNO)-1
         DO  J=ISTART,IEND

            IELNO=IELEMNO(J)

            DO K=1,nNodesPerElement

               IMATCH=0
               NN=IJK((IELNO-1)*nNodesPerElement+K)

               DO LL=1,ICOMP
                 IF(NN.EQ.INCOMP(LL))THEN
                    IMATCH=1
                    GOTO 30
                 ENDIF
               ENDDO

   30          IF(IMATCH.EQ.0)THEN

                 ICOMP=ICOMP+1
                 INCOMP(ICOMP)=NN

               ENDIF

            ENDDO
              
            ! FBAR PATCH --------------------------------------------------------------------------!
            ! the forces on this node are coupled to the displacements of the nodes of the elements
            ! that are in the same patch with each element tied to this node
            patchID = patchIDfromElemID(IELNO)
            do patchElemIdx = 1,NElementsInPatch(patchID)
               patchElemID = ElementsInPatch(patchID,patchElemIdx)

               do K=1,nNodesPerElement   ! traverse the nodes connected to this nNodesPerElement through elements, mark them as coupled

                  IMATCH=0
                  kNode=IJK((patchElemID-1)*nNodesPerElement+K)

                  DO LL=1,ICOMP
                    IF(kNode.EQ.INCOMP(LL))THEN
                       IMATCH=1
                       GOTO 40
                    ENDIF
                  ENDDO

   40             IF(IMATCH.EQ.0)THEN

                    ICOMP=ICOMP+1
                    INCOMP(ICOMP)=kNode

                  ENDIF

               enddo
               
            enddo
            ! FBAR PATCH --------------------------------------------------------------------------!
           
         ENDDO

         DO J=1,nNodes
            DO K=1,ICOMP
             IF(J.EQ.INCOMP(K))THEN
                DO L=1,nDOFs
                   INDX=INDX+1
                ENDDO
             ENDIF
            ENDDO
         ENDDO


      ENDDO

      NNZ_LOC=INDX-1

      END SUBROUTINE
      
      SUBROUTINE fBar_INIT_MSR_SLU(IJK,IPROC,N_UPDATE,IN_ST,IROWPTR, & 
                                   IBINDX,NNZ_LOC, &
                                   IELEMNO,NINDX)

      IMPLICIT none
      
      integer, intent(in) :: IJK(nNodesPerElement*nElements),IPROC,IN_ST,N_UPDATE, NNZ_LOC
      INTEGER, intent(out) :: IROWPTR(N_UPDATE+1)
      integer, intent(out):: IBINDX(NNZ_LOC)
      integer, intent(in) :: IELEMNO(:), NINDX(:)
      
!***********************************************************************
! C     IJK(I)
! C     N_UPDATE(I)
! C     IUPDATE(I)-CONTAINS ROWS TO BE UPDATED BY THIS PROCESSOR
! C     IROWPTR(O) - stores the data-index of first elements of each row. 
! C                  IROWPTR(:) = index of 1st el ol 1st row, index of 1st el of 2nd row... 
! C                  usage: k = IROWPTR(i) to IROWPTR(i+1)-1, VAL(k) traverses elements on row i
! C     IBINDX(O)-ARRAY CONTAINING NON-ZERO COLUMN INDICES FOR EACH ROW
! C     CONSISTENT WITH AZTEC
! C     VAL(O)-ARRAY WITH ELEMENTS INITIALIZED TO ZERO
      ! IBINDX has all the non-zero matrix coefficients for the equations of this CPU.
      ! IBINDX is ordered. (all DOFS coupled to the 1st DOF of this CPU, all DOFS coupled to the 2nd DOF of this CPU, all DOFS coupled to the 3rd DOF of this CPU, ... )
      ! INDX: total # of couplings from this CPU, size of IBINDX
!***********************************************************************
      ! locals
      INTEGER:: INCOMP(nNodes)
      integer:: IMATCH,IEND,IELNO,kNode
      INTEGER:: ICOMP,NNO,K
      INTEGER:: ISTART, NN,INDX,LL,I,J,L
      integer :: patchID,patchElemIdx,patchElemID
      
      if (.not.initializedFBarPatches) then
         write(*,*) 'Called fBar_INIT_MSR_SLU before fBar patches were initialized'
         stop
      endif

      INDX=1

      DO I=1,N_UPDATE         ! traverse all the nodes/DOFs handled by this processor

         NNO=(IN_ST+I-1)/nDOFs+1
         ICOMP=0
         IROWPTR(I)=INDX         !INDX = index of the coefficient in the 1-D data array

          IF(NNO.EQ.1)THEN
            ISTART=1
          ELSE
            ISTART=NINDX(NNO-1)  ! NINDX(i): the index of connectivity data of node a, in IELEMNO. IELEMNO(NINDX(i) + ...) = connected elemIDs
          ENDIF
          IEND=NINDX(NNO)-1

          DO  J=ISTART,IEND   ! traverse all elements connected to the node NNO. DOFs (I) of the node NNO are handled by this processor.
                              ! will now the node coupling information to determine positions of non-zero coefficients

              IELNO=IELEMNO(J) 

              DO K=1,nNodesPerElement

                 IMATCH=0
                 NN=IJK((IELNO-1)*nNodesPerElement+K)    ! found coupling: node NN is coupled to node NNO, (through element IELNO)

                 DO LL=1,ICOMP               ! check if already added to the list of couplings for this DOF (IN_ST+I)
                    IF(NN.EQ.INCOMP(LL))THEN
                       IMATCH=1
                       GOTO 30
                    ENDIF
                 ENDDO

  30             IF(IMATCH.EQ.0)THEN

                    ICOMP=ICOMP+1
                    INCOMP(ICOMP)=NN         ! add the coupled node to the list of couplings for this DOF (IN_ST+I)

                 ENDIF

             ENDDO
             

            ! FBAR PATCH --------------------------------------------------------------------------!
            ! the forces on this node are coupled to the displacements of the nodes of the elements
            ! that are in the same patch with each element tied to this node
            patchID = patchIDfromElemID(IELNO)
            do patchElemIdx = 1,NElementsInPatch(patchID)
               patchElemID = ElementsInPatch(patchID,patchElemIdx)

               do K=1,nNodesPerElement   ! traverse the nodes connected to this node through elements, mark them as coupled

                  IMATCH=0
                  kNode=IJK((patchElemID-1)*nNodesPerElement+K)

                  DO LL=1,ICOMP
                    IF(kNode.EQ.INCOMP(LL))THEN
                       IMATCH=1
                       GOTO 40
                    ENDIF
                  ENDDO

  40              IF(IMATCH.EQ.0)THEN

                    ICOMP=ICOMP+1
                    INCOMP(ICOMP)=kNode

                  ENDIF

               enddo
               
            enddo
            ! FBAR PATCH --------------------------------------------------------------------------!
             
          ENDDO
          
          ! ICOMP is the total # of nodes coupled to this DOF (IN_ST+I)
          ! INCOMP(:) is the list of nodes coupled to this DOF (IN_ST+I)
          
          ! traverse all the DOFs (nNodes -> (J-1)*nDOFs+L) and for the ones that are coupled to this DOF (INCOMP), mark the indices of non-zero coefficients in K(:,:) on IBINDX(:)
          DO J=1,nNodes
             DO K=1,ICOMP
                IF(J.EQ.INCOMP(K))THEN ! coupled node. mark all the DOFs of this node. store their EQN #s in order into IBINDX
                   DO L=1,nDOFs
                      IBINDX(INDX)=(J-1)*nDOFs+L
                      INDX=INDX+1
                   ENDDO
                ENDIF
             ENDDO
          ENDDO
          
          ! now all DOFs coupled to the DOF (IN_ST+I-1), are added to the end of IBINDX

      ENDDO
      
      IROWPTR(N_UPDATE+1)=INDX
      ! deniz - modified this line for a more conventional indexing structure
      ! OLD: IROWPTR(N_UPDATE+1)=INDX-1
      ! NEW: IROWPTR(N_UPDATE+1)=INDX
      ! now we can simply use:
      ! IS=IROWPTR(I)
      ! IE=IROWPTR(I+1)-1
      ! without the correction: IF(I.EQ.N_UPDATE)IE=IROWPTR(I+1)..
      
      END SUBROUTINE
      
      
      SUBROUTINE partitionElements(nElements,NPROCS,PROCID, & 
                                   NELST,NELEND,NELNOX)
      implicit none
      integer, intent(in) :: nElements, NPROCS, PROCID
      integer, intent(out) :: NELST,NELEND,NELNOX
      
      integer :: I1, I2, IA

      if (nElements < NPROCS) then ! each proc gets one.
         if (PROCID+1 <= nElements) then
            NELST = PROCID+1
            NELEND = PROCID+1
            NELNOX = 1
         else
            NELST = 0
            NELEND = 0
            NELNOX = 0
         endif
         return
      endif
      
      I1=nElements/NPROCS+1
      I2=nElements/NPROCS

      IA=nElements-I2*NPROCS

      IF((PROCID+1).LE.IA)THEN
         NELST=PROCID*I1+1
         NELNOX=I1
      ELSE
         NELST=IA*I1+(PROCID-IA)*I2+1
         NELNOX=I2
      ENDIF

      ! range of the elements handled by each processor: NELST -- NELEND
      NELEND=NELST+NELNOX-1
      END SUBROUTINE

      ! parallel calculation of FBar stabilizer factors, gradient matrices:
      ! ksi_local,G_local,VolTau_local
      !--------------------------------------------------
      subroutine calcFBarStabilizers(uT,uTau, & ! displacements over cycle, time increment info
                  IJK,G0XYZ,                  & ! mesh,parallelism etc.
                  calcStiffnessMatrix,success_out)
      implicit none
      include 'mpif.h'
      ! arguments
      real(8), intent(in) :: uT(nNodes*nDOFs),uTau(nNodes*nDOFs)
      integer, intent(in) :: ijk(nNodesPerElement*nElements)
      real(8), intent(in) :: g0xyz(nDOFs*nNodes)
      logical, intent(in) :: calcStiffnessMatrix
      logical, intent(out) :: success_out
      ! locals
      integer, parameter :: IPT=1 ! # of gauss pts
      integer :: totalVolumeRVE0
      REAL(8) ::XBAR_ST(6,12,1),XG_ST(9,12,1), &
                DVOLTAU,DVOL0,DVOLT, &
                DFGRD0(3,3),DFGRD1(3,3), &
                SHPTAU_ST(3,4)

      integer :: i,j,k,ii,nn,IERROR
      integer :: patchID,patchLocIdx,elemIdx,elemLocIdx,elemID,nodeID
      integer :: iProc,iNode
      real(8) :: uT_Element(nDOFs,nNodesPerElement)
      real(8) :: uTau_Element(nDOFs,nNodesPerElement)
      real(8) :: du_Element(nDOFs,nNodesPerElement)
      real(8) :: coords_Element(nDOFs,nNodesPerElement)
      real(8) :: patchVolumes(3,nPatches),patchVolume(3)
      real(8), allocatable :: detF_loc_temp(:,:,:)  ! (2,MaxElemInPatch,nPatches_PROC)
      real(8), allocatable :: detF_all(:,:,:)       ! (2,MaxElemInPatch,nPatches) only in master proc
      real(8), allocatable :: ksi_all(:,:)        ! (2,nElements) only in master proc
      real(8), allocatable :: G_local_temp(:,:,:,:)
      real(8), allocatable :: G_all(:,:,:,:)
      real(8), allocatable :: VolTau_loc_temp(:,:) ! (MaxElemInPatch,nPatches_PROC)
      real(8), allocatable :: VolTau_all(:,:)
      !real(8) :: G_buffer(3*4*MaxElemInPatch),bufferTest(100)
      logical :: success
      integer :: istatus(MPI_STATUS_SIZE)
      integer :: sizeRepositories(0:nProcs-1)
      integer :: shiftRepositories(0:nProcs-1)
      real(8) :: DET_FT, DET_FTAU, DET_F(2)
      integer, parameter :: TAG_G_TRANSFER = 300, TAG_ksi_TRANSFER = 301

      ! implementation
      success_out = .true. ! any error in the operations below will set this to false
      
      if (.not.initializedFBarPatches) then
         success_out = .false.
         return
      endif
      
      ! loop over elements in 'my' patches, extract displacement fields at times t and tau
      ! calculate their G matrices and F_bar stabilizers
      allocate(detF_loc_temp(2,MaxElemInPatch,nPatches_PROC(idMyProc)))
      detF_loc_temp = 0.d0

      if (calcStiffnessMatrix) then
         allocate(G_local_temp(3,4,MaxElemInPatch,nPatches_PROC(idMyProc)))
         allocate(VolTau_loc_temp(MaxElemInPatch,nPatches_PROC(idMyProc)))
         G_local_temp = 0.d0
         VolTau_loc_temp = 0.d0
      endif
      
      patchVolumes(:,:) = 0.d0
      do patchID = staPatch_PROC(idMyProc), endPatch_PROC(idMyProc)
         patchLocIdx = patchID - staPatch_PROC(idMyProc) + 1
         do elemIdx = 1,NElementsInPatch(patchID)
            elemID = ElementsInPatch(patchID,elemIdx)
            
            ! uElement <-- uT (extract element nodal displacements)
            uT_Element = 0.d0
            uTau_Element = 0.d0

            do iNode=1,nNodesPerElement
               nodeID=ijk(nNodesPerElement*(elemID-1)+iNode) 
               uT_Element(1:3,iNode) = uT((nodeID-1)*nDOFs+1:(nodeID-1)*nDOFs+3)
               uTau_Element(1:3,iNode) = uTau((nodeID-1)*nDOFs+1:(nodeID-1)*nDOFs+3)
               coords_Element(1:3,iNode) = G0XYZ((nodeID-1)*nDOFs+1:(nodeID-1)*nDOFs+3)
            enddo
            
            du_Element = uTau_Element - uT_Element
            
            ! calculate gradient and strain matrices for the element
            CALL calcElemGMATandVol(XBAR_ST,XG_ST,DVOLTAU, &
                coords_Element,uTau_Element,du_Element,DFGRD0,DFGRD1, &
                DVOL0,DVOLT,DET_FTAU,DET_FT,SHPTAU_ST,success)
                
            if (.not.success) then
               success_out = .false.
               exit
            endif
                
            ! sum contributions to init and current volumes of patches
            patchVolumes(timeTau,patchID) = patchVolumes(timeTau,patchID) + DVOLTAU
            patchVolumes(timeT,patchID) = patchVolumes(timeT,patchID) + DVOLT
            patchVolumes(time0,patchID) = patchVolumes(time0,patchID) + DVOL0
            detF_loc_temp(timeTau,elemIdx,patchLocIdx) = DET_FTAU
            detF_loc_temp(timeT,elemIdx,patchLocIdx) = DET_FT
                           
            ! save G matrix and element volume to repositories
            if (calcStiffnessMatrix) then
               G_local_temp(:,:,elemIdx,patchLocIdx) = SHPTAU_ST(:,:)
               VolTau_loc_temp(elemIdx,patchLocIdx) = DVOLTAU
            endif
            
         enddo
         if (.not.success_out) exit
      enddo
      
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,success_out,1, & 
                         MPI_LOGICAL, & 
                         MPI_LAND,MPI_COMM_WORLD,IERROR)
                         
      if (.not.success_out) then            
         if (allocated(G_local_temp)) deallocate(G_local_temp)
         if (allocated(VolTau_loc_temp)) deallocate(VolTau_loc_temp)
         deallocate(detF_loc_temp)
         return
      endif

      ! allreduce patch volumes,
      ! and collect detF's in the master processor
      ! (master processor will calculate the hydrostatic multipliers)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,patchVolumes,3*nPatches, & 
                         MPI_DOUBLE_PRECISION, & 
                         MPI_SUM,MPI_COMM_WORLD,IERROR)
      if (idMyProc==0) allocate(detF_all(2,MaxElemInPatch,nPatches))
      if (idMyProc==0) allocate(ksi_all(2,nElements))
      sizeRepositories(0:nProcs-1) = 2*MaxElemInPatch*nPatches_PROC(0:nProcs-1)
      shiftRepositories(0:nProcs-1) = 2*MaxElemInPatch*(staPatch_PROC(0:nProcs-1)-1)
      CALL MPI_Gatherv(detF_loc_temp,sizeRepositories(idMyProc),MPI_DOUBLE_PRECISION, &
                       detF_all,     sizeRepositories(:),shiftRepositories(:),MPI_DOUBLE_PRECISION, &
                       0,MPI_COMM_WORLD,IERROR)

      ! now we have the patch volumes and individual element detFtau's
      ! we can calculate hydrostatic multipliers (CPU-cheap task)
      if (idMyProc==0) then
         totalVolumeRVE0 = 0.d0 ! debug/code verification 
         do patchID=1,nPatches
            do elemIdx = 1,NElementsInPatch(patchID)
               elemID = ElementsInPatch(patchID,elemIdx)
               
               patchVolume(:) = patchVolumes(:,patchID)
               DET_F(:) = detF_all(:,elemIdx,patchID)
               ! calculate hydrostatic multiplier (FBar)
               ksi_all(timeTau,elemID) &
                  = (patchVolume(timeTau)/patchVolume(time0)/DET_F(timeTau))**(1.d0/3.d0)
               ksi_all(timeT,elemID) &
                  = (patchVolume(timeT)/patchVolume(time0)/DET_F(timeT))**(1.d0/3.d0)
            enddo
            totalVolumeRVE0 = totalVolumeRVE0 + patchVolumes(time0,patchID)
         enddo
      endif
      
      ! if this is the first call, allocate the processor-based data structures that will store the FBar stabilizers
      if (.not.allocated(ksi_local)) allocate(ksi_local(2,nElems_PROC(idMyProc)))
      if (.not.allocated(G_local)) allocate(G_local(3,4,MaxElemInPatch,nElems_PROC(idMyProc)))
      if (.not.allocated(VolTau_local)) allocate(VolTau_local(MaxElemInPatch,nElems_PROC(idMyProc)))
      
      ! distribute hydrostatic multipliers of elements
      ! to their parent processors (according to UMAT/element partitioning)
      if (idMyProc==0) then ! send ksi_all --> ksi_local
         ! master proc - first update your own local array
         ksi_local(:,:) = ksi_all(:,staElem_PROC(0):endElem_PROC(0))
         ! send to others
         do iProc = 1,nProcs-1
            CALL MPI_Send(ksi_all(:,staElem_PROC(iProc):endElem_PROC(iProc)), &
                          2*nElems_PROC(iProc),MPI_DOUBLE_PRECISION, &
                          iProc,TAG_ksi_TRANSFER,MPI_COMM_WORLD,IERROR)
         enddo
      else
         CALL MPI_Recv(ksi_local(:,:), &
                       2*nElems_PROC(idMyProc),MPI_DOUBLE_PRECISION, &
                       0,TAG_ksi_TRANSFER,MPI_COMM_WORLD,ISTATUS,IERROR)
      endif

      CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
      if (idMyProc==0) deallocate(detF_all)
      if (idMyProc==0) deallocate(ksi_all)
      

      ! All calculations are done. At this final stage, 
      ! redistribute the information to the processors that will need them
      ! Gather patch-based G_matrix and element volume repositories in master proc
      ! then distribute them to the CPU's depending on UMAT/Element partitioning
      if (calcStiffnessMatrix) then
         ! -- (i) G_matrix --
         ! large global repository temporarily allocated in master proc (~400MB for 500K mesh)
         if (idMyProc==0) allocate(G_all(3,4,MaxElemInPatch,nPatches))
         sizeRepositories(0:nProcs-1) = 3*4*MaxElemInPatch*nPatches_PROC(0:nProcs-1)
         shiftRepositories(0:nProcs-1) = 3*4*MaxElemInPatch*(staPatch_PROC(0:nProcs-1)-1)
         CALL MPI_Gatherv(G_local_temp,sizeRepositories(idMyProc),MPI_DOUBLE_PRECISION, &
                          G_all,     sizeRepositories(:),shiftRepositories(:),MPI_DOUBLE_PRECISION, &
                          0,MPI_COMM_WORLD,IERROR)
         deallocate(G_local_temp)  ! done with this. all data is in G_all now

         if (idMyProc==0) then  ! send G_all --> G_local
            
            ! master proc - first update your own local repository
            do elemID = staElem_PROC(0), endElem_PROC(0)
               elemLocIdx = elemID - staElem_PROC(idMyProc) + 1
               patchID = patchIDfromElemID(elemID)
               G_local(:,:,:,elemLocIdx) = G_all(:,:,:,patchID)
            enddo
            ! now send to all other procs
            do iProc = 1, nProcs-1
               do elemID = staElem_PROC(iProc), endElem_PROC(iProc)
                  patchID = patchIDfromElemID(elemID)
                  CALL MPI_Send(G_all(:,:,:,patchID), &
                                3*4*MaxElemInPatch,MPI_DOUBLE_PRECISION, &
                                iProc,TAG_G_TRANSFER,MPI_COMM_WORLD,IERROR)
                  ! MPI_Send, returns as soon as the buffer is read (does not wait for recv. as in MPI_SSend)
               enddo
            enddo
            
         else                 ! receive G_local <-- G_all
            do elemID = staElem_PROC(idMyProc), endElem_PROC(idMyProc)
               elemLocIdx = elemID - staElem_PROC(idMyProc) + 1
               patchID = patchIDfromElemID(elemID) ! I'm receiving the G of the elems of this patch
               CALL MPI_Recv(G_local(:,:,:,elemLocIdx), &
                             3*4*MaxElemInPatch,MPI_DOUBLE_PRECISION, &
                             0,TAG_G_TRANSFER,MPI_COMM_WORLD,ISTATUS,IERROR)
            enddo
         endif
         CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
         if (idMyProc==0) deallocate(G_all) ! deallocate large global repository
         
         ! -- (ii) element volumes  --
         if (idMyProc==0) allocate(VolTau_all(MaxElemInPatch,nPatches))
         sizeRepositories(0:nProcs-1) = (1)*MaxElemInPatch*nPatches_PROC(0:nProcs-1)
         shiftRepositories(0:nProcs-1) = (1)*MaxElemInPatch*(staPatch_PROC(0:nProcs-1)-1)
         CALL MPI_Gatherv(VolTau_loc_temp,sizeRepositories(idMyProc),MPI_DOUBLE_PRECISION, &
                          VolTau_all,     sizeRepositories(:),shiftRepositories(:),MPI_DOUBLE_PRECISION, &
                          0,MPI_COMM_WORLD,IERROR)
         deallocate(VolTau_loc_temp)  ! done with this. all data is in VolTau_all now
         
         if (idMyProc==0) then  ! send VolTau_all --> VolTau_local
            
            ! master proc - first update your own local repository
            do elemID = staElem_PROC(0), endElem_PROC(0)
               elemLocIdx = elemID - staElem_PROC(idMyProc) + 1
               patchID = patchIDfromElemID(elemID)
               VolTau_local(:,elemLocIdx) = VolTau_all(:,patchID)
            enddo
            ! now send to all other procs
            do iProc = 1, nProcs-1
               do elemID = staElem_PROC(iProc), endElem_PROC(iProc)
                  patchID = patchIDfromElemID(elemID)
                  CALL MPI_Send(VolTau_all(:,patchID), &
                                (1)*MaxElemInPatch,MPI_DOUBLE_PRECISION, &
                                iProc,TAG_G_TRANSFER,MPI_COMM_WORLD,IERROR)
                  ! MPI_Send, returns as soon as the buffer is read (does not wait for recv. as in MPI_SSend)
               enddo
            enddo
            
         else                 ! receive VolTau_local <-- VolTau_all
            do elemID = staElem_PROC(idMyProc), endElem_PROC(idMyProc)
               elemLocIdx = elemID - staElem_PROC(idMyProc) + 1
               patchID = patchIDfromElemID(elemID) ! I'm receiving the G of the elems of this patch
               CALL MPI_Recv(VolTau_local(:,elemLocIdx), &
                             (1)*MaxElemInPatch,MPI_DOUBLE_PRECISION, &
                             0,TAG_G_TRANSFER,MPI_COMM_WORLD,ISTATUS,IERROR)
            enddo
         endif
         CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
         if (idMyProc==0) deallocate(VolTau_all) ! deallocate large global repository
         
      endif

      deallocate(detF_loc_temp)

      end subroutine
      
      
      end module
      
