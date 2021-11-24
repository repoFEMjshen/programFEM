      ! grains module
      !     analyses grains within an polycrystalline mesh and provides the following structures
      !     for use in CPFE codes and post-processing tools:
      !     grain orientation, size and phase 
      !     list of elements within each grain
      !     list of nodes within each grain
      !     grain-grain boundaries, and lists of nodes on those boundaries
      !     grain centroids, volumes
      !     Voronoi-tesellated areas of nodes at the grain-grain boundaries
      !     Lists of nodes with respect to domain:
      !        list of nodeIDs on each face and each corner of the domain
      !        spatial range of the domain in each dimension
      !     
      !  DEPENDENCIES:
      !  meshTools module
      
      MODULE grains

      integer, parameter :: MAXN_GRAINS=600 !3000
      integer, parameter :: MAXN_DOM_BND_NODES=15000 !15000
      integer, parameter :: MAXN_NODES_PER_GRAIN_BND=500 ! 6000
      integer, parameter :: MAXN_NODES_PER_GRGR_BND=200 !500
      integer, parameter :: MAX_GRAIN_NEIGHBORS=70 !100
      integer, parameter :: MAXN_GR_GR_BND=MAX_GRAIN_NEIGHBORS*MAXN_GRAINS
      integer, parameter :: MAX_GR_PER_NODE=10
      
      integer, allocatable :: grainElements(:)
      integer, allocatable :: grainIDelem(:)
      integer, allocatable :: grainIDnode(:)
      integer, allocatable :: grainListsForNodes(:,:)
      integer, allocatable :: nGrainListsForNodes(:)
      integer :: nGrainsOnBoundary
      logical, allocatable :: isGrainOnBoundary(:)
      
      real(8), allocatable :: grainCentroids(:,:)
      real(8), allocatable :: grainVolumes(:)
      
      integer, allocatable :: grainNodes(:)        ! list of nodes within each grain (compact array)
      integer, allocatable :: grainNodesIndices(:) ! indices for the starting positions of node lists for each grain. grainNodesIndices(iGrain-1)=start pos. for iGrain.
      integer, allocatable :: IJK_grainNodes(:,:)
      integer, allocatable :: nGrainNodes(:)
      integer :: nTotGrainNodes

      integer :: nGrains
      logical :: grainsImported,grainSizesImported
      logical :: grainPhasesImported,textureImported
      logical :: identifiedGrainNodes,identifiedBoundaryNodes
      real(8) :: grainSize(MAXN_GRAINS)
      integer :: grainPhase(MAXN_GRAINS)
      real(8) :: grainTexture(3,MAXN_GRAINS)
      integer :: nGrainElements(MAXN_GRAINS)
      integer :: grainElementIndx(MAXN_GRAINS)
      integer :: grainNeighborList(MAX_GRAIN_NEIGHBORS,MAXN_GRAINS)
      integer :: nGrainNeighbors(MAXN_GRAINS)
      
      
      integer :: domainBoundaryNodes(MAXN_DOM_BND_NODES)
      integer :: nDomainBoundaryNodes
      integer, allocatable :: grainBoundaryNodes(:)
      integer, allocatable :: nGrainBoundaryNodes(:)
      integer, allocatable :: grainBoundaryNodesIndx(:)
      integer :: gr_grBoundaryNodes(MAXN_GR_GR_BND, & 
                                    MAXN_NODES_PER_GRGR_BND)
      real(8) :: gr_grBoundaryNodeEffAreas(MAXN_GR_GR_BND, & 
                                           MAXN_NODES_PER_GRGR_BND)
      real(8) :: gr_grBoundaryAreas(MAXN_GR_GR_BND)
      integer :: nGrGrBoundaryNodes(MAXN_GR_GR_BND)
      integer :: gr_grBoundaryIDs(MAXN_GRAINS,MAXN_GRAINS)
      integer :: nGr_grBoundaries

      integer :: cornerNodes(8)
      integer :: faceNodes(6,MAXN_DOM_BND_NODES/6) ! stores nodeIDs given (faceID,nodeIndx)
      integer :: nFaceNodes(6)               ! nodes on each face of the cubic domain
      
      contains

      SUBROUTINE grains_Initialize(NX,NELX)
         implicit none
         integer, intent(in) :: NX,NELX
         grainsImported  = .FALSE.
         grainSizesImported = .FALSE.
         grainPhasesImported = .FALSE.
         textureImported = .FALSE.
         identifiedGrainNodes = .FALSE.
         identifiedBoundaryNodes = .False.
         allocate(grainElements(NELX))
         allocate(grainIDelem(NELX))
         allocate(grainIDnode(NX))
      END SUBROUTINE grains_Initialize
      
      SUBROUTINE grains_Destruct()
         implicit none
         if (allocated(grainCentroids)) deallocate(grainCentroids)
         if (allocated(grainVolumes)) deallocate(grainVolumes)
         if (allocated(grainElements)) deallocate(grainElements)
         if (allocated(grainIDelem)) deallocate(grainIDelem)
         if (allocated(grainIDnode)) deallocate(grainIDnode)
         if (allocated(grainNodes)) deallocate(grainNodes)
         if (allocated(grainNodesIndices)) deallocate(grainNodesIndices)
         if (allocated(grainListsForNodes)) deallocate(grainListsForNodes)
         if (allocated(nGrainListsForNodes)) deallocate(nGrainListsForNodes)
         if (allocated(nGrainNodes)) deallocate(nGrainNodes)
         if (allocated(isGrainOnBoundary)) deallocate(isGrainOnBoundary)
         if (allocated(IJK_grainNodes)) deallocate(IJK_grainNodes)
         if (allocated(grainBoundaryNodes)) deallocate(grainBoundaryNodes)
         if (allocated(nGrainBoundaryNodes)) deallocate(nGrainBoundaryNodes)
         if (allocated(grainBoundaryNodesIndx)) deallocate(grainBoundaryNodesIndx)
      END SUBROUTINE grains_Destruct

      SUBROUTINE createSingleGrain(NX,NELX)
      implicit none
      integer, intent(in)  :: NX,NELX
      ! OUT:
      !  grainIDelem(NELX)
      !  nGrains
      !  grainElements(NELX)
      !  nGrainElements(1)
      !  grainElementIndx(1)

      ! locals
      integer :: lastIndex, iIndex, staIndex, endIndex
      integer :: ElemID
      nGrains=1
      
      CALL grains_Initialize(NX,NELX)
      
      do ElemID=1,NELX
         grainIDelem(ElemID)=1
         grainElements(ElemID)=ElemID
      enddo
      nGrainElements(1)=NELX
      grainElementIndx(1)=NELX+1
      END SUBROUTINE
      
      ! in:
      !  grainElementIndx
      ! out: 
      !  nGrains
      !  grainIDelem(NELX)
      !  grainElements(NELX): grain1: elemID elemID elemID | grain2: elemID elemID ...
      !  grainElementIndx(MAXN_GRAINS) : indices for grainElements() indicating where the element list for each grain starts
      !  nGrainElements(MAXN_GRAINS)
      SUBROUTINE readGrains(strGrainFile,NX,NELX)
      implicit none
      
      character(len=*), intent(in) :: strGrainFile
      integer, intent(in)  :: NX,NELX

      !locals
      integer :: nElements,I,J,N,ISTART,IEND,INC,error
      integer :: lastIndex, iIndex, staIndex, endIndex
      integer :: ElemID
      logical :: fileExists
      
      nGrains = 0
      grainsImported = .FALSE.
      
      inquire(file=strGrainFile,exist=fileExists)
      if(.NOT.fileExists) RETURN

      open(51,file=trim(strGrainFile))
      READ(51,*,iostat=error) nGrains
      if(error.NE.0) goto 417

      IF(nGrains.GT.MAXN_GRAINS)THEN
         WRITE(*,*)'readGrains(): increase MAXN_GRAINS to ',nGrains
         STOP
      ENDIF
      
      CALL grains_Initialize(NX,NELX)

      ! read the element list for each grain
      lastIndex=0
      grainIDelem(:)=0
      DO I=1,nGrains    ! read in the list of element for each grain
         read(51,*,iostat=error) nElements
         read(51,*,iostat=error) (grainElements(lastIndex+J),J=1,nElements)
         if(error.NE.0) goto 417
         lastIndex=lastIndex+nElements
         grainElementIndx(I)=lastIndex+1
         nGrainElements(I)=nElements
      ENDDO
      
      ! construct element-to-grain connectivity
      staIndex = 1
      do I=1,nGrains
         if (I.NE.1) staIndex = grainElementIndx(I-1)
         endIndex = grainElementIndx(I)-1
         do iIndex=staIndex,endIndex
            ElemID = grainElements(iIndex)
            grainIDelem(ElemID) = I
         enddo
      enddo
      
      if (nGrains.NE.0) grainsImported = .TRUE.
      
      close(51)
      
      RETURN
      !ERROR HANDLING
417      WRITE(*,*) 'grains.inp -- corrupt data file.'
         STOP
      END SUBROUTINE readGrains
      
      SUBROUTINE readGrainsFromElementList(strGrainFile,NX,NELX)
!     constructs grainElements(:)
!	   this array stores the grain->element connection info
!	   g1e1 g1e2 g1e3 g2e1 g2e2 g2e3 g2e4 .. .. .. g3e1 ...		   gxex  0 0 0
!	   1	            ^grainElementIndx(1) 	     ^grainElementIndx(2) 	^grainElementIndx(nGrains)
      implicit none

      character(len=*), intent(in) :: strGrainFile
      integer, intent(in)  :: NX,NELX
      ! OUT:
      !  grainIDelem(NELX)
      !  nGrains
      !  grainElements(NELX)
      !  nGrainElements(MAXN_GRAINS)
      !  grainElementIndx(MAXN_GRAINS)
      !locals
      integer :: I,J,ISTART,IEND,INC
      integer :: lastIndex,iIndex,staIndex,endIndex
      integer :: elemID, nElemsInFile, myGrainID, error
      logical :: found
      integer :: grainList(MAXN_GRAINS)
      integer :: elementList(NELX)
      
      nGrains=0
      grainsImported = .FALSE.
      
      open(51,file=trim(strGrainFile))
      READ(51,*) nElemsInFile

      IF(nElemsInFile.GT.NELX)THEN
         WRITE(*,*) 'readGrainsFromElementList():'
         WRITE(*,*) 'number of elements in grain file is larger than the input file'
         STOP
      ENDIF
      
      CALL grains_Initialize(NX,NELX)
      
      nGrainElements(:)=0
      grainIDelem(:)=0
      ! read element-to-grain connectivity
      staIndex = 1
      do I=1,nElemsInFile
         READ(51,*,iostat=error) elemID, grainIDelem(elemID)
         elementList(I) = elemID
         if(error.NE.0) exit
         ! add the grainID to grain list
         found=.FALSE.
         do J=1,nGrains
            if(grainList(J).EQ.grainIDelem(elemID)) then
               found=.TRUE.
               exit
            endif
         enddo
         if(.NOT.found) then
            nGrains=nGrains+1
            grainList(nGrains)=grainIDelem(elemID)
         endif
         nGrainElements(grainIDelem(elemID))= &
               nGrainElements(grainIDelem(elemID)) + 1
         !complain, if too many grains
         IF(nGrains.GT.MAXN_GRAINS)THEN
            WRITE(*,*) 'readGrainsFromElementList(): INSUFFICIENT MEMORY-MAXN_GRAINS'
            STOP
         ENDIF
      enddo
      close(51)
      if(error.NE.0) then
         WRITE(*,*) 'elementGrains.inp -- corrupt data file.'
         STOP
      endif

      ! construct the element list for each grain
      lastIndex=0
      grainElementIndx(:)=0
      do myGrainID=1,nGrains
         do J=1,nElemsInFile
            elemID=elementList(J)
            if(grainIDelem(elemID).EQ.myGrainID)then
               lastIndex=lastIndex+1
               grainElements(lastIndex)=elemID
            endif
         enddo
         grainElementIndx(myGrainID)=lastIndex+1
      enddo
      
      if (nGrains.NE.0) grainsImported = .TRUE.

      RETURN
      END SUBROUTINE
      
!     reads texture from file and 
!     constructs grainTexture(3,:)
      SUBROUTINE readTexture(strTextureFile,fileExistsOut)
      implicit none

      character(len=*), intent(in) :: strTextureFile
      logical, optional, intent(out) :: fileExistsOut
      
      logical :: fileExists      
      
      !locals
      integer :: I,J
      integer :: lastIndex,iIndex,staIndex,endIndex
      integer :: elemID,grainID,nGrainsInFile, error
      
      grainTexture(:,:)=0.d0
      
      inquire(file=trim(strTextureFile),exist=fileExists)
      if(present(fileExistsOut)) fileExistsOut = fileExists
      if(.NOT.fileExists) then
         textureImported = .FALSE.
         return
      endif

      open(51,file=trim(strTextureFile))
      READ(51,*) nGrainsInFile

      IF(.not.grainsImported)THEN
         WRITE(*,*) 'readTexture() called before grains imported. call readGrains first.'         
         STOP
      ENDIF
      
      IF(nGrainsInFile.GT.nGrains)THEN
         WRITE(*,*) 'number of grains in texture file does not match the number of grains imported'
         STOP
      ENDIF
      
      ! read euler angles for each grain
      do I=1,nGrainsInFile
         READ(51,*,iostat=error) grainID, grainTexture(:,grainID)
         if(error.NE.0) exit
      enddo
      close(51)
      if(error.NE.0) then
         WRITE(*,*) trim(strTextureFile)//' -- corrupt data file.'
         STOP
      endif

      textureImported = .TRUE.

      RETURN
      END SUBROUTINE
      

      SUBROUTINE readTexture_direct(strTextureFile,nGrainsInFile,success)
      implicit none

      character(len=*), intent(in) :: strTextureFile
      logical, optional, intent(out) :: success
      integer, intent(out) :: nGrainsInFile
      
      logical :: fileExists      
      
      !locals
      integer :: I,J
      integer :: lastIndex,iIndex,staIndex,endIndex
      integer :: elemID,grainID, error
      
      grainTexture(:,:)=0.d0
      if(present(success)) success = .false.
      
      inquire(file=trim(strTextureFile),exist=fileExists)
      if(.NOT.fileExists) return

      open(51,file=trim(strTextureFile))
      read(51,*) nGrainsInFile
      
      ! read euler angles for each grain
      do I=1,nGrainsInFile
         READ(51,*,iostat=error) grainID, grainTexture(:,grainID)
         if(error.NE.0) exit
      enddo
      close(51)
      if(error.NE.0) then
         WRITE(*,*) trim(strTextureFile)//' -- corrupt data file.'
         return
      endif

      if(present(success)) success = .true.
      nGrains = nGrainsInFile ! update module nGrains 
      textureImported = .true.

      END SUBROUTINE
      
      SUBROUTINE readGrainPropsCMRL(strTextureFile,success)
      implicit none

      character(len=*), intent(in) :: strTextureFile
      logical, optional, intent(out) :: success
      
      !locals
      logical :: fileExists
      integer :: I,J
      integer :: lastIndex,iIndex,staIndex,endIndex
      integer :: elemID,grainID,nGrainsInFile,nGrainProps,error
      
      grainTexture(:,:)=0.d0
      grainPhase(:) = 0
      grainSize(:) = 0.d0
      if(present(success)) success = .false.
      
      inquire(file=trim(strTextureFile),exist=fileExists)
      if(.NOT.fileExists) return

      open(51,file=trim(strTextureFile))
      read(51,*) nGrainsInFile, nGrainProps
      
      if (nGrainProps /= 5) then
         WRITE(*,*) trim(strTextureFile)//' -- unexpected file format. nGrainProps /= 5'
         return
      endif
      
      ! read euler angles for each grain
      do I=1,nGrainsInFile
         READ(51,*,iostat=error) grainID, nGrainProps, grainPhase(grainID), grainTexture(:,grainID), grainSize(grainID)
         if(error.NE.0) exit
      enddo
      close(51)
      if(error.NE.0) then
         WRITE(*,*) trim(strTextureFile)//' -- corrupt data file. read error:',error
         return
      endif

      if(present(success)) then
         success = .true.
         nGrains = nGrainsInFile
         grainSizesImported = .true.
         textureImported = .true.
         grainPhasesImported = .true.
      endif

      END SUBROUTINE
	  
      ! IDENTIFY THE NODES IN INDIVIDUAL GRAINS
      SUBROUTINE identifyGrainNodes(NX,NELX,IJK)
      implicit none

      integer, intent(in) :: NX,NELX
      integer, intent(in) :: IJK(4*NELX)
      ! locals
      integer :: iElemID,iLocNode,iGloNode,I,J
      integer :: iGrainID, iNodeIndex, staIndex, endIndex
      integer :: NFACES,NNODES
      integer, allocatable :: nodeIDtoGrainNodeID_temp(:,:)
      integer, allocatable :: grainNodes_temp(:,:)
      
      NFACES = 4
      NNODES = 4
      
      allocate(grainNodes_temp(NX,nGrains))
      allocate(nodeIDtoGrainNodeID_temp(NX,nGrains))
      allocate(nGrainNodes(nGrains))
      allocate(grainListsForNodes(NX,MAX_GR_PER_NODE))
      allocate(nGrainListsForNodes(NX))
      allocate(IJK_grainNodes(4,NELX))

      nGrainNodes(:)=0
      grainNodes_temp(:,:)=0
      nodeIDtoGrainNodeID_temp(:,:)=0
      grainListsForNodes(:,:) = 0 
      nGrainListsForNodes(:) = 0
      identifiedGrainNodes = .false.
      
      do iElemID=1,NELX
         iGrainID=grainIDelem(iElemID)
         do iLocNode=1,4
            iGloNode=IJK(NNODES*(iElemID-1)+iLocNode)
            if(nodeIDtoGrainNodeID_temp(iGloNode,iGrainID)==0) then
               ! add node to the list of the grain
               iNodeIndex = nGrainNodes(iGrainID) + 1
               nGrainNodes(iGrainID) = iNodeIndex
               grainNodes_temp(iNodeIndex,iGrainID) = iGloNode
               ! add grain to the list of the node
               nGrainListsForNodes(iGloNode) = nGrainListsForNodes(iGloNode) + 1
               if(nGrainListsForNodes(iGloNode).GT.MAX_GR_PER_NODE) then
                  write(*,*) 'increase MAX_GR_PER_NODE'
                  stop
               endif
               grainListsForNodes(iGloNode,nGrainListsForNodes(iGloNode)) = iGrainID
               nodeIDtoGrainNodeID_temp(iGloNode,iGrainID)=iNodeIndex
            endif
         enddo
      enddo
      
      ! determine the size of the compact array (total number of grain nodes)
      nTotGrainNodes = 0
      do iGrainID=1,nGrains
         nTotGrainNodes = nTotGrainNodes + nGrainNodes(iGrainID)
      enddo
      ! create the compact array of grain nodes
      allocate(grainNodes(nTotGrainNodes))
      allocate(grainNodesIndices(nGrains))
      staIndex = 1
      do iGrainID=1,nGrains
         endIndex = staIndex + nGrainNodes(iGrainID) - 1
         grainNodes(staIndex:endIndex) = grainNodes_temp(1:nGrainNodes(iGrainID),iGrainID)
         grainNodesIndices(iGrainID) = endIndex + 1
         staIndex = staIndex + nGrainNodes(iGrainID)
      enddo
      
      ! construct the element->grainNodeID connectivity matrix
      do iElemID=1,NELX
         do iLocNode=1,4
            iGloNode = IJK(NNODES*(iElemID-1)+iLocNode)
            iGrainID = grainIDelem(iElemID)
            iNodeIndex = nodeIDtoGrainNodeID_temp(iGloNode,iGrainID)
            staIndex = 1
            if (iGrainID /= 1) staIndex = grainNodesIndices(iGrainID-1)
            IJK_grainNodes(iLocNode,iElemID) = staIndex - 1 + iNodeIndex
         enddo
      enddo
      
      ! deallocate large array
      deallocate(grainNodes_temp)
      deallocate(nodeIDtoGrainNodeID_temp)
      
      identifiedGrainNodes = .true.
      
      END SUBROUTINE
      
      SUBROUTINE getGrainIDFromGrainNodeID(grainNodeID,grainID)
      implicit none
      integer, intent(in) :: grainNodeID
      integer, intent(out):: grainID
      
      grainID = 0
      if (.not.identifiedGrainNodes) return
      do grainID=1,nGrains
         if (grainNodesIndices(grainID) > grainNodeID) exit
      enddo
      
      END SUBROUTINE

      SUBROUTINE readGrainSizes(strGrainFile,fileExistsOut)
         implicit none
         character(len=*), intent(in) :: strGrainFile
         logical, optional, intent(out) :: fileExistsOut
         
         logical :: fileExists
         integer :: iGrain, nEntries, I, error

         grainSize(:) = 0.D0
         grainSizesImported = .FALSE.

         inquire(file=strGrainFile,exist=fileExists)
         if(present(fileExistsOut)) fileExistsOut = fileExists
         if (fileExists) then
            open(51,file=strGrainFile)
            read(51,*) nEntries
            do I=1,nEntries
               read(51,*,iostat=error) iGrain,grainSize(iGrain)
               if (error > 0) then
                  write(*,*) 'file read error while reading grain sizes. error #',error, &
                             'file:',trim(strGrainFile)
               elseif( error < 0) then
                  write(*,*) 'end of grain size list before all grains are read. # of entries in file header:',nEntries,&
                             'file:',trim(strGrainFile)
               endif
               if(iGrain.GT.nGrains)then
                  write(*,*) 'unidentified grain ID ',iGrain, & 
                             ' in grainSizes.inp.'
                  write(*,*) '# of grains in model:',nGrains
                  STOP
               endif
            enddo
            close(51)
         else  ! assign default grain sizes
            grainSize(:) = 10.D0
         endif
         
         if(fileExists) grainSizesImported = .TRUE.
         
      END SUBROUTINE      
      
      SUBROUTINE readGrainSizes_direct(strGrainFile,nGrainsInFile,success)
      implicit none
      character(len=*), intent(in) :: strGrainFile
      logical, optional, intent(out) :: success
      integer, intent(out) :: nGrainsInFile
      
      logical :: fileExists
      integer :: iGrain, I, error

      grainSize(:) = 0.D0
      if(present(success)) success = .false.

      inquire(file=strGrainFile,exist=fileExists)
      if (.not.fileExists) return
   
      open(51,file=strGrainFile)
      read(51,*) nGrainsInFile
      
      do I=1,nGrainsInFile
         read(51,*,iostat=error) iGrain,grainSize(iGrain)
         if (error > 0) then
            write(*,*) 'file read error while reading grain sizes. error #',error, &
                       'file:',trim(strGrainFile)
         elseif( error < 0) then
            write(*,*) 'end of grain size list before all grains are read. # of entries in file header:',nGrainsInFile,&
                       'file:',trim(strGrainFile)
         endif
         if(iGrain.GT.nGrainsInFile)then
            write(*,*) 'unidentified grain ID ',iGrain, & 
                       ' in grainSizes.inp.'
            write(*,*) '# of grains in model:',nGrainsInFile
            STOP
         endif
      enddo
      close(51)
      
      if(present(success)) success = .true.
      grainSizesImported = .TRUE.
      nGrains = nGrainsInFile
      
      END SUBROUTINE      

      SUBROUTINE readGrainPhases(strGrainFile,fileExistsOut)
      implicit none
      character(len=*), intent(in) :: strGrainFile
      logical, optional, intent(out) :: fileExistsOut
      
      logical :: fileExists
      integer :: iGrain, nEntries, I, error

      ! assign default grain phase
      grainPhase(:) = 1
      grainPhasesImported = .FALSE.

      inquire(file=strGrainFile,exist=fileExists)
      if(present(fileExistsOut)) fileExistsOut = fileExists
      if (fileExists) then
         open(51,file=strGrainFile)
         read(51,*) nEntries
         do I=1,nEntries
            read(51,*,iostat=error) iGrain,grainPhase(iGrain)
            if (error > 0) then
               write(*,*) 'file read error while reading grain phases. error #',error, &
                          'file:',trim(strGrainFile)
            elseif( error < 0) then
               write(*,*) 'end of grain phase list before all grains are read. # of entries in file header:',nEntries,&
                          'file:',trim(strGrainFile)
            endif
            if(iGrain.GT.nGrains)then
               write(*,*) 'unidentified grain ID ',iGrain, & 
                          ' in',strGrainFile
               write(*,*) '# of grains in model:',nGrains
               STOP
            endif
         enddo
         close(51)
         
         grainPhasesImported = .TRUE.
      endif
      
      END SUBROUTINE
      

      SUBROUTINE readGrainPhases_direct(strGrainFile,nGrainsInFile,success)
      implicit none

      character(len=*), intent(in) :: strGrainFile
      integer, intent(out) :: nGrainsInFile
      logical, optional, intent(out) :: success
      
      logical :: fileExists
      integer :: iGrain, I, error

      ! assign default grain phase
      grainPhase(:) = 1
      if(present(success)) success = .false.

      inquire(file=strGrainFile,exist=fileExists)
      if (.not.fileExists) return
    
      if (.not.fileExists) return
   
      open(51,file=strGrainFile)
      read(51,*) nGrainsInFile
            
      do I=1,nGrainsInFile
         read(51,*,iostat=error) iGrain,grainPhase(iGrain)
         if (error > 0) then
            write(*,*) 'file read error while reading grain phases. error #',error, &
                       'file:',trim(strGrainFile)
         elseif( error < 0) then
            write(*,*) 'end of grain phase list before all grains are read. # of entries in file header:',nGrainsInFile,&
                       'file:',trim(strGrainFile)
         endif
         if(iGrain.GT.nGrainsInFile)then
            write(*,*) 'unidentified grain ID ',iGrain, & 
                       ' in',strGrainFile
            write(*,*) '# of grains in model:',nGrainsInFile
            STOP
         endif
      enddo
      close(51)
      
      if(present(success)) success = .true.
      grainPhasesImported = .TRUE.
      nGrains=nGrainsInFile
      
      END SUBROUTINE
      
      
      
      ! IDENTIFY THE NODES AT 
      ! 1) DOMAIN BOUNDARY
      ! 2) GRAIN BOUNDARIES
      ! 3) GRAIN-GRAIN BOUNDARIES
      ! ********************************************
      ! following arrays are constructed:
      ! grainIDnode(nodeID) =
      !     0: node on domain boundary.
      !   > 0: node inside a grain. returned value is the grain ID.
      !   < 0: node on a grain-grain boundary. returned value is the gr-gr boundary ID
      ! domainBoundaryNodes(nodeIndex) = nodeID
      !        stores the IDs of nodes that are on the domain boundary
      ! grainBoundaryNodes(grainID,nodeIndex) = nodeID
      !        stores the IDs of nodes on individual grain boundaries
      ! gr_grBoundaryIDs(grain1_ID,grain2_ID) = gr-gr-BoundaryID
      !        returns the grain-grain boundary ID, given two grains.
      !        returns 0, if the two grains do not share a boundary
      ! gr_grBoundaryNodes(gr-gr-BoundaryID,nodeIndex) = nodeID
      !        stores the IDs of nodes on individual gr-gr boundaries
      !
      !  faceNodes(face,i) - list of nodeIDs on each face
      !  nFaceNodes(face)  - number of nodes on each face
      !  cornerNodes(1:8)  - nodes at the corner of the cubic domain
      !  domainRange(3,2)  - min and max coordinates of nodes at each direction
      ! assumes TETRAHEDRAL elements !  
      ! ********************************************
      SUBROUTINE identifyBoundaryNodes(G0XYZ,NELX,NX,IJK,elemFaceNeighbors,domainRange)
      
      use meshtools
      implicit none

      integer, intent(in) :: elemFaceNeighbors(NELX,4)
      real(8), intent(in) :: G0XYZ(NX*3)
      integer, intent(in) :: NELX,NX
      integer, intent(in) :: IJK(NELX*4)
      real(8), intent(out) :: domainRange(3,2)

      ! locals
      integer :: iElemID,jElemID,iFace,iLocNode,iGloNode,I,J
      integer :: iGrainID,jGrainID, grainID 
      integer :: sGrainID,lGrainID
      integer :: grgrBndID
      integer :: NFACES,NNODES
      integer :: nGrGrBndNodes
      logical :: found
      logical :: nodeOnDomainBoundary(NX)
      logical :: nodeOnGrGrBoundary(NX)

      integer :: idxBndNode, idxGrain, nodeID, iNode, iDIM, MDIM
      real(8) :: domainLen(3)
      logical :: nodeAlreadyAdded(NX,3)
      real(8) :: pos, tol(3), xNodes(3,4), xNode(3), faceArea, faceAreaVector(3)
      real(8) :: nodeDiagonalPos(8),cornerNodeDiagonalPos(8),diagTransform(8,3)
      real(8) :: elemFaceArea
      integer :: nTotGrainBoundaryNodes,staGrBndNode,endGrBndNode
      integer, allocatable :: grainBoundaryNodes_temp(:,:)
      
      ! assuming TET
      NFACES = 4
      NNODES = 4     
      
      identifiedBoundaryNodes = .false.
      
      nodeOnDomainBoundary(:)=.FALSE.
      nodeOnGrGrBoundary(:)=.FALSE.
      grainIDnode(:)=0
      nDomainBoundaryNodes=0
      domainBoundaryNodes(:)=0
      allocate(grainBoundaryNodes_temp(MAXN_NODES_PER_GRAIN_BND,nGrains))
      allocate(nGrainBoundaryNodes(nGrains))
      allocate(grainBoundaryNodesIndx(nGrains))
      grainBoundaryNodes_temp(:,:) = 0
      nGrainBoundaryNodes(:) = 0    
      grainBoundaryNodesIndx(:) = 0
      gr_grBoundaryIDs(:,:)=0
      nGr_grBoundaries=0
      gr_grBoundaryNodes(:,:)=0
      gr_grBoundaryNodeEffAreas(:,:)=0.d0
      nGrGrBoundaryNodes(:)=0
      gr_grBoundaryAreas(:) = 0.d0
      
      grainNeighborList(:,:) = 0
      nGrainNeighbors(:) = 0
      
      nodeAlreadyAdded(:,:)=.FALSE.
      nFaceNodes(:)=0
      cornerNodes(:)=0

      MDIM=3
      do I=1,MDIM
         ! initialize to a value
         domainRange(I,2) = G0XYZ(I)
         domainRange(I,1) = G0XYZ(I)
         ! find minimum and maximum coordinates
         do iNode=1,NX
            domainRange(I,2) = MAX(domainRange(I,2),  &
                                   G0XYZ(MDIM*(iNode-1)+I))
            domainRange(I,1) = MIN(domainRange(I,1),  &
                                   G0XYZ(MDIM*(iNode-1)+I))
         enddo
      enddo
      domainLen(:)=domainRange(:,2)-domainRange(:,1)
      tol(:)=0.05*domainLen(:)
      
      if(nGrains.EQ.0) return
      
      do iElemID=1,NELX
         ! get element nodal positions
         CALL getElemXNodes(iElemID,IJK,G0XYZ,4,xNodes)
         do iFace=1,NFACES
            jElemID=elemFaceNeighbors(iElemID,iFace)
            if (jElemID.EQ.0) then           ! THIS FACE IS ON THE DOMAIN BOUNDARY
               ! determine the direction it faces
               CALL faceNormal(iFace,xNodes,iDIM)
               ! isolate the 3 nodes on this face:
               do iLocNode = 1,NNODES
                  if (iLocNode.EQ.iFace) cycle  !this node is not on the face
                  ! iFace=1: nodes 2,3,4. iFace=2: nodes 1,3,4 etc.
                  ! see how elemFaceNeighbors is constructed in calcFaceNeighborsTET
                  iGloNode = IJK((iElemID-1)*NNODES+iLocNode)
                  ! search if already added to a domain-face
                  if (.NOT.nodeAlreadyAdded(iGloNode,iDIM)) then
                     pos = xNodes(iDIM,iLocNode)
                     if (pos+tol(iDIM).GT.domainRange(iDIM,2).AND.   & 
                         pos-tol(iDIM).LT.domainRange(iDIM,2)) then
                        ! node is on the lower bound edge in direction iDIM
                        nFaceNodes(iDIM*2)=nFaceNodes(iDIM*2)+1
                        faceNodes(iDIM*2,nFaceNodes(iDIM*2)) = iGloNode
                        nodeAlreadyAdded(iGloNode,iDIM)=.TRUE.
                     elseif (pos+tol(iDIM).GT.domainRange(iDIM,1).AND.   & 
                             pos-tol(iDIM).LT.domainRange(iDIM,1)) then
                        ! node is on the upper bound edge in direction iDIM
                        nFaceNodes(iDIM*2-1)=nFaceNodes(iDIM*2-1)+1     
                        faceNodes(iDIM*2-1,nFaceNodes(iDIM*2-1))=iGloNode
                        nodeAlreadyAdded(iGloNode,iDIM)=.TRUE.
                     endif
                  endif
                  found=.FALSE.
                  do I=1,nDomainBoundaryNodes
                     if(domainBoundaryNodes(I).EQ.iGloNode) then
                        found=.TRUE.
                        exit
                     endif
                  enddo
                  if(.NOT.found) then        ! add to the list of domain boundary nodes
                     nDomainBoundaryNodes=nDomainBoundaryNodes+1
                     domainBoundaryNodes(nDomainBoundaryNodes)   & 
                                                   = iGloNode
                     grainIDnode(iGloNode) = 0
                     nodeOnDomainBoundary(iGloNode) = .TRUE.
                  endif
                  if(nDomainBoundaryNodes.EQ.MAXN_DOM_BND_NODES)then
                     write(*,*) 'INCREASE MAXN_DOM_BND_NODES'
                     STOP
                  endif
                                             ! also add to the list of grain boundary nodes
                  iGrainID = grainIDelem(iElemID)
                  ! search if already added
                  found=.FALSE.
                  do I=1,nGrainBoundaryNodes(iGrainID)
                     if(grainBoundaryNodes_temp(I,iGrainID).EQ.iGloNode) then
                        found=.TRUE.
                        exit
                     endif
                  enddo
                  if(.NOT.found) then
                     if(nGrainBoundaryNodes(iGrainID) >= &
                        MAXN_NODES_PER_GRAIN_BND)then
                        write(*,*) 'MAXN_NODES_PER_GRAIN_BND'
                        STOP
                     endif
                     nGrainBoundaryNodes(iGrainID) =   & 
                                     nGrainBoundaryNodes(iGrainID)+1
                     grainBoundaryNodes_temp(nGrainBoundaryNodes(iGrainID), &
                                        iGrainID)  & 
                                      = iGloNode
                  endif
               enddo
            else                             ! THIS FACE IS WITHIN THE DOMAIN
               iGrainID = grainIDelem(iElemID)
               jGrainID = grainIDelem(jElemID)
               if(iGrainID.EQ.0.OR.jGrainID.EQ.0) cycle
               if(MAX(iGrainID,jGrainID).GT.MAXN_GRAINS)then
                     write(*,*) 'INCREASE MAXN_GRAINS'
                     STOP
               endif
               if (iGrainID.NE.jGrainID) then! THIS FACE IS ON A GRAIN-GRAIN BOUNDARY
               
                  if(jGrainID==0) then
                     write(*,*) 'corrupt structure: an element has a grain ID=0 attached to it'
                     write(*,*) 'element ID',jElemID
                     stop
                  endif
               
                  ! add j-Grain to the list of neighbors of the i-Grain
                  ! search if already added
                  found=.false.
                  do I=1,nGrainNeighbors(iGrainID)
                     if (grainNeighborList(I,iGrainID)==jGrainID) then
                        found = .true.
                        exit
                     endif
                  enddo
                  if (.not.found) then ! add jGrain to the list of neighbors of iGrain
                     nGrainNeighbors(iGrainID) = nGrainNeighbors(iGrainID) + 1
                     if (nGrainNeighbors(iGrainID) >= MAX_GRAIN_NEIGHBORS) then
                        write(*,*) 'increase MAX_GRAIN_NEIGHBORS'
                        STOP
                     endif
                     grainNeighborList(nGrainNeighbors(iGrainID),iGrainID) = jGrainID
                  endif
                  
                  ! calculate the Voronoi/effective area of the nodes on this grain-grain boundary face element
                  CALL faceNormalVector(iFace, xNodes, faceAreaVector)
                  faceArea=dsqrt(dot_product(faceAreaVector,faceAreaVector)) / 2.d0
                  
                  ! now add the nodes on the face to the list of grain-grain boundary nodes
                  ! isolate the 3 nodes on this face:
                  do iLocNode = 1,NNODES
                     if (iLocNode.EQ.iFace) cycle  !this node is not on the face
                      
                                             ! add to the list of grain boundary nodes
                     iGloNode = IJK((iElemID-1)*NNODES+iLocNode)
                     
                     ! search if already added
                     found=.FALSE.
                     do I=1,nGrainBoundaryNodes(iGrainID)
                        if(grainBoundaryNodes_temp(I,iGrainID)==iGloNode) then
                           found=.TRUE.
                           exit
                        endif
                     enddo
                     if(.NOT.found) then
                        nGrainBoundaryNodes(iGrainID) =   & 
                                        nGrainBoundaryNodes(iGrainID)+1
                        grainBoundaryNodes_temp(nGrainBoundaryNodes(iGrainID), &
                                           iGrainID) &
                                         = iGloNode
                     endif

                                             ! add it to the list of gr-gr boundary nodes
                     sGrainID=MIN(iGrainID,jGrainID)
                     lGrainID=MAX(iGrainID,jGrainID)
                     if(gr_grBoundaryIDs(sGrainID,lGrainID).EQ.0) then
                        ! create new ID for gr-gr boundary
                        nGr_grBoundaries=nGr_grBoundaries+1
                        gr_grBoundaryIDs(sGrainID,lGrainID) =    & 
                                                   nGr_grBoundaries
                        gr_grBoundaryIDs(lGrainID,sGrainID) =    & 
                                                   nGr_grBoundaries
                     endif
                     grgrBndID = gr_grBoundaryIDs(sGrainID,lGrainID)
                     if(grgrBndID.GT.MAXN_GR_GR_BND)then
                        write(*,*) 'INCREASE MAXN_GR_GR_BND'
                        STOP
                     endif
                          
                     ! search if this node is added to the gr-gr boundary list already
                     found=.FALSE.
                     do I=1,nGrGrBoundaryNodes(grgrBndID)
                        if(gr_grBoundaryNodes(grgrBndID,I)   & 
                                             .EQ.iGloNode) then
                           found=.TRUE.
                           ! add 1/3 of the face area to the effective area of the node within this gr-gr bnd
                           gr_grBoundaryNodeEffAreas(grgrBndID,I)   & 
                              = gr_grBoundaryNodeEffAreas(grgrBndID,I) &
                              + faceArea / 3.d0
                           exit
                        endif
                     enddo
                     if(.NOT.found) then
                        nGrGrBoundaryNodes(grgrBndID) =   & 
                                 nGrGrBoundaryNodes(grgrBndID) + 1
                        if(nGrGrBoundaryNodes(grgrBndID).GT.   & 
                                 MAXN_NODES_PER_GRGR_BND)then
                        write(*,*) 'MAXN_NODES_PER_GRGR_BND'
                        STOP
                        endif
                        gr_grBoundaryNodes(grgrBndID,   & 
                                        nGrGrBoundaryNodes(grgrBndID))   & 
                                      = iGloNode
                        ! add 1/3 of the face area to the effective area of the node within this gr-gr bnd
                        gr_grBoundaryNodeEffAreas(grgrBndID,   & 
                                     nGrGrBoundaryNodes(grgrBndID))   & 
                          = gr_grBoundaryNodeEffAreas(grgrBndID,   & 
                                   nGrGrBoundaryNodes(grgrBndID)) &
                          + faceArea / 3.d0
                        if(.NOT.nodeOnDomainBoundary(iGloNode)) then
                           grainIDnode(iGloNode) = -1*grgrBndID
                           nodeOnGrGrBoundary(iGloNode) = .TRUE.
                        endif
                     endif
                  enddo   

                  ! add element face area to grain-grain boundary area       
                  ! at this point, an ID must have been assigned to this gr-gr Boundary: gr_grBoundaryIDs(iGrainID,jGrainID)
                  ! to avoid double-counting of areas, only the grain with lower ID calculates the area:
                  if (iGrainID < jGrainID) then
                     gr_grBoundaryAreas(gr_grBoundaryIDs(iGrainID,jGrainID)) = &
                        gr_grBoundaryAreas(gr_grBoundaryIDs(iGrainID,jGrainID)) + faceArea
                  endif
                  
               else                          ! THIS FACE IS WITHIN A GRAIN
                  ! mark all face nodes as 'inside a grain'
                  ! if one of the nodes of this face happen to lie on a grain-grain boundary
                  ! or on the domain boundary, and has not yet been detected as such (nodeOnDomainBoundary or nodeOnGrGrBoundary)
                  ! it will be temporarily marked here as 'inside a grain' ( grainIDnode(iGloNode) = iGrainID  )
                  ! but this node will be correctly marked as 'on grain boundary' or 'on domain boundary' when it is 
                  ! traversed from the face that is actually on the boundary. (see the above two if-conditions)
                  do iLocNode = 1,NNODES
                     if (iLocNode.EQ.iFace) cycle  !this node is not on the face
                      
                                             ! add to the list of grain boundary nodes
                     iGloNode = IJK((iElemID-1)*NNODES+iLocNode)
                     if(.NOT.nodeOnDomainBoundary(iGloNode).AND.   & 
                        .NOT.nodeOnGrGrBoundary(iGloNode))   & 
                        grainIDnode(iGloNode) = iGrainID 
                  enddo
               endif
            endif
         enddo
      enddo
      
      ! compress grainBoundaryNode list into a compact array
      nTotGrainBoundaryNodes = 0
      do iGrainID=1,nGrains
         nTotGrainBoundaryNodes = nTotGrainBoundaryNodes + nGrainBoundaryNodes(iGrainID)
         grainBoundaryNodesIndx(iGrainID) = nTotGrainBoundaryNodes + 1
      enddo
      allocate(grainBoundaryNodes(nTotGrainBoundaryNodes))
      do iGrainID=1,nGrains
         staGrBndNode = 1
         if (iGrainID /= 1) staGrBndNode = grainBoundaryNodesIndx(iGrainID-1)
         endGrBndNode = grainBoundaryNodesIndx(iGrainID)-1
         grainBoundaryNodes(staGrBndNode:endGrBndNode) &
                          = grainBoundaryNodes_temp(1:nGrainBoundaryNodes(iGrainID),iGrainID)
      enddo
      deallocate(grainBoundaryNodes_temp)
      
      ! Identify corner nodes
      do I=1,8
         do J=1,3
            if (BTEST(I-1,J-1)) then
               diagTransform(I,J) = -1
            else
               diagTransform(I,J) = +1
            endif
         enddo
      enddo

      cornerNodeDiagonalPos(:) = MATMUL(diagTransform,G0XYZ(1:3))
      cornerNodes(:)=1
      do iNode=1,NX
         xNode(:) = G0XYZ(MDIM*(iNode-1)+1:MDIM*(iNode-1)+3)
         nodeDiagonalPos = MATMUL(diagTransform,xNode)
         do I=1,8
            if(nodeDiagonalPos(I).LT.cornerNodeDiagonalPos(I)) then
               cornerNodes(I)=iNode
               cornerNodeDiagonalPos(I)=nodeDiagonalPos(I)
            endif
         enddo
      enddo
      
      ! identify grains on the domain boundary      
      allocate(isGrainOnBoundary(nGrains))
      isGrainOnBoundary(:) = .false.
      do idxBndNode=1,nDomainBoundaryNodes
         nodeID = domainBoundaryNodes(idxBndNode)
         do idxGrain=1,nGrainListsForNodes(nodeID)
            grainID = grainListsForNodes(nodeID,idxGrain)
            isGrainOnBoundary(grainID) = .true.
         enddo
      enddo
      nGrainsOnBoundary = 0
      do grainID=1,nGrains
         if (isGrainOnBoundary(grainID)) nGrainsOnBoundary = nGrainsOnBoundary + 1
      enddo
      
      ! success
      identifiedBoundaryNodes = .true.

      END SUBROUTINE
      
      SUBROUTINE getGrainIdxForNode(iGloNode,iGrain,grainIdx,error)
      implicit none
      integer, intent(in) :: iGloNode, iGrain
      integer, intent(out):: grainIdx
      integer, intent(out):: error
      
      logical :: found
      found = .FALSE.
      
      grainIdx = 0
      error = 1
      if (.not.identifiedGrainNodes) return
      
      do grainIdx=1,MAX_GR_PER_NODE
         if(grainListsForNodes(iGloNode,grainIdx)==iGrain) then
            found = .TRUE.
            error = 0
            exit
         endif
      enddo

      END SUBROUTINE
      
      SUBROUTINE calcGrainCentriods(elemCentroids,elemVolumes,NELX)
      implicit none
      
      real(8), intent(in) :: elemCentroids(3,NELX)
      real(8), intent(in) :: elemVolumes(NELX)
      integer, intent(in) :: NELX
      ! locals
      integer :: iGrain,staIndex,endIndex,elemIdx,elemID
      real(8) :: elemVol
      if (.not.grainsImported) return
      
      if (.not.allocated(grainCentroids)) allocate(grainCentroids(3,nGrains))
      if (.not.allocated(grainVolumes)) allocate(grainVolumes(nGrains))
      
      grainCentroids = 0.d0
      grainVolumes = 0.d0
      
      do iGrain=1,nGrains
         staIndex = 1
         if (iGrain.NE.1) staIndex = grainElementIndx(iGrain-1)
         endIndex = grainElementIndx(iGrain)-1
         do elemIdx=staIndex,endIndex
            elemID = grainElements(elemIdx)
            elemVol = elemVolumes(elemID)
            grainCentroids(:,iGrain) = grainCentroids(:,iGrain) + &
                                       elemVol*elemCentroids(:,elemID)
            grainVolumes(iGrain) = grainVolumes(iGrain) + elemVol
         enddo
         
         grainCentroids(:,iGrain) = grainCentroids(:,iGrain) / grainVolumes(iGrain)
      enddo
      
      END SUBROUTINE
      

      SUBROUTINE calcAvgGrainMisorientations(grainMisAngle,grainMisInd,grainRot,cAxis)
      implicit none
      real(8), intent(out) ::  grainMisAngle(:),grainMisInd(:)
      real(8), intent(in)  :: grainRot(:,:,:)
      real(8), intent(in)  :: cAxis(3)
      ! locals
      real(8), parameter :: PI = 3.14159265359d0
      integer :: iGrainID, jGrainID, grainIdx
      real(8) :: grainSurface
      real(8) :: thisBoundaryArea,thisMisAngle,thisMisInd
      
      grainMisAngle(:) = 0.d0
      grainMisInd(:) = 0.d0
      
      if (.not.grainsImported .or. .not.textureImported .or. .not.identifiedBoundaryNodes) then
         grainMisInd(:) = -1.d0
         return
      endif
      if (size(grainMisInd) < nGrains .or. &
          size(grainRot,1) /= 3 .or. &
          size(grainRot,2) /= 3 .or. &
          size(grainRot,3) < nGrains) then     ! allocated to correct size?
         grainMisInd(:) = -1.d0
         return
      endif
      
      do iGrainID=1,nGrains
         
         grainSurface = 0.d0
         do grainIdx = 1,nGrainNeighbors(iGrainID)
            jGrainID = grainNeighborList(grainIdx,iGrainID)

            thisBoundaryArea = gr_grBoundaryAreas(gr_grBoundaryIDs(iGrainID,jGrainID))
            grainSurface = grainSurface + thisBoundaryArea
            thisMisInd = DOT_PRODUCT(MATMUL(grainRot(:,:,iGrainID),cAxis(:)), &
                                     MATMUL(grainRot(:,:,jGrainID),cAxis(:)))
            thisMisInd = dabs(thisMisInd)
            ! avoid PFE due to dot product slightly exceeding 1 for identical rotation matrices:
            if (thisMisInd > 1.d0 .and. thisMisInd < 1.000001d0) thisMisInd = 1.d0 
            thisMisAngle = dacos(thisMisInd) * (180.d0/PI)
            grainMisInd(iGrainID) = grainMisInd(iGrainID) + thisMisInd*thisBoundaryArea
            grainMisAngle(iGrainID) = grainMisAngle(iGrainID) + thisMisAngle*thisBoundaryArea
         enddo
         grainMisInd(iGrainID) = grainMisInd(iGrainID) / grainSurface
         grainMisInd(iGrainID) = 1.d0 - grainMisInd(iGrainID)
         grainMisAngle(iGrainID) = grainMisAngle(iGrainID) / grainSurface
                     
      enddo
      
      END SUBROUTINE
      

      SUBROUTINE calcGrGrMisorientations(grgrMisAngle,grgrMisInd,grainRot,cAxis)
      implicit none
      real(8), intent(out) ::  grgrMisAngle(:),grgrMisInd(:)
      real(8), intent(in)  :: grainRot(:,:,:)
      real(8), intent(in)  :: cAxis(3)
      ! locals
      real(8), parameter :: PI = 3.14159265359d0
      integer :: grgrBoundaryID, grainIdx, iGrainID, jGrainID
      real(8) :: thisBoundaryArea,thisMisAngle,thisMisInd
      
      grgrMisAngle(:) = 0.d0
      grgrMisInd(:) = 0.d0
      
      if (.not.grainsImported .or. .not.textureImported .or. .not.identifiedBoundaryNodes) then
         grgrMisInd(:) = -1.d0
         return
      endif
      if (size(grgrMisAngle) < nGr_grBoundaries .or. &
          size(grgrMisInd) < nGr_grBoundaries .or. &
          size(grainRot,1) /= 3 .or. &
          size(grainRot,2) /= 3 .or. &
          size(grainRot,3) < nGrains) then     ! allocated to correct size?
         grgrMisInd(:) = -1.d0
         return
      endif
      
      do iGrainID=1,nGrains
         
         do grainIdx = 1,nGrainNeighbors(iGrainID)
            jGrainID = grainNeighborList(grainIdx,iGrainID)
            grgrBoundaryID = gr_grBoundaryIDs(iGrainID,jGrainID)
            thisBoundaryArea = gr_grBoundaryAreas(grgrBoundaryID)
            thisMisInd = DOT_PRODUCT(MATMUL(grainRot(:,:,iGrainID),cAxis(:)), &
                                     MATMUL(grainRot(:,:,jGrainID),cAxis(:)))
            thisMisInd = dabs(thisMisInd)
            ! avoid PFE due to dot product slightly exceeding 1 for identical rotation matrices:
            if (thisMisInd > 1.d0 .and. thisMisInd < 1.000001d0) thisMisInd = 1.d0 
            thisMisAngle = dacos(thisMisInd) * (180.d0/PI)
            grgrMisInd(grgrBoundaryID) = 1.d0 - thisMisInd
            grgrMisAngle(grgrBoundaryID) = thisMisAngle
         enddo
                     
      enddo
      
      END SUBROUTINE
      
      END MODULE grains
