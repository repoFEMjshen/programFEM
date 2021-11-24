      PROGRAM calcPatches
      
      use meshtools
      !use options
      use grains
      use nodalRecovery, only : createPatchMatrices
      
      implicit none
      
      include 'PARDIS.H'
      
      integer :: NX, NELX
      integer, parameter :: NODE = 4
      real(8) :: G0XYZ(MAXCRD)
      integer :: IJK(MNELX)
      integer :: IELEMNO(MAXNODELEM*MAXNODE),NINDX(MAXNODE)
      integer, parameter :: MAX_ALLOWED_PATCH_ELEMENTS = 2000
      integer, parameter :: MAX_ALLOWED_PATCH_NODES = 1000
      character(len=30) :: meshName
      logical, parameter :: DEBUG_patches = .false.
      logical :: mustIncludeConnectedElements
      
      ! element-element connectivity
      integer, allocatable:: elemFaceNeighbors(:,:)
      ! non-local patches
      integer, allocatable:: listPatchElementsNonlocal(:,:)
      integer, allocatable:: listPatchNodesNonlocal(:,:)
      integer, allocatable:: listPatchGlobalNodesNonLocal(:,:)
      integer, allocatable:: nPatchElementsNonLocal(:)
      integer, allocatable:: nPatchNodesNonLocal(:)
      real(8), allocatable:: distPatchElementsNonLocal(:,:)
      ! local patches
      integer, allocatable:: listPatchElementsLocal(:,:)
      integer, allocatable:: listPatchNodesLocal(:,:)
      integer, allocatable:: listPatchGlobalNodesLocal(:,:)
      integer, allocatable:: nPatchElementsLocal(:)
      integer, allocatable:: nPatchNodesLocal(:)
      integer :: nElemsInThisPatch, minElemsInThisPatch
      ! patchID-nodeID mappings
!      integer, allocatable :: nPatchesInGrain(:)
!      integer, allocatable :: patchGrainID(:)
      integer, allocatable :: patchIDFromNodeID(:)
      integer, allocatable :: nodeIDFromPatchID(:)
      ! patch calculation
      real(8), allocatable :: P(:)
      real(8), allocatable :: A(:,:)
      real(8) :: A_norm
      real(8), allocatable :: b(:,:) ! (4,nQuantDOF)
      integer :: p_deg_patch,nPolyTerms,nQuantDOF_patch
      real(8) :: xPatchGauss(3,MAX_ALLOWED_PATCH_ELEMENTS)
      real(8) :: xPatchRange(3,2),xPatchCenter(3)
      real(8), allocatable :: gaussValues_dummy(:,:)
      
      integer :: nMaxPatchesPerGrain, MAX_N_PATCHES
      integer :: nMaxPatchElementsLocal_actual, nMaxPatchElementsNonLocal_actual
      integer :: nMaxPatchNodesLocal_actual, nMaxPatchNodesNonLocal_actual

      real(8) :: xNode(3)
      real(8) :: rcond, rcondMin
      integer :: powCond
      real(8), allocatable :: work(:)
      integer, allocatable :: iwork(:)
      
      integer, allocatable :: patchIDs(:,:)
      integer :: nPatchIDs, patchID

      real(8) :: maxPatchRadius
      real(8), allocatable :: elemCenters(:,:)
      real(8), allocatable :: elemVol(:)
      real(8) :: xNodes(3,4), elemJacMat(3,3), elemJac
      real(8) :: distVec(3), dist, weight, normalizationConstant
      
      logical :: success, successReadNeighbors, successReadNodeConnect
      
      ! variables for domain geometry
      real(8) :: domainRange(3,2), domainLen(3)
      
      integer :: nodeIdx,staNodeIdx,endNodeIdx
      integer :: elemNode,staIndex,endIndex,elemIdx,thisIndex
      integer :: i,iElem, jElem, nElem, iNode, nodeID, iGrain, elemID, n
      integer :: patchNodeIdx
      logical :: found
      
      integer*4 :: lapack_INFO
      
      ! read in the mesh and connectivity
      CALL READGM(G0XYZ,IJK,NX,NELX,meshName)      

      write(*,*) 'Mesh imported:'
      write(*,*) '# of elements: ', NELX
      write(*,*) '# of nodes: ', NX
      
      
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
         write(*,*) 'saved connectivity data file: elementNeighbors.dat'
      endif

      ! IDENTIFY THE NODES IN INDIVIDUAL GRAINS
      CALL identifyGrainNodes(NX,NELX,IJK)

      ! IDENTIFY THE NODES AT 
      ! 1) DOMAIN BOUNDARY, DOMAIN CORNERS
      ! 2) GRAIN BOUNDARIES
      ! 3) GRAIN-GRAIN BOUNDARIES
      CALL identifyBoundaryNodes(G0XYZ,NELX,NX,IJK,elemFaceNeighbors,domainRange)

      ! domainRange0: initial domain range
      ! domainRange: may evolve in time
      domainLen(:) = dabs(domainRange(:,2) - domainRange(:,1))
      
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
      allocate(elemCenters(3,NELX))
      do iElem=1,NELX
         CALL getElemXNodes(iElem,IJK,G0XYZ,NODE,xNodes)
         CALL getElemXCenter(iElem,IJK,G0XYZ,NODE,elemCenters(:,iElem))
         CALL calcJacTET4(xNodes, elemJacMat, elemJac, elemVol(iElem))   
      enddo

      write(*,*) 'enter max radius of each patch:'
      read(*,*) maxPatchRadius

      nMaxPatchesPerGrain = 0
      nPatchIDs = 0
      do iGrain=1,nGrains
         if (nGrainNodes(iGrain) > nMaxPatchesPerGrain) &
            nMaxPatchesPerGrain = nGrainNodes(iGrain)
         nPatchIDs = nPatchIDs + nGrainNodes(iGrain)
      enddo
      allocate(patchIDs(nMaxPatchesPerGrain,nGrains))
      patchIDs = 0
      
      write(*,*) 'number of patches: ',nPatchIDs
      
      allocate(distPatchElementsNonLocal(MAX_ALLOWED_PATCH_ELEMENTS,nPatchIDs))
      allocate(listPatchElementsNonLocal(MAX_ALLOWED_PATCH_ELEMENTS,nPatchIDs))
      allocate(listPatchNodesNonLocal(MAX_ALLOWED_PATCH_NODES,nPatchIDs))
      allocate(listPatchGlobalNodesNonLocal(MAX_ALLOWED_PATCH_NODES,nPatchIDs))
      allocate(nPatchElementsNonLocal(nPatchIDs))
      allocate(nPatchNodesNonLocal(nPatchIDs))
      
      nPatchElementsNonLocal = 0
      listPatchElementsNonlocal = 0
      nPatchNodesNonLocal = 0
      listPatchNodesNonLocal = 0
      listPatchGlobalNodesNonLocal = 0
      distPatchElementsNonLocal = 0.d0
      
      allocate(listPatchElementsLocal(MAX_ALLOWED_PATCH_ELEMENTS,nPatchIDs))
      allocate(listPatchNodesLocal(MAX_ALLOWED_PATCH_NODES,nPatchIDs))
      allocate(listPatchGlobalNodesLocal(MAX_ALLOWED_PATCH_NODES,nPatchIDs))
      allocate(nPatchElementsLocal(nPatchIDs))
      allocate(nPatchNodesLocal(nPatchIDs))
      
      allocate(nodeIDFromPatchID(nPatchIDs))
      allocate(patchIDFromNodeID(NX))
      
      nPatchElementsLocal = 0
      listPatchElementsLocal = 0
      nPatchNodesLocal = 0
      listPatchNodesLocal = 0
      listPatchGlobalNodesLocal = 0
      
      write(*,*) 'require connected elements in the non-local patch?'
      read(*,*) mustIncludeConnectedElements
      
      write(*,*) 'max allowed digits of precision loss? 1: strict -- 10: loose. 0: accept anything'
      read(*,*) powCond
      rcondMin = (1.d-1)**powCond
      
      write(*,*) '.. with respect to what degree of polynomials?'
      read(*,*) p_deg_patch


      if (DEBUG_patches) open(768,file='debugPatches.txt')
      
      write(*,*) 'finding elements for each patch...'
      if (DEBUG_patches) write(768,*) 'finding elements for each patch...'
      patchID = 0
      do iGrain=1,nGrains
         
         staNodeIdx = 1
         if (iGrain.NE.1) staNodeIdx = grainNodesIndices(iGrain-1)
         endNodeIdx = grainNodesIndices(iGrain)-1
         
         do nodeIdx=staNodeIdx,endNodeIdx
         
            patchID = patchID + 1
            patchIDs(nodeIdx-staNodeIdx+1,iGrain) = patchID
            
            iNode=grainNodes(nodeIdx)
            xNode(1:3) = G0XYZ(3*(iNode-1)+1:3*(iNode-1)+3)
            
            ! add to the patchID --> global node ID mapping
            nodeIDFromPatchID(patchID) = iNode
            
            if (DEBUG_patches) write(768,*) 'patch ID',patchID
            if (DEBUG_patches) write(768,*) 'node',iNode,xNode(:),'grain:',iGrain
            
            ! ******* CREATE THE LOCAL PATCH AROUND THE NODE iNode
            
            ! add the central node to the beginning of the list of patch nodes
            nPatchNodesLocal(patchID) = 1
            listPatchGlobalNodesLocal(1,patchID) = iNode
            listPatchNodesLocal(1,patchID) = patchID
            
            ! traverse the elements around the central node
            staIndex=1
            if(iNode/=1) staIndex = NINDX(iNode-1)
            endIndex = NINDX(iNode)-1
               
            do elemIdx=staIndex,endIndex 
               elemID=IELEMNO(elemIdx)
            
               if (grainIDelem(elemID) /= iGrain) cycle   ! only elements in the same grain
               
               nPatchElementsLocal(patchID) = nPatchElementsLocal(patchID) + 1  ! add to the list of patch-elements for this node
               thisIndex = nPatchElementsLocal(patchID)
               listPatchElementsLocal(thisIndex,patchID) = elemID
               
               ! add the nodes of the element to the list of nodes on the patch
               do elemNode = 1,NODE
                  nodeID = IJK((elemID-1)*4+elemNode)
                  thisIndex = nPatchNodesLocal(patchID) + 1
                  if (any(listPatchGlobalNodesLocal(1:thisIndex-1,patchID) == nodeID)) cycle
                  listPatchGlobalNodesLocal(thisIndex,patchID) = nodeID ! add only if not on the list already
                  nPatchNodesLocal(patchID) = thisIndex
               enddo
         
            enddo
            
            if (DEBUG_patches) then
               write(768,*) 'local patch got #elems:',nPatchElementsLocal(patchID)
               write(768,*) listPatchElementsLocal(1:nPatchElementsLocal(patchID),patchID)
            endif
            if (DEBUG_patches) write(768,*) 'local patch got #nodes:',nPatchNodesLocal(patchID)
            
            
            ! ******* CREATE THE NON-LOCAL PATCH AROUND THE NODE iNode
            ! add the central node to the beginning of the list of patch nodes
            nPatchNodesNonLocal(patchID) = 1
            listPatchGlobalNodesNonLocal(1,patchID) = iNode
            listPatchNodesNonLocal(1,patchID) = patchID
            ! initialize the nonlocal patch with the local patch.
            ! this makes sure that elements connected to the node are in the non-local patch
            ! irrespective of their spatial distance to the node
            if (mustIncludeConnectedElements) then 
               nPatchElementsNonLocal(patchID) = nPatchElementsLocal(patchID)
               listPatchElementsNonLocal(:,patchID) = listPatchElementsLocal(:,patchID)
               do elemIdx = 1,nPatchElementsLocal(patchID)  ! also calculate and assign distances, as it is required for sorting
                  iElem = listPatchElementsNonLocal(elemIdx,patchID)
                  distVec(:) = elemCenters(:,iElem) - xNode(:)
                  dist = dsqrt(DOT_PRODUCT(distVec,distVec))
                  distPatchElementsNonLocal(elemIdx,patchID) = dist
               enddo
            endif
            ! the other nodes on the patch will be added after the actual patch size is determined
            ! now search through all the elements in the grain to find non-immediate neighbors of the node
            staIndex = 1
            if (iGrain.NE.1) staIndex = grainElementIndx(iGrain-1)
            endIndex = grainElementIndx(iGrain)-1
            do elemIdx=staIndex,endIndex
               jElem = grainElements(elemIdx)
               
               if (mustIncludeConnectedElements) then ! local patch elements were already added. skip, if this is one of them.
                  if (any(listPatchElementsLocal(1:nPatchElementsLocal(patchID),patchID) == jElem)) cycle
               endif
               
               distVec(:) = elemCenters(:,jElem) - xNode(:)
               dist = dsqrt(DOT_PRODUCT(distVec,distVec))
               
               if (dist < maxPatchRadius) then
               
                  nPatchElementsNonLocal(patchID) = nPatchElementsNonLocal(patchID) + 1
                  thisIndex = nPatchElementsNonLocal(patchID)
                  if (thisIndex > MAX_ALLOWED_PATCH_ELEMENTS) then
                     write(*,*) 'increase MAX_ALLOWED_PATCH_ELEMENTS'
                     stop
                  endif
                  listPatchElementsNonlocal(thisIndex,patchID) = jElem
                  ! xPatchElementsNonLocal(thisIndex,patchID) = elemCenters(:,jElem) ----> 
                     ! ---> potential optimization. (so that you don't re-calculate each time while creating patches)
                  distPatchElementsNonLocal(thisIndex,patchID) = dist
               endif
            enddo
         
            write(*,'(I0,A,I0,A,I0,A,I0)') thisIndex,' elements in patchID ',patchID,' - node ', &
                                         iNode,' grain ',iGrain      
            if (DEBUG_patches) write(768,*) 'non local #elements:',nPatchElementsNonLocal(patchID)                                         
            if (DEBUG_patches) write(768,*) listPatchElementsNonlocal(1:thisIndex,patchID)
            
         enddo
         
      enddo
      
      if (DEBUG_patches) flush(768)
      
      if(patchID /= nPatchIDs) then
         write(*,*) 'assertion errror'
         stop
      endif
      
      write(*,*) 'sorting elements within each non-local patch...'
      if (DEBUG_patches) write(768,*) 'sorting elements within each non-local patch...'
      
      do patchID=1,nPatchIDs
   
         ! sort patch elements w.r.t. distance
         if (DEBUG_patches) write(768,*) 'patchID',patchID
         if (DEBUG_patches) write(768,*) 'BEFORE dist:', distPatchElementsNonLocal(1:nPatchElementsNonLocal(patchID),patchID)
         if (DEBUG_patches) write(768,*) 'BEFORE list:', listPatchElementsNonLocal(1:nPatchElementsNonLocal(patchID),patchID)
         
         if (mustIncludeConnectedElements) then ! do not include locally-connected elements in the sorting to ensure they remain in the patch
            CALL Sort(distPatchElementsNonLocal(nPatchElementsLocal(patchID)+1:nPatchElementsNonLocal(patchID),patchID), &
                      listPatchElementsNonlocal(nPatchElementsLocal(patchID)+1:nPatchElementsNonLocal(patchID),patchID), &
                                                           nPatchElementsNonLocal(patchID)-nPatchElementsLocal(patchID))            
         else
            CALL Sort(distPatchElementsNonLocal(:,patchID),listPatchElementsNonlocal(:,patchID), &
                   nPatchElementsNonLocal(patchID))
         endif
                   
         if (DEBUG_patches) write(768,*) 'AFTER  dist:', distPatchElementsNonLocal(1:nPatchElementsNonLocal(patchID),patchID)
         if (DEBUG_patches) write(768,*) 'AFTER  list:', listPatchElementsNonLocal(1:nPatchElementsNonLocal(patchID),patchID)
      
      enddo
      
      if (DEBUG_patches) flush(768)
      
      write(*,*) 'finding smallest patch sizes with well-conditioned matrices...'
      if (DEBUG_patches) write(768,*) 'finding smallest patch sizes with well-conditioned matrices...'
      nQuantDOF_patch = 1
      !p_deg_patch = 2   ! quadratic patch
                        ! 10 coefficients in the 3D Quad polynomial
      nPolyTerms = (p_deg_patch+3)*(p_deg_patch+2)*(p_deg_patch+1)/6
                        
      allocate(P(nPolyTerms))
      allocate(A(nPolyTerms,nPolyTerms))
      allocate(b(nPolyTerms,nQuantDOF_patch)) 
      allocate(gaussValues_dummy(NELX,nQuantDOF_patch))
      gaussValues_dummy = 0.d0
      ! find needed work-space          
       allocate(work(10*3*nPolyTerms))
       allocate(iwork(10*nPolyTerms))
      
      do patchID=1,nPatchIDs
         
         if (DEBUG_patches) write(768,*) 'patchID',patchID
         
         ! save element gauss pt centers in this patch
         do elemIdx = 1,nPatchElementsNonLocal(patchID)
            xPatchGauss(:,elemIdx) = elemCenters(:,listPatchElementsNonLocal(elemIdx,patchID))
         enddo
         
         if (DEBUG_patches) &
            CALL printMatrix(xPatchGauss(:,1:nPatchElementsNonLocal(patchID)), &
                          3,nPatchElementsNonLocal(patchID),   &
                          'xPatchGauss',768)
                          
         ! NORMALIZE
         xPatchCenter = 0.d0
         do elemIdx = 1,nPatchElementsNonLocal(patchID)
            xPatchCenter(:)=xPatchCenter(:)+xPatchGauss(:,elemIdx)
         enddo
         xPatchCenter(:) = xPatchCenter(:) / nPatchElementsNonLocal(patchID)
         
         do elemIdx = 1,nPatchElementsNonLocal(patchID)
            xPatchGauss(:,elemIdx)=xPatchGauss(:,elemIdx)-xPatchCenter(:)
         enddo
         !calculate the spatial range of the patch
         !small or skewed range results in ill-conditioned matrix
         xPatchRange(:,1) = xPatchGauss(:,1) ! stores the min
         xPatchRange(:,2) = xPatchGauss(:,1) ! stores the max
         do elemIdx = 1,nPatchElementsNonLocal(patchID)
            xPatchRange(1,1)= MIN(xPatchGauss(1,elemIdx),xPatchRange(1,1)) ! X
            xPatchRange(2,1)= MIN(xPatchGauss(2,elemIdx),xPatchRange(2,1)) ! Y
            xPatchRange(3,1)= MIN(xPatchGauss(3,elemIdx),xPatchRange(3,1)) ! Z
            xPatchRange(1,2)=MAX(xPatchGauss(1,elemIdx),xPatchRange(1,2))  ! X
            xPatchRange(2,2)=MAX(xPatchGauss(2,elemIdx),xPatchRange(2,2))  ! Y
            xPatchRange(3,2)=MAX(xPatchGauss(3,elemIdx),xPatchRange(3,2))  ! Z
         enddo

         ! NORMALIZE all gauss point positions
         do elemIdx =1,nPatchElementsNonLocal(patchID)
            xPatchGauss(:,elemIdx) = &
             -1 + 2*(xPatchGauss(:,elemIdx)-xPatchRange(:,1)) &
                     / (xPatchRange(:,2)-xPatchRange(:,1))
         enddo
         if (DEBUG_patches) &
            CALL printMatrix(xPatchGauss(:,1:nPatchElementsNonLocal(patchID)), &
                          3,nPatchElementsNonLocal(patchID),   &
                          'xPatchGauss NORMALIZED',768)
         
         ! find the optimal-minimum patch size
         found = .false.
         minElemsInThisPatch = nPolyTerms
         if (mustIncludeConnectedElements) then
            ! include the  locally-connected elements, and grow the patch as necessary
            minElemsInThisPatch = nPatchElementsLocal(patchID)
         endif
         do nElemsInThisPatch=minElemsInThisPatch, nPatchElementsNonLocal(patchID)
         
            if (DEBUG_patches) write(768,*) 'trying nElemsInThisPatch=',nElemsInThisPatch
         
            ! construct a from the first nElemsInThisPatch .. elements
            
            
            CALL createPatchMatrices(A,b,P,nPolyTerms,nQuantDOF_patch,p_deg_patch, &
                                     listPatchElementsNonLocal(:,patchID), &
                                     nElemsInThisPatch,xPatchGauss, &
                                     distPatchElementsNonLocal(:,patchID),gaussValues_dummy, &
                                     NX,NELX,.false.)            
                                     
            if (DEBUG_patches) CALL printMatrix(A,nPolyTerms,nPolyTerms,   &
                          'matrixA',768)
            
            ! test conditioning
            ! first calculate the norm of A
            CALL calcNorm(A_norm,A,nPolyTerms)
            
            ! first get the Cholesky factorization of A (symmetric positive definite)
            CALL dpotrf('U',nPolyTerms,A,nPolyTerms,lapack_INFO)
            if(lapack_INFO.NE.0 .and. powCond/=0) then
               if (DEBUG_patches) write(768,*) 'dpotrf error', lapack_INFO,' increasing patch size...'
               cycle
            endif
            
            call dpocon('U',nPolyTerms,A,nPolyTerms,A_norm,rcond, &
                        work, iwork,lapack_INFO)
            if(lapack_INFO.NE.0 .and. powCond/=0) then
               if (DEBUG_patches) write(768,*) 'dpocon error', lapack_INFO,' increasing patch size...'
               cycle
            endif
            
           if (DEBUG_patches)  write(768,*) 'rcond',rcond
            
            !call dposv('U',nPolyTerms,nQuantDOF_patch,A,nPolyTerms,b,nPolyTerms, &
            !                                                lapack_INFO)
            
            if (rcond > rcondMin .or. powCond==0) then 
               found = .true.
               exit
            endif
            
         enddo
         
         if (DEBUG_patches) flush(768)
         
         if (found) then
            write(*,*) 'patch size',nPatchElementsNonLocal(patchID),'--->',nElemsInThisPatch
            nPatchElementsNonLocal(patchID) = nElemsInThisPatch
            write(*,*) 'rcond',rcond
            if (DEBUG_patches) write(768,*) 'rcondFOUND !!'
            if (DEBUG_patches) flush(768)
         else
            write(*,*) 'cannot reduce patch size with available elements!'
            write(*,*) nElemsInThisPatch,minElemsInThisPatch, nPatchElementsNonLocal(patchID)
            if (DEBUG_patches) write(768,*) 'cannot reduce patch size with available elements!'
            if (DEBUG_patches) flush(768)
            stop
         endif
         
         ! now that the actual number of elements in the non-local patches are determined,
         ! construct the list of nodes in the non-local patch
         do elemIdx=1,nPatchElementsNonLocal(patchID)
            
            elemID = listPatchElementsNonLocal(elemIdx,patchID)

            ! add the nodes of the element to the list of nodes on the patch
            do elemNode = 1,NODE
               nodeID = IJK((elemID-1)*4+elemNode)
               thisIndex = nPatchNodesNonLocal(patchID) + 1
               if (any(listPatchGlobalNodesNonLocal(1:thisIndex-1,patchID) == nodeID)) cycle
               listPatchGlobalNodesNonLocal(thisIndex,patchID) = nodeID ! add only if not on the list already
               nPatchNodesNonLocal(patchID) = thisIndex
            enddo
         
         enddo

      enddo
      
      ! now create the list of patchNodeIDs for each local patch
      do iGrain=1,nGrains
               
         ! initialize the temporary inverse mapping array
         patchIDFromNodeID = 0

         staNodeIdx = 1
         if (iGrain.NE.1) staNodeIdx = grainNodesIndices(iGrain-1)
         endNodeIdx = grainNodesIndices(iGrain)-1
         
         do nodeIdx=staNodeIdx,endNodeIdx
         
            iNode=grainNodes(nodeIdx)
            patchID = patchIDs(nodeIdx-staNodeIdx+1,iGrain)
                  
            ! add to the temporary array for inverse mapping: global Node ID --> patchID (only valid for this grain)
            patchIDFromNodeID(iNode) = patchID
         enddo
         
         ! will use the inverse mapping patchIDFromNodeID, which is valid for this grain only.
         do nodeIdx=1,nGrainNodes(iGrain)
         
            patchID = patchIDs(nodeIdx,iGrain)
                  
            ! add to the temporary array for inverse mapping: global Node ID --> patchID (only valid for this grain)
            ! local patches
            do patchNodeIdx=1,nPatchNodesLocal(patchID)
               nodeID = listPatchGlobalNodesLocal(patchNodeIdx,patchID)
               listPatchNodesLocal(patchNodeIdx,patchID) = patchIDFromNodeID(nodeID)         
            enddo
            ! non-local patches
            do patchNodeIdx=1,nPatchNodesNonLocal(patchID)
               nodeID = listPatchGlobalNodesNonLocal(patchNodeIdx,patchID)
               listPatchNodesNonLocal(patchNodeIdx,patchID) = patchIDFromNodeID(nodeID)         
            enddo
         enddo
         
      enddo
      
      ! write to file
      nMaxPatchElementsLocal_actual = MAXVAL(nPatchElementsLocal(1:nPatchIDs))
      nMaxPatchElementsNonLocal_actual = MAXVAL(nPatchElementsNonLocal(1:nPatchIDs))
      nMaxPatchNodesLocal_actual = MAXVAL(nPatchNodesLocal(1:nPatchIDs))
      nMaxPatchNodesNonLocal_actual = MAXVAL(nPatchNodesNonLocal(1:nPatchIDs))
      
      open(unit=101,file='listPatches.dat')
      write(101,*) nPatchIDs
      write(101,*) 'patchIDs'
      write(101,*) NX,NELX,nGrains,NODE
      write(101,*) nMaxPatchesPerGrain
      write(101,*) nMaxPatchElementsLocal_actual
      write(101,*) nMaxPatchNodesLocal_actual
      write(101,*) nMaxPatchElementsNonLocal_actual
      write(101,*) nMaxPatchNodesNonLocal_actual
      do iGrain=1,nGrains
         write(101,*) nGrainNodes(iGrain), (patchIDs(nodeIdx,iGrain),nodeIdx=1,nGrainNodes(iGrain))
         write(101,*) (nodeIDFromPatchID(patchIDs(nodeIdx,iGrain)),nodeIdx=1,nGrainNodes(iGrain))
      enddo
      write(101,*) 'listPatchElementsLocal'
      do patchID=1,nPatchIDs
         write(101,*) nPatchElementsLocal(patchID)
         write(101,*) (listPatchElementsLocal(elemIdx,patchID),elemIdx=1,nPatchElementsLocal(patchID))
         write(101,*) nPatchNodesLocal(patchID)
         write(101,*) (listPatchNodesLocal(nodeIdx,patchID),nodeIdx=1,nPatchNodesLocal(patchID))
         write(101,*) (listPatchGlobalNodesLocal(nodeIdx,patchID),nodeIdx=1,nPatchNodesLocal(patchID))
      enddo
      write(101,*) 'listPatchElementsNonLocal'
      do patchID=1,nPatchIDs
         write(101,*) nPatchElementsNonLocal(patchID)
         write(101,*) (listPatchElementsNonLocal(elemIdx,patchID),elemIdx=1,nPatchElementsNonLocal(patchID))
         !write(101,*) (distPatchElementsNonLocal(elemIdx,patchID),elemIdx=1,nPatchElementsNonLocal(patchID))
         write(101,*) nPatchNodesNonLocal(patchID)
         write(101,*) (listPatchNodesNonLocal(nodeIdx,patchID),nodeIdx=1,nPatchNodesNonLocal(patchID))
         write(101,*) (listPatchGlobalNodesNonLocal(nodeIdx,patchID),nodeIdx=1,nPatchNodesNonLocal(patchID))
      enddo 
      write(101,*) 'isBoundaryNode'
      do patchID=1,nPatchIDs
         write(101,*) (grainIDnode(nodeIDFromPatchID(patchID)) <= 0) ! if <=0 on grain or domain boundary. if > 0 in the bulk of a grain.
      enddo
      
      close(101)
      if (DEBUG_patches) close(768)
      
      write(*,*) 'listPatches.dat -- created'
      
      deallocate(elemFaceNeighbors)
      
      deallocate(listPatchElementsNonlocal,nPatchElementsNonLocal,distPatchElementsNonLocal)
      deallocate(listPatchNodesNonlocal,nPatchNodesNonLocal)
      deallocate(listPatchElementsLocal,nPatchElementsLocal)
      deallocate(listPatchNodesLocal,nPatchNodesLocal)
      deallocate(patchIDs)
      deallocate(elemCenters,elemVol)
      
      deallocate(P,A,b)
      deallocate(gaussValues_dummy)
      deallocate(work,iwork)
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
         write(*,*) 'Increase MAXNODE to ', NX, MAXNODE
         WRITE(*,*)'Increase MAXNODE to ', NX, MAXNODE
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
      

      SUBROUTINE  Sort(x_ref, x, Size)
         IMPLICIT  NONE
         real(8), INTENT(INOUT)                :: x_ref(1:Size)
         INTEGER, INTENT(INOUT)                :: x(1:Size)
         INTEGER, INTENT(IN)                   :: Size
         INTEGER                               :: i
         INTEGER                               :: Location
         real(8) :: tempReal
         integer :: tempInteger

         DO i = 1, Size-1			! except for the last
            Location = MINLOC(x_ref(i:Size), 1)+i-1	! find min from this to last
            ! swap this and the minimum
            tempReal = x_ref(i)
            x_ref(i) = x_ref(Location)
            x_ref(Location) = tempReal
            ! swap this and the minimum
            tempInteger = x(i)
            x(i) = x(Location)
            x(Location) = tempInteger
         END DO
      END SUBROUTINE  Sort
      
      SUBROUTINE calcNorm(A_norm,A,n)
      implicit none
      real(8), intent(out) :: A_norm
      real(8), intent(in) :: A(n,n)
      integer, intent(in) :: n
      !locals
      real(8) :: row_sum
      integer :: i,j
      
      A_norm = 0.d0
      do i=1,n
         row_sum = 0.d0
         do j=1,n
            row_sum = row_sum + dabs(A(i,j))
         enddo
         if (row_sum > A_norm) then
            A_norm = row_sum
         endif
      enddo
      END SUBROUTINE
      
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
      
      SUBROUTINE MATINV3_UEL(A,AI,DET)
      IMPLICIT NONE
      REAL(8), intent(in) :: A(3,3)
      REAL(8), intent(out):: AI(3,3)
      REAL(8), intent(out):: DET

      DET=(A(1,1)*A(2,2)*A(3,3)-A(1,1)*A(2,3)*A(3,2)-A(2,1)*A(1,2) & 
            *A(3,3)+A(2,1)*A(1,3)*A(3,2)+A(3,1)*A(1,2)*A(2,3)-A(3,1) & 
           *A(1,3)*A(2,2))

      AI(1,1) =  ( A(2,2)*A(3,3)-A(2,3)*A(3,2))/DET
      AI(1,2) = -( A(1,2)*A(3,3)-A(1,3)*A(3,2))/DET
      AI(1,3) =  ( A(1,2)*A(2,3)-A(1,3)*A(2,2))/DET
      AI(2,1) = -( A(2,1)*A(3,3)-A(2,3)*A(3,1))/DET
      AI(2,2) =  ( A(1,1)*A(3,3)-A(1,3)*A(3,1))/DET
      AI(2,3) = -( A(1,1)*A(2,3)-A(1,3)*A(2,1))/DET
      AI(3,1) =  ( A(2,1)*A(3,2)-A(2,2)*A(3,1))/DET
      AI(3,2) = -( A(1,1)*A(3,2)-A(1,2)*A(3,1))/DET
      AI(3,3) =  ( A(1,1)*A(2,2)-A(1,2)*A(2,1))/DET
      RETURN
      END SUBROUTINE