!  getElemXNodes
!  getElemXCenter
!  calcFaceAreaTET
!  faceNormal
!  getLargestEdge
!  checkInside
!  NRM3
!  assignElementsToBoxes
!  getDomainRange
!  getDomainLen
!  findContainingElement
!  calcJacTET --- unsupported
!  calcJacTET4
!  calcFrobeniusNorm
!  calcCurlTensor
!  calcGradTensorTET
!  approximateDisplField
!  calcDistanceNormCoeffs
!  calcDistanceMatrix

      module meshTools

      private :: cross

      contains
      

      SUBROUTINE readFaceNeighborsTET(strNeighborFile,elemFaceNeighbors,NELX,success)
      implicit none
      
      character(len=*), intent(in) :: strNeighborFile      
      integer, intent(out) :: elemFaceNeighbors(NELX,4)
      logical, intent(out) :: success
      integer, intent(in)  :: NELX
      !Local
      logical :: fileExists
      integer :: error, iElem, nElem_File

      success = .false.
      elemFaceNeighbors(:,:)=0
      
      inquire(file=trim(strNeighborFile),exist=fileExists)
      if(.NOT.fileExists) return
      
      open(716,file=trim(strNeighborFile))
      read(716,*,iostat=error) nElem_File
      if(error.NE.0 .or. nElem_File.ne.NELX) then
         close(716)
         return
      endif
      do iElem=1,NELX
         read(716,*,iostat=error) elemFaceNeighbors(iElem,1:4)
         if(error.NE.0) then
            close(716)
            return
         endif
      enddo
      
      !success!
      success = .true.
      close(716)
      return
      
      end subroutine

      
      SUBROUTINE writeFaceNeighborsTET(strNeighborFile,elemFaceNeighbors,NELX)
      implicit none
      
      character(len=*), intent(in) :: strNeighborFile      
      integer, intent(in) :: elemFaceNeighbors(NELX,4)
      integer, intent(in)  :: NELX
      !Local
      integer :: iElem
            
      open(716,file=trim(strNeighborFile))
      write(716,*) NELX
      do iElem=1,NELX
         write(716,*) elemFaceNeighbors(iElem,1:4)
      enddo
      close(716)
      return
      
      end subroutine

      SUBROUTINE calcFaceNeighborsTET(elemFaceNeighbors,IJK,IELEMNO,NINDX, &
                                      NELX,NODE,MAXNODELEM,MAXNODE)

      implicit none
      ! arguments
      integer, intent(out):: elemFaceNeighbors(NELX,4)
      integer, intent(in) :: IJK(4*NELX)
      integer, intent(in) :: IELEMNO(MAXNODELEM*MAXNODE),NINDX(MAXNODE)
      integer, intent(in) :: MAXNODELEM,MAXNODE
      integer, intent(in) :: NELX,NODE
      !Local
      INTEGER:: I,J,II,JJ,neighFace,elementFirst,elementLast
      integer:: neighCount, nodeCount, iNode, nodeID, iElem, elemID, nElem_File
      integer:: NFACES
      logical :: neighborNode(NODE),fileExists

      !Initialize
      NFACES=4
      elemFaceNeighbors(:,:)=0
      neighFace=0
      
      DO I=1,NELX
         neighCount=0
         do iNode=1,NODE
            nodeID = IJK((I-1)*NODE+iNode)
            if(nodeID==1) then 
               elementFirst = 1
            else
               elementFirst = NINDX(nodeID-1)
            endif
            elementLast = NINDX(nodeID)-1
            
            DO iElem=elementFirst,elementLast
               J = IELEMNO(iElem)
               if(J==I) cycle
               if (elemFaceNeighbors(I,1)==J) cycle
               if (elemFaceNeighbors(I,2)==J) cycle
               if (elemFaceNeighbors(I,3)==J) cycle
               if (elemFaceNeighbors(I,4)==J) cycle
               
               ! count the number of common nodes
               ! if 3, they're sharing a common face.
               nodeCount=0
               neighborNode(:)=.FALSE.
               DO II=1,NODE
                  DO JJ=1,NODE
                     IF(IJK((I-1)*NODE+II).EQ.IJK((J-1)*NODE+JJ))THEN
                        neighborNode(II)=.TRUE.
                        nodeCount=nodeCount+1
                     ENDIF
                  ENDDO
               ENDDO
               if(nodeCount.EQ.3)then
                  neighCount=neighCount+1
                  ! DETERMINE WHICH FACE of ELEMENT (I) 
                  ! the connection is through
                  ! then do elemFaceNeighbors(I,iFACE)=J
                  if(.NOT.neighborNode(1)) neighFace=1
                  if(.NOT.neighborNode(2)) neighFace=2
                  if(.NOT.neighborNode(3)) neighFace=3
                  if(.NOT.neighborNode(4)) neighFace=4
                  
                  elemFaceNeighbors(I,neighFace)=J
               endif
               
               !if number of neighbors found = number of faces. skip to the next element
               if(neighCount.EQ.NFACES) exit
                  
            ENDDO
            if(neighCount.EQ.NFACES) exit
         enddo
      ENDDO
      END SUBROUTINE
            
      SUBROUTINE calcFaceNeighborsTET_FAST(elemFaceNeighbors,elemsInBox,nElemsInBox,nBoxes, &
                                           domainRange,IJK,G0XYZ,NELX,NX,NODE)
      implicit none
      ! arguments
      integer, intent(out):: elemFaceNeighbors(NELX,4)
      integer, intent(in) :: elemsInBox(:,:,:,:)
      integer, intent(in) :: nElemsInBox(:,:,:)
      integer, intent(in) :: nBoxes(3)
      integer, intent(in) :: IJK(4*NELX)
      real(8), intent(in) :: G0XYZ(NX*3),domainRange(3,2)
      integer, intent(in) :: NELX,NX,NODE
      !Local
      INTEGER:: I,J,II,JJ,neighFace
      integer:: neighCount, nodeCount
      integer:: NFACES
      logical :: neighborNode(NODE)
      real(8) :: xElemCenter(3), lenBoxes(3), domainLen(3)
      integer :: iBox(3), elemIdx

      !Initialize
      NFACES=4
      elemFaceNeighbors(:,:)=0
      neighFace=0
      
      domainLen(:)=domainRange(:,2)-domainRange(:,1)
      do I=1,3
      lenBoxes(I)=domainLen(I)/nBoxes(I)
      enddo
      
      DO I=1,NELX
         neighCount=0
         
         

         CALL getElemXCenter(I,IJK,G0XYZ,NODE,xElemCenter)
         ! find the virtual box that contains the point
         iBox(:) = (xElemCenter(:)-domainRange(:,1))/lenBoxes(:) + 1

         
         do elemIdx = 1,nElemsInBox(iBox(1),iBox(2),iBox(3))
            J = elemsInBox(elemIdx,iBox(1),iBox(2),iBox(3))

            IF(J.NE.I)THEN
               ! count the number of common nodes
               ! if 3, they're sharing a common face.
               nodeCount=0
               neighborNode(:)=.FALSE.
               DO II=1,NODE
                  DO JJ=1,NODE
                     IF(IJK((I-1)*NODE+II).EQ.IJK((J-1)*NODE+JJ))THEN
                        neighborNode(II)=.TRUE.
                        nodeCount=nodeCount+1
                     ENDIF
                  ENDDO
               ENDDO
               if(nodeCount.EQ.3)then
                  neighCount=neighCount+1
                  ! DETERMINE WHICH FACE of ELEMENT (I) 
                  ! the connection is through
                  ! then do elemFaceNeighbors(I,iFACE)=J
                  if(.NOT.neighborNode(1)) neighFace=1
                  if(.NOT.neighborNode(2)) neighFace=2
                  if(.NOT.neighborNode(3)) neighFace=3
                  if(.NOT.neighborNode(4)) neighFace=4
                  
                  elemFaceNeighbors(I,neighFace)=J
               endif
               
               !if number of neighbors found = number of faces. skip to the next element
               if(neighCount.EQ.NFACES) exit
                  
            ENDIF
         ENDDO
      ENDDO
      END SUBROUTINE

      SUBROUTINE calcEdgeNeighborsTET(elemEdgeNeighbors,IJK,NELX,NODE)
      ! arguments
      integer, intent(out):: elemEdgeNeighbors(NELX*6)
      integer, intent(in) :: IJK(4*NELX)
      integer, intent(in) :: NELX,NODE
      !Local
      INTEGER:: I,J,II,JJ
      integer:: neighCount, nodeCount
      integer:: NEDGES

      !Initialize
      NEDGES=6
      elemEdgeNeighbors(:)=0
      
      DO I=1,NELX
         neighCount=0
         DO J=1,NELX
            IF(J.NE.I)THEN
               ! count the number of common nodes
               ! if 2 or 3, they're sharing a common edge
               nodeCount=0
               DO II=1,NODE
                  DO JJ=1,NODE
                     IF(IJK((J-1)*NODE+JJ).EQ.IJK((I-1)*NODE+II))THEN
                        nodeCount=nodeCount+1
                     ENDIF
                  ENDDO
               ENDDO
               if(nodeCount.GE.2)then
                  neighCount=neighCount+1
                  elemEdgeNeighbors((I-1)*NEDGES+neighCount)=J
               endif
               
               !if number of neighbors found = number of edges. skip to the next element
               if(neighCount.EQ.NEDGES) exit
                  
            ENDIF
         ENDDO
      ENDDO
      END SUBROUTINE
      
      ! test subroutine for passing array arguments with assumed size
      subroutine subr(arr1)
      real(8), intent(out) :: arr1(:)
      integer :: sza
      
      sza = size(arr1,1)
      
      arr1(1:sza) = 2.0
      
      write(*,*) sza
      write(*,*) arr1
      
      end subroutine
      
      SUBROUTINE getElemXNodes(elemID,IJK,GXYZ,NODE,xNodes)
      implicit none
      integer, intent(in) :: elemID
      integer, intent(in) :: IJK(:)
      real(8), intent(in) :: GXYZ(:)
      integer, intent(in) :: NODE
      real(8), intent(out):: xNodes(3,NODE)
      integer :: NX_MAX,NDOF_MAX
      integer :: iNode,iGloNode    
      
      NX_MAX = size(IJK)
      NDOF_MAX = size(GXYZ)
      
      do iNode=1,NODE
         iGloNode=IJK((elemID-1)*NODE+iNode)
         xNodes(1:3,iNode)=GXYZ((iGloNode-1)*3+1:iGloNode*3)
      enddo
      
      END SUBROUTINE
      
      SUBROUTINE getElemFaceXNodes(elemID,iFace,IJK,GXYZ,NODE,xFaceNodes)
      implicit none
      integer, intent(in) :: elemID, iFace
      integer, intent(in) :: IJK(:)
      real(8), intent(in) :: GXYZ(:)
      integer, intent(in) :: NODE
      real(8), intent(out):: xFaceNodes(3,3)
      integer :: NX_MAX,NDOF_MAX
      integer :: iNode,iFaceNode,iGloNode    
      !Local
      real(8) :: xNodes(3,4)
      integer :: edgeV(4,3)
      real(8):: edge1(3),edge2(3)
            
      NX_MAX = size(IJK)
      NDOF_MAX = size(GXYZ)

      
      do iNode=1,NODE
         iGloNode=IJK((elemID-1)*NODE+iNode)
         xNodes(1:3,iNode)=GXYZ((iGloNode-1)*3+1:iGloNode*3)
      enddo
      
      iFaceNode=0
      do iNode = 1,NODE
         if (iNode==iFace) cycle ! exclude the apex node across the face
         iFaceNode=iFaceNode+1
         xFaceNodes(:,iFaceNode)=xNodes(:,iNode)
      enddo

      
      END SUBROUTINE
      
      SUBROUTINE getElemXCenter(elemID,IJK,GXYZ,NODE,xElemCenter)
      implicit none
      integer, intent(in) :: elemID
      integer, intent(in) :: IJK(:)
      real(8), intent(in) :: GXYZ(:)
      integer, intent(in) :: NODE
      real(8), intent(out):: xElemCenter(3)
      
      real(8) :: xElemNodes(3,NODE)
      integer :: NX_MAX,NDOF_MAX
      integer :: i
      
      NX_MAX = size(IJK)
      NDOF_MAX = size(GXYZ)
      
      xElemCenter(:)=0.D0
      CALL getElemXNodes(elemID,IJK,GXYZ,NODE,xElemNodes)
      do i=1,NODE
         xElemCenter(:) = xElemCenter(:) + xElemNodes(:,i)
      enddo
      xElemCenter(:) = xElemCenter(:) / NODE
      END SUBROUTINE
      
      SUBROUTINE calcElemFaceAreas(xNodes,faceAreas)
      implicit none
      real(8), intent(in) :: xNodes(3,4)
      real(8), intent(out):: faceAreas(4)
      
      integer :: iFace
      
      faceAreas(:)=0.D0
      
      do iFace = 1,4
         CALL calcFaceAreaTET(iFace,xNodes,faceAreas(iFace))
      enddo
      
      END SUBROUTINE

      SUBROUTINE calcFaceAreaTET(iFace,xNodes,faceArea)
      implicit none
      integer, intent(in) :: iFace
      real(8), intent(in) :: xNodes(3,4)
      real(8), intent(out):: faceArea
      ! locals
      real(8) :: v1(3), v2(3)
      integer :: edgeNodes(3)
      
      edgeNodes = 0
      
      if(iFace.EQ.1) then
         edgeNodes(1)=2
         edgeNodes(2)=3
         edgeNodes(3)=4
      elseif(iFace.EQ.2) then
         edgeNodes(1)=1
         edgeNodes(2)=3
         edgeNodes(3)=4
      elseif(iFace.EQ.3) then
         edgeNodes(1)=1
         edgeNodes(2)=2
         edgeNodes(3)=4
      elseif(iFace.EQ.4) then
         edgeNodes(1)=1
         edgeNodes(2)=2
         edgeNodes(3)=3
      endif   
      ! calculate the two edge vectors of the triangular face
      v1(:) = xNodes(:,edgeNodes(2))-xNodes(:,edgeNodes(1))
      v2(:) = xNodes(:,edgeNodes(3))-xNodes(:,edgeNodes(1))
      ! calculate the magnitude of the cross product, |v1 x v2|
      faceArea = (+v1(2)*v2(3)-v1(3)*v2(2))**2 +   &
                 (-v1(1)*v2(3)+v1(3)*v2(1))**2 +   &
                 (+v1(1)*v2(2)-v1(2)*v2(1))**2
      faceArea = 0.5 * SQRT(faceArea)
      END SUBROUTINE
      
      SUBROUTINE faceNormal(iFace, xElemNodes, dir)
      implicit none
      integer, intent(in) :: iFace
      real(8), intent(in) :: xElemNodes(3,4)
      integer, intent(out):: dir
      
      !Local
      real(8):: nFace(3)
      
      call faceNormalVector(iFace, xElemNodes, nFace)
      
      nFace(:)=abs(nFace(:))
      dir = MAXLOC(nFace(:),1)
      END SUBROUTINE

      
      SUBROUTINE faceNormalVector(iFace, xElemNodes, nFace, xFaceCenter)
      implicit none
      integer, intent(in) :: iFace
      real(8), intent(in) :: xElemNodes(3,4)
      real(8), intent(out):: nFace(3)
      real(8), intent(out), optional :: xFaceCenter(3)
      
      !Local
      integer :: edgeV(4,3)
      real(8):: edge1(3),edge2(3)
      
      edgeV(1,1) = 2
      edgeV(1,2) = 3
      edgeV(1,3) = 4
      
      edgeV(2,1) = 1
      edgeV(2,2) = 4
      edgeV(2,3) = 3
      
      edgeV(3,1) = 1
      edgeV(3,2) = 2
      edgeV(3,3) = 4
      
      edgeV(4,1) = 1
      edgeV(4,2) = 3
      edgeV(4,3) = 2

      edge1(:) = xElemNodes(:,edgeV(iFace,3)) - xElemNodes(:,edgeV(iFace,1))
      edge2(:) = xElemNodes(:,edgeV(iFace,2)) - xElemNodes(:,edgeV(iFace,1))
      
      CALL cross(nFace,edge1,edge2)
      
      if (present(xFaceCenter)) then
         xFaceCenter = (xElemNodes(:,edgeV(iFace,1)) + xElemNodes(:,edgeV(iFace,2)) + xElemNodes(:,edgeV(iFace,3))) / 3.d0
      endif
      
      END SUBROUTINE
      
      SUBROUTINE getLargestEdge(xNodes,maxEdgeLength)
      implicit none
      real(8), intent(in) :: xNodes(3,4)
      real(8), intent(out) :: maxEdgeLength
      real(8) :: edgeVector(3)
      
      maxEdgeLength=0.D0
      
      edgeVector(:)=xNodes(:,4)-xNodes(:,3)
      maxEdgeLength=MAX(maxEdgeLength,   &
       SQRT(edgeVector(1)**2.0+edgeVector(2)**2.0+edgeVector(3)**2.0))

      edgeVector(:)=xNodes(:,4)-xNodes(:,2)
      maxEdgeLength=MAX(maxEdgeLength,   &
       SQRT(edgeVector(1)**2.0+edgeVector(2)**2.0+edgeVector(3)**2.0))

      edgeVector(:)=xNodes(:,4)-xNodes(:,1)
      maxEdgeLength=MAX(maxEdgeLength,   &
       SQRT(edgeVector(1)**2.0+edgeVector(2)**2.0+edgeVector(3)**2.0))

      edgeVector(:)=xNodes(:,3)-xNodes(:,2)
      maxEdgeLength=MAX(maxEdgeLength,   &
       SQRT(edgeVector(1)**2.0+edgeVector(2)**2.0+edgeVector(3)**2.0))

      edgeVector(:)=xNodes(:,3)-xNodes(:,1)
      maxEdgeLength=MAX(maxEdgeLength,   &
       SQRT(edgeVector(1)**2.0+edgeVector(2)**2.0+edgeVector(3)**2.0))

      edgeVector(:)=xNodes(:,2)-xNodes(:,1)
      maxEdgeLength=MAX(maxEdgeLength,   &
       SQRT(edgeVector(1)**2.0+edgeVector(2)**2.0+edgeVector(3)**2.0))     
      END SUBROUTINE

      ! this subroutine splits the domain into virtual 'boxes', and assigns each element to the box that contains the element's geometric center
      ! this data structure is useful, when one needs to find an element at a particular region in space, by limiting the number of elements to search through.
      

      SUBROUTINE GETVOL(XNDS,VOL)
      !In
      REAL(8):: XNDS(3,4)
      !Out
      REAL(8):: VOL
      !Local
      REAL(8):: X1(3),X2(3),X3(3),X4(3)
      REAL(8):: XJACT(3,3)

      X1=XNDS(1:3,1)
      X2=XNDS(1:3,2)
      X3=XNDS(1:3,3)
      X4=XNDS(1:3,4)

      XJACT(1,:)=X1(1:3)-X4(1:3)
      XJACT(2,:)=X2(1:3)-X4(1:3)
      XJACT(3,:)=X3(1:3)-X4(1:3)

      VOL=XJACT(1,1)*(XJACT(2,2)*XJACT(3,3)-XJACT(2,3)*XJACT(3,2))
      VOL=VOL-XJACT(1,2)*(XJACT(2,1)*XJACT(3,3)-XJACT(2,3)*XJACT(3,1))
      VOL=VOL+XJACT(1,3)*(XJACT(2,1)*XJACT(3,2)-XJACT(2,2)*XJACT(3,1))
      VOL=VOL/6.D0
      END SUBROUTINE


      SUBROUTINE CHECKINSIDE(xElemNodes,xPoint,inside)
      !In
      REAL(8), intent(in) :: xElemNodes(3,4),xPoint(3)
      !Out
      logical, intent(out) :: inside
      !Local
      real(8) :: nFace(3), xFaceCenter(3)

      inside=.TRUE.  ! assume inside
      
      do iFace=1,4
      
         ! check if outside
         CALL faceNormalVector(iFace, xElemNodes, nFace, xFaceCenter)
         if (DOT_PRODUCT(xPoint-xFaceCenter,nFace) * &
             DOT_PRODUCT(xElemNodes(:,iFace)-xFaceCenter,nFace) < 0.d0) then   ! point and the apex are on different sides of the face
      
            inside=.FALSE.      ! mark as outside and leave
            exit
            
         endif
      enddo
      
      END SUBROUTINE
      
      SUBROUTINE assignElementsToBoxes(elemsInBox,nElemsInBox,nBoxes,   &
                                       domainRange, accordingToNodes,   &
                                       G0XYZ,NELX,NX,IJK,ierror)
                                       
      implicit none
      
      integer, intent(out):: elemsInBox(:,:,:,:)
      integer, intent(out):: nElemsInBox(:,:,:)
      integer, intent(in) :: nBoxes(3)
      real(8), intent(out) :: domainRange(3,2)
      logical, intent(in) :: accordingToNodes   ! or else, according to centroid
      real(8), intent(in) :: G0XYZ(NX*3)
      integer, intent(in) :: NELX,NX
      integer, intent(in) :: IJK(NELX*4)
      integer, intent(out):: ierror
      
      ! -- locals
      
      real(8) :: domainLen(3), tol(3), lenBoxes(3)
      real(8) :: xElemCenter(3), xElemNodes(3,4)
      
      integer :: iElem, iNode, iBox(3), newIndex, NODE
      integer :: elementAlreadyAssignedTo(4,3)
      
      integer :: MAX_ELEMS_IN_BOX
      integer :: MAX_NBOXES(3)
      MAX_ELEMS_IN_BOX = size(elemsInBox,1)
      MAX_NBOXES(1) = size(elemsInBox,2)
      MAX_NBOXES(2) = size(elemsInBox,3)
      MAX_NBOXES(3) = size(elemsInBox,4)

      NODE = 4
      ierror = 0
      
      if (MAX_NBOXES(1) .LT. nBoxes(1).OR. &
          MAX_NBOXES(2) .LT. nBoxes(2).OR. &
          MAX_NBOXES(3) .LT. nBoxes(3)) then
          
         ierror = 1  ! allocate more boxes
         RETURN
      endif

      ! -- implementation

      CALL getDomainRange(G0XYZ,NX,domainRange,domainLen)
      lenBoxes(:)=domainLen(:)/nBoxes(:)
      tol(:)=0.01*lenBoxes(:)
      
      nElemsInBox(:,:,:)=0
      elemsInBox(:,:,:,:)=0
      do iElem = 1,NELX
         elementAlreadyAssignedTo(:,:) = 0
         if(accordingToNodes) then
            CALL getElemXNodes(iElem,IJK,G0XYZ,NODE,xElemNodes)
            do iNode = 1,NODE
               ! -- find the virtual box that contains the point
               iBox(:) = (xElemNodes(:,iNode)-domainRange(:,1))/lenBoxes(:) + 1
               ! -- nodes on the boundary might go out of range. fix:
               if(iBox(1).GT.nBoxes(1)) iBox(1) = nBoxes(1)
               if(iBox(2).GT.nBoxes(2)) iBox(2) = nBoxes(2)
               if(iBox(3).GT.nBoxes(3)) iBox(3) = nBoxes(3)
               if(iBox(1).EQ.0) iBox(1) = 1
               if(iBox(2).EQ.0) iBox(2) = 1
               if(iBox(3).EQ.0) iBox(3) = 1
               ! -- add element to the box
               if(ALL(iBox(:).EQ.elementAlreadyAssignedTo(1,:)).OR. &
                  ALL(iBox(:).EQ.elementAlreadyAssignedTo(2,:)).OR. &
                  ALL(iBox(:).EQ.elementAlreadyAssignedTo(3,:)).OR. &
                  ALL(iBox(:).EQ.elementAlreadyAssignedTo(4,:))) then
                  ! this element already was assigned to this box. skip to next node.
                  cycle
               endif
               newIndex = nElemsInBox(iBox(1),iBox(2),iBox(3)) + 1
               if (newIndex.GT.MAX_ELEMS_IN_BOX) then
                  ierror = 2 ! allocate larger element list
                  return
               endif
               nElemsInBox(iBox(1),iBox(2),iBox(3)) = &
                  nElemsInBox(iBox(1),iBox(2),iBox(3)) + 1
               elemsInBox(newIndex,iBox(1),iBox(2),iBox(3)) = iElem
               elementAlreadyAssignedTo(iNode,:)=iBox(:)
            enddo
         else  !accordingToCentroid
            CALL getElemXCenter(iElem,IJK,G0XYZ,NODE,xElemCenter)
            ! find the virtual box that contains the point
            iBox(:) = (xElemCenter(:)-domainRange(:,1))/lenBoxes(:) + 1
            ! -- nodes on the boundary might go out of range. fix:
            if(iBox(1).GT.nBoxes(1)) iBox(1) = nBoxes(1)
            if(iBox(2).GT.nBoxes(2)) iBox(2) = nBoxes(2)
            if(iBox(3).GT.nBoxes(3)) iBox(3) = nBoxes(3)
            if(iBox(1).EQ.0) iBox(1) = 1
            if(iBox(2).EQ.0) iBox(2) = 1
            if(iBox(3).EQ.0) iBox(3) = 1
            ! -- add element to the box
            newIndex = nElemsInBox(iBox(1),iBox(2),iBox(3)) + 1
            if (newIndex.GT.MAX_ELEMS_IN_BOX) then
               ierror = 2 ! allocate larger element list
               return
            endif
            nElemsInBox(iBox(1),iBox(2),iBox(3)) = &
               nElemsInBox(iBox(1),iBox(2),iBox(3)) + 1
            elemsInBox(newIndex,iBox(1),iBox(2),iBox(3)) = iElem            
         endif
      enddo
      
      END SUBROUTINE
      
      SUBROUTINE getDomainRange(G0XYZ,NX,domainRange,domainLen)
      implicit none
      real(8), intent(in) :: G0XYZ(NX*3)
      integer, intent(in) :: NX
      real(8), intent(out) :: domainRange(3,2), domainLen(3)
      
      integer :: MDIM, I, iNode
      
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

      END SUBROUTINE

      SUBROUTINE getFacePositionsAVG(domainRange,domainLen,GXYZ,  &
                 faceNodes,nFaceNodes)
      implicit none
      real(8), intent(out) :: domainRange(3,2)
      real(8), intent(out) :: domainLen(3)
      real(8), intent(in) :: GXYZ(*)
      integer, intent(in) :: faceNodes(:,:)
      integer, intent(in) :: nFaceNodes(6)
      
      integer :: iNode,nodeIdx,iFace,iDOF,iDir
      integer :: szFN
      real(8) :: facePosMax, facePosMin
      
      szFN = size(faceNodes,1)
      
      do iDir = 1,3
      iFace = (iDir-1)*2+2
      facePosMax = 0.0
      do nodeIdx = 1,nFaceNodes(iFace)
         iNode = faceNodes(iFace,nodeIdx)
         iDOF = (iNode-1)*3+iDir
         facePosMax = facePosMax + GXYZ(iDOF)
      enddo
      facePosMax = facePosMax / nFaceNodes(iFace)

      iFace = (iDir-1)*2+1
      facePosMin = 0.0
      do nodeIdx = 1,nFaceNodes(iFace)
         iNode = faceNodes(iFace,nodeIdx)
         iDOF = (iNode-1)*3+iDir
         facePosMin = facePosMin + GXYZ(iDOF)
      enddo
      facePosMin = facePosMin / nFaceNodes(iFace)
      
      domainRange(iDir,1) = facePosMin
      domainRange(iDir,2) = facePosMax
      domainLen(iDir) = facePosMax -  facePosMin
         
      enddo   
      END SUBROUTINE
      
      SUBROUTINE findContainingElement(xPoint,elemID, &
                  NELX,NX,IJK,G0XYZ,                  &
                  elemsInBox,nElemsInBox,nBoxes,domainRange)
      
      implicit none
      
      real(8), intent(in) :: xPoint(3)
      integer, intent(out):: elemID
      integer, intent(in) :: NELX,NX
      integer, intent(in) :: IJK(NELX*4)
      real(8), intent(in) :: G0XYZ(NX*3)
      integer, intent(in) :: elemsInBox(:,:,:,:)
      integer, intent(in) :: nElemsInBox(:,:,:)
      integer, intent(in) :: nBoxes(3)
      real(8), intent(in) :: domainRange(3,2)

      ! -- locals
      
      integer :: iBox(3),elemIdx
      real(8) :: lenBoxes(3),domainLen(3)
      logical :: inside, found
      real(8) :: xNodes(3,4)
      integer :: NODE
      
      NODE = 4

      domainLen(:)=domainRange(:,2)-domainRange(:,1)
      lenBoxes(:)=domainLen(:)/nBoxes(:)
      
      iBox(:) = (xPoint(:)-domainRange(:,1))/lenBoxes(:) + 1
      
      found = .FALSE.
      do elemIdx = 1,nElemsInBox(iBox(1),iBox(2),iBox(3))
         elemID = elemsInBox(elemIdx,iBox(1),iBox(2),iBox(3))
         CALL getElemXNodes(elemID,IJK,G0XYZ,NODE,xNodes)

         CALL checkInside(xNodes,xPoint,inside)
         if (inside) then
            found = .TRUE.
            EXIT
         endif
      enddo
      
      if(.NOT.found) elemID = 0
      
      END SUBROUTINE
     
      ! calculates the element jacobian
      ! input: 
      !     xNodes(3,NODES)   : spatial positions of the nodes
      !  output:
      !     JAC(3,3)          : jacobian correspt. to natural->spatial mapping
      !     DVOL              : volume of the element
      SUBROUTINE calcJACTET(xNodes,JAC)
      implicit none

      REAL(8), intent(out):: JAC(3,3)
      REAL(8), intent(in) :: xNodes(3,4)
      
      real(8) :: XJACINV(3,3)
      real(8):: edge1(3),edge2(3),edge3(3)
      INTEGER:: I,J, INODE, IPT

      DO INODE=1,3
         edge1(:)=xNodes(:,INODE)-xNodes(:,4)
         edge2(:)=xNodes(:,INODE)-xNodes(:,4)
         edge3(:)=xNodes(:,INODE)-xNodes(:,4)
      ENDDO

      ! jacobian corresponding to the natural to spatial mapping of the element
      ! JAC_ij = d(x_i)/d(s_j)
      JAC(:,1)=edge1(:)
      JAC(:,2)=edge2(:)
      JAC(:,3)=edge3(:)
     
      RETURN
      END SUBROUTINE
      
      SUBROUTINE calcJacTET4(xElemNodes, elemJacMat, elemJac, elemVol)
      implicit none
      real(8), intent(in) :: xElemNodes(3,4)
      real(8), intent(out):: elemJacMat(3,3), elemJac, elemVol
      real(8):: XX(3),YY(3),ZZ(3)
      integer :: iCoord, iDir

      do iCoord=1,3
         do iDir=1,3
            elemJacMat(iCoord,iDir) =     &
               xElemNodes(iCoord,iDir)-xElemNodes(iCoord,4)
         enddo
      enddo
      
      elemJac= &
          elemJacMat(1,1)*elemJacMat(2,2)*elemJacMat(3,3)   &
         -elemJacMat(1,1)*elemJacMat(2,3)*elemJacMat(3,2)   &
         -elemJacMat(2,1)*elemJacMat(1,2)*elemJacMat(3,3)   &
         +elemJacMat(2,1)*elemJacMat(1,3)*elemJacMat(3,2)   &
         +elemJacMat(3,1)*elemJacMat(1,2)*elemJacMat(2,3)   &
         -elemJacMat(3,1)*elemJacMat(1,3)*elemJacMat(2,2)
     
      elemVol = elemJac/6
      
      END SUBROUTINE
      
      SUBROUTINE calcElemVol(elemID, elemVol, IJK, GXYZ, NODE)
      implicit none
      real(8), intent(out):: elemVol
      integer, intent(in) :: elemID
      integer, intent(in) :: IJK(:)
      real(8), intent(in) :: GXYZ(:)
      integer, intent(in) :: NODE
      
      real(8) :: xElemNodes(3,4), elemJacMat(3,3), elemJac
      integer :: iCoord, iDir
      
      CALL getElemXNodes(elemID,IJK,GXYZ,NODE,xElemNodes)

      do iCoord=1,3
         do iDir=1,3
            elemJacMat(iCoord,iDir) =     &
               xElemNodes(iCoord,iDir)-xElemNodes(iCoord,4)
         enddo
      enddo
      
      elemJac= &
          elemJacMat(1,1)*elemJacMat(2,2)*elemJacMat(3,3)   &
         -elemJacMat(1,1)*elemJacMat(2,3)*elemJacMat(3,2)   &
         -elemJacMat(2,1)*elemJacMat(1,2)*elemJacMat(3,3)   &
         +elemJacMat(2,1)*elemJacMat(1,3)*elemJacMat(3,2)   &
         +elemJacMat(3,1)*elemJacMat(1,2)*elemJacMat(2,3)   &
         -elemJacMat(3,1)*elemJacMat(1,3)*elemJacMat(2,2)
     
      elemVol = elemJac/6      
      
      END SUBROUTINE
      
      SUBROUTINE calcGradTensorTET(jacInvElem,nodalValuesTensor,gradTensor)
      implicit none
      real(8), intent(in) :: jacInvElem(3,3)
      real(8), intent(in) :: nodalValuesTensor(3,3,4)
      real(8), intent(out) :: gradTensor(3,3,3)
      real(8) :: gradN(3,4)
      integer :: iNode,i,j,k,l
      gradN(:,:) = 0.D0
      gradN(1,1)=1.D0
      gradN(2,2)=1.D0
      gradN(3,3)=1.D0
      gradN(1:3,4)=-1.D0
      
      gradTensor(:,:,:)=0.D0
      do i=1,3
         do j=1,3
            do k=1,3
               do iNode=1,4   !contract
                  do l=1,3    !contract
      gradTensor(i,j,k) = gradTensor(i,j,k)+ &
                     gradN(l,iNode)*jacInvElem(l,k)* &
                     nodalValuesTensor(i,j,iNode)
                  enddo
               enddo
            enddo
         enddo
      enddo
      
      END SUBROUTINE

!-----------------------------------------------------------------------
      SUBROUTINE calcCurlTensor(gradTensor,curlTensor)
      ! OLD CONVENTION:
      !  (See D.S.Chandrasekharaiah, 1990, Stokes Type Theorem for Tensors)
      !  CROSS PRODUCT OF vector v with 2nd ORDER TENSOR A IS DEFINED AS:
      !  (v x A).c =def= v x ((A^T).c)    for any vector, c.
      !  this leads to:
      !  (v x A)_ij = e_irs v_r A_js
      ! NEW CONVENTION: cross product isolates and hits on the first basis vector
      !  [a x (bc)]_ij = [(a x b)c]_ij = e_irs a_r b_s c_j
      implicit none
      real(8), intent(in) :: gradTensor(3,3,3)
      real(8), intent(out):: curlTensor(3,3)
      real(8) :: alter(3,3,3)
      integer::  i,j,k,r,s
      
      alter(:,:,:)=0.0d0
      
      alter(1,2,3)=1.0d0
      alter(2,3,1)=1.0d0
      alter(3,1,2)=1.0d0
      alter(2,1,3)=-1.0d0
      alter(1,3,2)=-1.0d0
      alter(3,2,1)=-1.0d0

      curlTensor(:,:)=0.0d0
      
      do i=1,3
         do j=1,3
            do s=1,3
               do r=1,3
      ! IN OLD CONVENTION: curlTensor(i,j)=curlTensor(i,j)+alter(i,r,s)*gradTensor(j,s,r)
      curlTensor(i,j)=curlTensor(i,j)+alter(i,r,s)*gradTensor(s,j,r)
               enddo
            enddo
         enddo
      enddo
    
      return
      END SUBROUTINE
      
      SUBROUTINE calcFrobeniusNorm(matrix,NCoeffs,norm)
      real(8), intent(in) :: matrix(NCoeffs)
      integer, intent(in) :: NCoeffs
      real(8), intent(out) :: norm
      integer :: i
      norm = 0.D0
      do i = 1,NCoeffs
         norm = norm + matrix(i)*matrix(i)
      enddo
      norm = sqrt(norm)
      END SUBROUTINE
      
      SUBROUTINE calcDistance(xNode1, xNode2, distance)
      implicit none
      real(8), intent(in) :: xNode1(3), xNode2(3)
      real(8), intent(out) :: distance
      
      distance = DSQRT(DOT_PRODUCT(xNode1(1:3)-xNode2(1:3),xNode1(1:3)-xNode2(1:3)))
            
      END SUBROUTINE
      
      SUBROUTINE cross(vc, v1, v2)
      implicit none
      real(8), intent(in) :: v1(3),v2(3)
      real(8), intent(out):: vc(3)
      
      vc(1)=+v1(2)*v2(3)-v1(3)*v2(2)
      vc(2)=-v1(1)*v2(3)+v1(3)*v2(1)
      vc(3)=+v1(1)*v2(2)-v1(2)*v2(1)
      
      END SUBROUTINE
      

      end module meshTools