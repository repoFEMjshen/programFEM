      PROGRAM correctConventionCMRL
      
      implicit none
      
      include 'PARDIS.H'
      
      integer :: NX, NELX
      real(8) :: G0XYZ(MAXCRD)
      real(8) :: nodePositions(3,MAXCRD/3)
      integer :: IJK(MNELX)
      integer :: NSET,NNSET(MAXNSET),NLIST(MAXNS,MAXNSET)
      integer :: NESET,NELSET(MAXELSET),NELLIST(MAXELS,MAXELSET)
      integer :: nElementSetsCMRL,iElemSet,iNodeSet
      character(len=30) :: meshName
      character(len=15) :: strGrain
      character(len=30) :: strElementType
      integer :: nElems_oldConvention
      logical :: fileExists,convertConvention
      logical :: swapIDs
      real(8) :: euler(3)
      
      integer :: i,iDummy,iGrain,nElem,iElem,iNode,nodeID
      integer :: iEdge,elemIdx,staIndex,endIndex,iIndex
      integer :: nTotalElements
      integer :: nShiftElemSetIDs,nShiftNodeSetIDs
      real(8) :: sampleLength
      integer :: loadAxis
      integer :: idNode3,idNode4 ! oldCMRL --> CMRL/ABQ requires swapping nodes 3/4
      real(8) :: xNodes(3,4)
      real(8) :: elemVol
      
      integer :: error
      character(len=150) :: lineStr
      
      ! read in the mesh and connectivity
      CALL readMesh_CMRL(G0XYZ,IJK,NX,NELX,strElementType)
      !CALL readElementSets_CMRL(NESET,NELSET,NELLIST)
      
      write(*,*) 'Mesh imported:'
      write(*,*) '# of elements: ', NELX
      write(*,*) '# of nodes: ', NX
      
      ! check if the node numbering follows oldCMRL or CMRL/ABAQUS convention
      nElems_oldConvention = 0
      do iElem=1,NELX
         do iNode=1,4
            nodeID = IJK((iElem-1)*4+iNode)
            xNodes(:,iNode) = G0XYZ((nodeID-1)*3+1:(nodeID-1)*3+3)
         enddo
         call GETVOL_oldCMRL(xNodes,elemVol)
         if (elemVol > 0.d0) nElems_oldConvention = nElems_oldConvention + 1
      enddo
      swapIDs = .false.
      if (nElems_oldConvention == 0) then
         write(*,*) 'this mesh already has the newCMRL/ABAQUS convention for node numbering.'
      elseif (nElems_oldConvention == NELX) then
         write(*,*) 'this mesh has the oldCMRL convention for node numbering'
         write(*,*) 'connectivity and faceIDs will be swapped ...'
         swapIDs = .true.
      else
         write(*,*) 'WARNING!'
         write(*,*) 'node numbering in this mesh is does not follow a consistent convention'
         stop      
      endif
      
      ! CMRL/ABAQUS convention has the nodes 3/4 swapped wrt oldCMRL convention.
      ! (this also requires swapping faceIDs 1/2 for consistency -- see READ_WRITE_BC())
      if (swapIDs) then
         do iElem=1,NELX
            idNode3 = IJK((iElem-1)*4+3)
            idNode4 = IJK((iElem-1)*4+4)
            IJK((iElem-1)*4+3) = idNode4
            IJK((iElem-1)*4+4) = idNode3
         enddo
      endif

      ! create list of nodes
!      write(*,*) 'creating MS_nodes.inp ... from the main input file'
!      open(105,file='MS_nodes_corrected.inp')
!      write(105,*) NX
!      do iNode=1,NX
!         nodePositions(:,iNode) = G0XYZ((iNode-1)*3+1:(iNode-1)*3+3)
!         write(105,*) iNode, nodePositions(:,iNode)
!      enddo
!      close(105)
!      write(*,*) 'created MS_nodes_corrected.inp'
   
      ! create list of nodes
      write(*,*) 'creating MS_connectivity.inp ... from the main input file'
      open(105,file='MS_connectivity_corrected.inp')
      write(105,*) NELX, strElementType
      do iElem=1,NELX
         write(105,*) iElem, IJK((iElem-1)*4+1:(iElem-1)*4+4)
      enddo
      close(105)
      write(*,*) 'created MS_connectivity_corrected.inp'

      
      END PROGRAM

      

!*********************************************************
      subroutine readMesh_CMRL(G0XYZ,IJK,NX,NELX,strElementType)
      
      implicit none
      include 'PARDIS.H'
      real(8), intent(out) :: G0XYZ(MAXCRD)
      integer, intent(out) :: IJK(MNELX)
      integer, intent(out) :: NX, NELX
      character(len=*), intent(out) :: strElementType
      !locals
      integer :: IJKNEL(MAXELN)
      real(8) :: xc(3)
      integer :: MDIM,NDFADD,NDF,NNODE,NGAUSS
      integer :: error
      logical :: fileExists
      character(len=50) :: fileName
      character(len=150) :: lineStr
      !dummies
      integer :: icyc,len,i,ii,j,jj,nnum
      integer :: iElem,iElemDummy,iNode,iNodeDummy,NEL,NEQ
      
      INQUIRE(file='MS_nodes.inp',exist=fileExists)
      if(.NOT.fileExists) then
         write(*,*) 'file not found: MS_nodes.inp'
         stop
      endif

      ! create list of nodes
      write(*,*) 'reading node positions: MS_nodes.inp ...'
      open(105,file='MS_nodes.inp')
      read(105,*) NX
      do iNodeDummy=1,NX
         read(105,*) iNode, G0XYZ((iNode-1)*3+1:(iNode-1)*3+3)
      enddo
      close(105)
      write(*,*) 'done.'
      
      INQUIRE(file='MS_connectivity.inp',exist=fileExists)
      if(.NOT.fileExists) then
         write(*,*) 'file not found: MS_connectivity.inp'
         stop
      endif

      ! create list of nodes
      write(*,*) 'reading connectivity matrix: MS_connectivity.inp ...'
      open(105,file='MS_connectivity.inp')
      read(105,*) NELX, strElementType
      do iElemDummy=1,NELX
         read(105,*) iElem, IJK((iElem-1)*4+1:(iElem-1)*4+4)
      enddo
      close(105)
      write(*,*) 'done.'
      
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

      SUBROUTINE GETVOL_oldCMRL(XNDS,VOL)
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