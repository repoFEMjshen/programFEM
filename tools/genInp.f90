      PROGRAM generateMesh
      implicit none
      real(8), parameter :: PI = 3.14159265359D0
      
      integer, parameter :: nslip=30, nslip_fam=5
      real(8) :: cmt(3,nslip),cst(3,nslip),sch(3,3,nslip)
      real(8) :: gmt(3,nslip),gst(3,nslip)
      real(8) :: rot(3,3)
      integer :: nslips(nslip_fam), isys
      integer :: dir
      
      real(8) :: Lx,Ly,Lz
      real(8) :: dx,dy,dz
      real(8) :: stackX, stackY, stackZ
      real(8) :: dummyReal, dummyTexture(3)
      integer :: nNx,nNy,nNz
      integer :: inclusionSizeBlocks, nInclusions
      integer :: inclusionGrainID, inclusionHardGrainID, inclusionSoftGrainID
      integer :: offsetFromCenter
      logical :: breakHardGrainIntoSmallerGrains
      logical :: conventionABQ_CMRL
      character :: charConventionABQ_CMRL
      character :: elemType !T: TET B:BRICK
      character :: problemType, relaxHalfFace, bicrystalLayerDirection, bicryOri
      integer :: NX,NELX,nNodes,nElemsPerBlock,nBlocks,nGrains,iGrain,nGrainsMedium
      real(8), allocatable :: grainTexture(:,:)
      real(8), allocatable :: grainSizes(:)
      integer :: xShift, yShift, zShift
      integer :: nDim
      integer :: iElem,nElem_sta,iOri,iVar
      integer :: iBlock,jBlock,kBlock
      integer :: i,j,k,iNode,iFace,ind
      integer, allocatable :: tBlockIndices(:,:)
      integer, allocatable :: elemPhase(:), elemGrain(:)
      integer :: nInitStateVars
      real(8) :: initStateVars(100)
      integer :: iMat,nMat,nMatProp
      real(8) :: defaultMatProp(50,8)
      integer :: matVar,iLine,nLine,nRemain
      character :: singleOri
      real(8) :: Ori0(3) = (/7.32999966D-02,0.717800021D+00,2.44009995D+00/)
      real(8) :: Ori1(3) = (/0.00000000D+00,1.57079632679D+00,0.52359878D+00/) 
      real(8) :: Ori2(3) = (/0.00000000D+00,0.00000000D+00,0.00000000D+00/)
      real(8) :: Ori_softHCP(3)   =(/0.0D+00,  PI/4.0,  PI/6.0/)       ! if HCP, soft for basal <a> tension in the Z-dir
      real(8) :: Ori_hardHCP(3)   =(/0.0D+00, 0.0D+00, 0.0D+00/)       ! upright. if HCP, hard for tension in the Z direction
      real(8) :: Ori_LeftTilt(3)  =(/0.0D+00, -PI/4.0,  PI/6.0/)
      real(8) :: Ori_RightTilt(3) =(/0.0D+00,  PI/4.0,  PI/6.0/)
      real(8) :: Ori_TiltOnSide(3)=(/0.0D+00,  PI/4.0, 0.0D+00/)       ! HPC tilted 45 deg on its side
      real(8) :: Ori_OnSide(3)    =(/0.0D+00,  PI/2.0,  PI/6.0/)       ! HPC on its side soft for prismatic <a>
      real(8) :: Ori_OnSideFlat(3)=(/0.0D+00,  PI/2.0, 0.0D+00/)       ! HPC on its side, prismatic plane on ground
      real(8) :: Ori_single(3)
      
      character :: creep_or_CSR
      integer :: loadAxis
      character(len=100) :: strTextureFile
      real(8), allocatable :: textureEntries(:,:)
      integer :: nTextureEntries, iTextureEntry
      real(8) :: randomNumber
      logical :: fileExists
      integer :: error
      
      integer, allocatable :: boundaryNodes(:,:) ! stores nodeIDs given (faceID,nodeIndx)
      integer :: nBoundaryNodes(6)               ! nodes on each face
      integer, allocatable :: nodeSetNodes(:,:)
      integer, allocatable :: nNodeSetNodes(:)
      integer :: iNodeSet,nNodeSets
      integer :: II, JJ
      integer :: getNodeIDbyCoords
      integer, allocatable :: FACE_ELS_1(:),FACE_ELS_2(:)
      integer :: N_FACE_ELS
      
      character(len=32) :: meshName
      
      CALL init_random_seed()
               
      if (iargc().GE.1) then
         CALL getarg(1, meshName)
      else
         write(*,*) 'Enter mesh filename:'
         read(*,*) meshName
      endif
      CALL system('mkdir '//TRIM(meshName))
      CALL system('cp ./options.inp '//TRIM(meshName)//'/options.inp')
      CALL system('cp ./props_Ti6242.dat '//TRIM(meshName)//'/props_Ti6242.dat')

      write(*,*) 'Enter domain dimensions (Lx, Ly, Lz):'
      read(*,*) Lx,Ly,Lz
      write(*,*)
      
      write(*,*) 'Enter number of elements along each direction (Nx, Ny, Nz):'
      read(*,*) nNx,nNy,nNz
      write(*,*)
      dx=Lx/nNx
      dy=Ly/nNy
      dz=Lz/nNz
      nNx=nNx+1
      nNy=nNy+1
      nNz=nNz+1
      NX=nNx*nNy*nNz
      
      write(*,*) 'Node spacings (dx, dy, dz):'
      write(*,'(3F10.6)') dx,dy,dz
      write(*,*) '# of elements:', NX
      write(*,*)
      

      charConventionABQ_CMRL=' '
      do while (charConventionABQ_CMRL.NE.'A'.AND.charConventionABQ_CMRL.NE.'D')
         write(*,*) 'what convention to use for numbering nodes/element faces?'
         write(*,*) 'A: ABAQUS/CMRL'
         write(*,*) 'D: Deniz''s code'
         read(*,'(A1)') charConventionABQ_CMRL
         CALL toUpperCase(charConventionABQ_CMRL) !convert to uppercase
      enddo
      if (charConventionABQ_CMRL=='A') conventionABQ_CMRL=.true.
      if (charConventionABQ_CMRL=='D') conventionABQ_CMRL=.false.
      
      creep_or_CSR=' '
      do while (creep_or_CSR.NE.'C'.AND.creep_or_CSR.NE.'S' &
           .AND.creep_or_CSR.NE.'M')
         write(*,*) 'Specify loading type:'
         write(*,*) 'C: Creep or traction-controlled cyclic loading.'
         write(*,*) 'S: Constant Strain Rate, or displacement '// &
                    'controlled cyclic loading:'
         write(*,*) 'M: minimal BCs. only suppress RBMs.'
         read(*,'(A1)') creep_or_CSR
         CALL toUpperCase(creep_or_CSR) !convert to uppercase
      enddo
      write(*,*)
      
      loadAxis=0
      do while (.not.(loadAxis >= 1 .and. loadAxis <= 3))
         write(*,*) 'load axis: 1=X, 2=Y, 3=Z.'
         read(*,*) loadAxis
      enddo
      write(*,*)
      
      open(101,file=trim(meshName)//'/loadtime.inp')
      write(101,*) '//specify whether cyclic (1) or creep/constant strain rate(2)'
      write(101,*) 2
      write(101,*) '//specify file name'
      write(101,*) meshName
      close(101)

      open(102,file=trim(meshName)//'/loadtime_creep.inp')
      if(creep_or_CSR.EQ.'C') then
         write(102,*) '//tramp'
         write(102,*) 10.d0
         write(102,*) '//creep (1), CSR (2), tract cyclic/history(3/5), displ. cyclic/history conrolled(4/6)'      
         write(102,*) 1
         write(102,*) '//screep/strain rate'
         write(102,*) 730
      else
         write(102,*) '//tramp'
         write(102,*) 0.d0
         write(102,*) '//creep (1), CSR (2), tract cyclic/history(3/5), displ. cyclic/history conrolled(4/6)'      
         write(102,*) 2
         write(102,*) '//screep/strain rate'
         write(102,*) 0.0005
      endif
      close(102)

      open(103,file=trim(meshName)//'/slu_proc_distrib.inp')
      write(103,*) 2, 2
      close(103)
      
      elemType='T'
      do while (elemType.NE.'B'.AND.elemType.NE.'T')
         write(*,*) 'Specify element type. (T)ET, or (B)RICK:'
         read(*,'(A1)') elemType
         CALL toUpperCase(elemType) !convert to uppercase
      enddo
      write(*,*)
      if (elemType.EQ.'b') elemType='B'
      if (elemType.EQ.'t') elemType='T'
      
      problemType=' '
      do while (problemType.NE.'S'.AND.problemType.NE.'I'.AND.problemType.NE.'B' &
           .AND.problemType.NE.'P'.AND.problemType.NE.'H')
         write(*,*) '(S)ingle crystal'
         write(*,*) '(I)nclusion of hard grain at center'
         write(*,*) '(B)icrystal'
         write(*,*) '(P)air of soft-hard grains embedded in a medium of brick-shaped grains'
         write(*,*) '(H)ard inclusion embedded in a medium of brick-shaped grains'
         read(*,'(A1)') problemType
         CALL toUpperCase(problemType) !convert to uppercase
      enddo
      write(*,*)
      
      if (problemType.EQ.'S') then
         singleOri=' '
         do while (singleOri.NE.'U'.AND.singleOri.NE.'S'.AND.singleOri.NE.'P' &
              .AND.singleOri.NE.'T'.AND.singleOri.NE.'E')
            write(*,*) 'Orientation of the single HCP grain?:'
            write(*,*) '(U)pright (i.e., hard)'
            write(*,*) '(S)oft (i.e., tilted 45 deg, basal corner down)'
            write(*,*) '(T)ilted 45 deg, prismatic plane down)'
            write(*,*) 'laid on the (P)rismatic plane'
            write(*,*) 'user defined (E)uler angles'
            read(*,'(A1)') singleOri
            CALL toUpperCase(singleOri) !convert to uppercase
         enddo
         if(singleOri.EQ.'U') Ori_single(:)=Ori_hardHCP(:)
         if(singleOri.EQ.'S') Ori_single(:)=Ori_softHCP(:)
         if(singleOri.EQ.'T') Ori_single(:)=Ori_TiltOnSide(:)
         if(singleOri.EQ.'P') Ori_single(:)=Ori_OnSideFlat(:)
         if(singleOri.EQ.'E') then
            write(*,*) 'Enter euler angles [phi 1, theta, phi 2] (in degrees) for the single crystal:'
            read(*,*) Ori_single(1:3)
            Ori_single(:) = Ori_single(:)*2.D0*PI/360.D0
            write(*,*)
         endif
      elseif(problemType.EQ.'I') then
         bicryOri=' '
         do while (bicryOri.NE.'D'.AND.bicryOri.NE.'E')
            write(*,*) 'Crystal orientations of the inner and outer grains:'
            write(*,*) '(D)efault hard/soft orientation'
            write(*,*) 'user defined (E)uler angles'
            read(*,'(A1)') bicryOri
            CALL toUpperCase(bicryOri) !convert to uppercase
         enddo
         write(*,*)
         if(bicryOri.EQ.'E')then
            write(*,*) 'Enter euler angles [phi 1, theta, phi 2] (in degrees) for the inner (''hard'') grain:'
            read(*,*) Ori_hardHCP(1:3)
            write(*,*) 'Enter euler angles [phi 1, theta, phi 2] (in degrees) for the outer (''soft'') grain:'
            read(*,*) Ori_softHCP(1:3)
            Ori_hardHCP(:) = Ori_hardHCP(:)*2.D0*PI/360.D0
            Ori_softHCP(:) = Ori_softHCP(:)*2.D0*PI/360.D0
         endif
         
            write(*,*) 'Laterally offset the hard grain from the center?'
            write(*,*) 'No: 0, otherwise: enter distance to center in number of blocks'
            read(*,*) offsetFromCenter
      elseif(problemType=='P' .or. problemType=='H') then
         write(*,*) 'enter the size of the inclusion(s) in number of blocks:'
         read(*,*) inclusionSizeBlocks
         
         write(*,*) 'break hard grain into smaller grains? (T/F)'
         read(*,*) breakHardGrainIntoSmallerGrains
         
         fileExists = .false.
         do while (.not.fileExists)
            write(*,*) 'to assign orientations to the surrounding medium, we need a list of euler angles.'
            write(*,*) 'enter file name:'
            write(*,*) '(no header, no grain ID column, list of euler angles, 3 entries per line):'
            read(*,*) strTextureFile
            inquire(file=trim(strTextureFile),exist=fileExists)
         enddo
         
         nTextureEntries=0
         open(120,file=trim(strTextureFile))
         do
            read(120,*,iostat=error) dummyTexture(1:3)
            if(error.EQ.-1) EXIT
            nTextureEntries=nTextureEntries+1
         enddo
         close(120)
         allocate(textureEntries(3,nTextureEntries))
         textureEntries(:,:) = 0.d0
         iTextureEntry = 0
         open(120,file=trim(strTextureFile))
         do
            read(120,*,iostat=error) textureEntries(:,iTextureEntry+1)
            if(error.EQ.-1) EXIT
            iTextureEntry=iTextureEntry+1
         enddo
         close(120)
         write(*,*) nTextureEntries, ' texture entries imported for surrounding medium.'
      endif
      
      ! calculate and print Schmidt factors      
      call inicrys_hcp(cst,cmt,nslips)
      if(problemType=='I') then
         open(120,file=trim(meshName)//'/schmidt.txt')
         ! INNER / HARD / INCLUSION GRAIN
         call getRotationMatrix(Ori_hardHCP(1),Ori_hardHCP(2),Ori_hardHCP(3),rot)
         do isys=1,nslip
            gst(:,isys) = MATMUL(rot,cst(:,isys))  ! slip directions
            gmt(:,isys) = MATMUL(rot,cmt(:,isys))  ! slip plane normals
         enddo
         call create_schmidt(sch,gst,gmt,nslips)
         write(120,*) 'FOR THE INNER GRAIN/INCLUSION'
         do dir=1,3
            call print_schmidt(sch,dir)
         enddo
         write(*,*) 'Schmidt Factors saved to file: schmidt.txt'
         
         ! OUTER / SOFT / MATRIX GRAIN
         call getRotationMatrix(Ori_softHCP(1),Ori_softHCP(2),Ori_softHCP(3),rot)
         do isys=1,nslip
            gst(:,isys) = MATMUL(rot,cst(:,isys))  ! slip directions
            gmt(:,isys) = MATMUL(rot,cmt(:,isys))  ! slip plane normals
         enddo
         call create_schmidt(sch,gst,gmt,nslips)
         write(120,*) 'FOR THE OUTER GRAIN/MATRIX'
         do dir=1,3
            call print_schmidt(sch,dir)
         enddo
         write(*,*) 'Schmidt Factors saved to file: schmidt.txt'
         close(120)
      elseif(problemType=='S') then
         open(120,file=trim(meshName)//'/schmidt.txt')
         call getRotationMatrix(Ori_single(1),Ori_single(2),Ori_single(3),rot)
         do isys=1,nslip
            gst(:,isys) = MATMUL(rot,cst(:,isys))  ! slip directions
            gmt(:,isys) = MATMUL(rot,cmt(:,isys))  ! slip plane normals
         enddo
         call create_schmidt(sch,gst,gmt,nslips)
         write(120,*) 'SINGLE'
         do dir=1,3
            call print_schmidt(sch,dir)
         enddo
         write(*,*) 'Schmidt Factors saved to file: schmidt.txt'      
         close(120)
      endif
      


      relaxHalfFace='N'
      if (problemType.EQ.'S') then
         relaxHalfFace='N'
         do while (relaxHalfFace.NE.'Y'.AND.relaxHalfFace.NE.'N')
            write(*,*) 'Partially relieve displ. boundary conditions on top and bottom surfaces? Y/N:'
            write(*,*) 'This will create stress gradients'
            read(*,'(A1)') relaxHalfFace
            CALL toUpperCase(relaxHalfFace) !convert to uppercase
         enddo
         write(*,*)
      endif
      
      open(104,file=trim(meshName)//'/'//trim(meshName)//'.inp')
      write(104,*) '//DIMENSION,ADDITIONAL DOF PER NODE'
      write(104,*) 3, 0
                  
! ***** write node information
      write(104,*) '//NODE DATA'
      write(104,*) NX
      iNode=0
      do k=1,nNz
         do j=1,nNy
            do i=1,nNx
               iNode=iNode+1
               write(104,'(I5,3F20.10)') iNode, (i-1)*dx,(j-1)*dy,(k-1)*dz
            enddo
         enddo
      enddo
      
! ***** write element-node connectivity information

      if(elemType.EQ.'T') then
         nNodes = 4
         nElemsPerBlock = 6 ! 6 TETs per block
         NELX=(nNx-1)*(nNy-1)*(nNz-1)*nElemsPerBlock
         nBlocks = (nNx-1)*(nNy-1)*(nNz-1)

         allocate(tBlockIndices(nElemsPerBlock,nNodes))
         ! the following matrix indicates the positioning of
         ! tetrahedral elements inside blocks
         if (conventionABQ_CMRL) then
            tBlockIndices(1,1)=0                !    O   4
            tBlockIndices(1,2)=1                !  O   O |   1Z-,3X+
            tBlockIndices(1,3)=nNx+1            !        3
            tBlockIndices(1,4)=nNx+nNx*nNy+1    !  1---2
            tBlockIndices(2,1)=0                !    O   3
            tBlockIndices(2,2)=1                !  4   O
            tBlockIndices(2,3)=nNx+nNx*nNy+1    !  |     O   2Y-
            tBlockIndices(2,4)=nNx*nNy          !  1---2
            tBlockIndices(3,1)=1                !        4
            tBlockIndices(3,2)=nNx*nNy          !  2---3     1Y-,3Z+,4X+
            tBlockIndices(3,3)=nNx*nNy+1        !      | O
            tBlockIndices(3,4)=nNx+nNx*nNy+1    !      1
            tBlockIndices(4,1)=0                !    3   O
            tBlockIndices(4,2)=nNx              !  O | O     2Z-,1X-,3Y+
            tBlockIndices(4,3)=nNx+nNx*nNy      !    2   4
            tBlockIndices(4,4)=nNx+1            !  1   O
            tBlockIndices(5,1)=0                !    2---4
            tBlockIndices(5,2)=nNx+nNx*nNy      !  3
            tBlockIndices(5,3)=nNx*nNy          !  |         1X-, 3Z+
            tBlockIndices(5,4)=nNx+nNx*nNy+1    !  1
            tBlockIndices(6,1)=0                !    3---4
            tBlockIndices(6,2)=nNx+1            !        |   3Y+
            tBlockIndices(6,3)=nNx+nNx*nNy      !        2
            tBlockIndices(6,4)=nNx+nNx*nNy+1    !  1
         else ! deniz/oldCMRL convention
            tBlockIndices(1,1)=0                !    O   3
            tBlockIndices(1,2)=1                !  O   O |   3Z-,1X+
            tBlockIndices(1,3)=nNx+nNx*nNy+1    !        4
            tBlockIndices(1,4)=nNx+1            !  1---2
            tBlockIndices(2,1)=0                !    O   4
            tBlockIndices(2,2)=1                !  3   O
            tBlockIndices(2,3)=nNx*nNy          !  |     O   4Y-
            tBlockIndices(2,4)=nNx+nNx*nNy+1    !  1---2
            tBlockIndices(3,1)=1                !        3
            tBlockIndices(3,2)=nNx*nNy          !  2---4     3Y-,1Z+,2X+
            tBlockIndices(3,3)=nNx+nNx*nNy+1    !      | O
            tBlockIndices(3,4)=nNx*nNy+1        !      1
            tBlockIndices(4,1)=0                !    4   O
            tBlockIndices(4,2)=nNx              !  O | O     4Z-,3X-,1Y+
            tBlockIndices(4,3)=nNx+1            !    2   3
            tBlockIndices(4,4)=nNx+nNx*nNy      !  1   O
            tBlockIndices(5,1)=0                !    2---3
            tBlockIndices(5,2)=nNx+nNx*nNy      !  4
            tBlockIndices(5,3)=nNx+nNx*nNy+1    !  |         3X-, 1Z+
            tBlockIndices(5,4)=nNx*nNy          !  1
            tBlockIndices(6,1)=0                !    4---3
            tBlockIndices(6,2)=nNx+1            !        |   1Y+
            tBlockIndices(6,3)=nNx+nNx*nNy+1    !        2
            tBlockIndices(6,4)=nNx+nNx*nNy      !  1
         endif

      elseif(elemType.EQ.'B') then
         nNodes = 8
         nElemsPerBlock = 1
         NELX=(nNx-1)*(nNy-1)*(nNz-1)*nElemsPerBlock
         nBlocks = (nNx-1)*(nNy-1)*(nNz-1)

         allocate(tBlockIndices(nElemsPerBlock,nNodes))      
         ! positioning of 1 brick element in an 8-node block:
         tBlockIndices(1,1)=0
         tBlockIndices(1,2)=1
         tBlockIndices(1,3)=3
         tBlockIndices(1,4)=4
         tBlockIndices(1,5)=5
         tBlockIndices(1,6)=6
         tBlockIndices(1,7)=7
         tBlockIndices(1,8)=8

      endif
                  
      write(104,*) '//ELEM DATA'      
      write(104,*) NELX
      
      iElem=0
      !loop over the 8-node blocks of elements
      do k=1,nNz-1
         do j=1,nNy-1
            do i=1,nNx-1
               iNode=(k-1)*nNx*nNy+(j-1)*nNx+i
               do II=1,nElemsPerBlock
                  iElem=iElem+1
                  write(104,'(10I8)') iElem, &
                     ((tBlockIndices(II,JJ)+iNode),JJ=1,nNodes)
               enddo
            enddo
         enddo
      enddo

! ***** NODESETS
      nBoundaryNodes(1)=nNy*nNz
      nBoundaryNodes(2)=nNx*nNz
      nBoundaryNodes(3)=nNx*nNy
      nBoundaryNodes(4)=nNy*nNz
      nBoundaryNodes(5)=nNx*nNz
      nBoundaryNodes(6)=nNx*nNy
      allocate(boundaryNodes(6,MAXVAL(nBoundaryNodes(:))))
      boundaryNodes(:,:)=0
      ! nodes on x-faces (pulled if loadAxis=1)
      iNode=0
      do k=1,nNz
         do j=1,nNy
            iNode=iNode+1
            boundaryNodes(1,iNode)=getNodeIDbyCoords(1,  j,k,nNx,nNy,nNz)
            boundaryNodes(4,iNode)=getNodeIDbyCoords(nNx,j,k,nNx,nNy,nNz)
         enddo
      enddo
      ! nodes on y-faces (pulled if loadAxis=2)
      iNode=0
      do k=1,nNz
         do i=1,nNx
            iNode=iNode+1
            boundaryNodes(2,iNode)=getNodeIDbyCoords(i,1,  k,nNx,nNy,nNz)
            boundaryNodes(5,iNode)=getNodeIDbyCoords(i,nNy,k,nNx,nNy,nNz)
         enddo
      enddo
      ! nodes on z-faces (pulled if loadAxis=3)
      iNode=0
      do j=1,nNy
         do i=1,nNx
            iNode=iNode+1
            boundaryNodes(3,iNode)=getNodeIDbyCoords(i,j,1,  nNx,nNy,nNz)
            boundaryNodes(6,iNode)=getNodeIDbyCoords(i,j,nNz,nNx,nNy,nNz)
         enddo
      enddo
      
      !assign nodeSets <-- top & bottom faces
      nNodeSets=2
      allocate(nodeSetNodes(nNodeSets,MAXVAL(nBoundaryNodes(:))))
      allocate(nNodeSetNodes(nNodeSets))
      if(relaxHalfFace.EQ.'Y') then
         nNodeSetNodes(1)=nNx*((nNy+1)/2)
         nNodeSetNodes(2)=nNx*((nNy+1)/2)
      else
         nNodeSetNodes(1)=nBoundaryNodes(loadAxis)
         nNodeSetNodes(2)=nBoundaryNodes(loadAxis+3)
      endif
      do iNode=1,nNodeSetNodes(1)
         nodeSetNodes(1,iNode)=boundaryNodes(loadAxis,iNode)
         nodeSetNodes(2,iNode)=boundaryNodes(loadAxis+3,iNode)
      enddo

      ! Default BCs:
      ! 2 nodesets: TOP and BOTTOM faces
      ! bottom face: keep u_z fixed
      ! top face: u_z = assigned by function
      ! necessary and sufficient BCs to restrain rigid body motions:
      ! node 1: keep u_x, u_y fixed
      ! node nNx: keep u_y fixed
      write(104,*) '//NODESETS'
      write(104,*) nNodeSets
      do iNodeSet=1,nNodeSets
         write(104,*) iNodeSet,0
         write(104,*) nNodeSetNodes(iNodeSet),  &
          (nodeSetNodes(iNodeSet,iNode),iNode=1,nNodeSetNodes(iNodeSet))
      enddo

! ***** element sets
      allocate(FACE_ELS_1((nNx-1)*(nNy-1)))
      allocate(FACE_ELS_2((nNx-1)*(nNy-1)))
      N_FACE_ELS = 0
      if (loadAxis==3) then
         do i=1,nNx-1
            do j=1,nNy-1
               N_FACE_ELS=N_FACE_ELS+1
               
               if (conventionABQ_CMRL) then
                  FACE_ELS_1(N_FACE_ELS)=((nNx-1)*(nNy-1)*(nNz-2)+(j-1)*(nNx-1)+(i-1))*6 + 3
                  FACE_ELS_2(N_FACE_ELS)=((nNx-1)*(nNy-1)*(nNz-2)+(j-1)*(nNx-1)+(i-1))*6 + 5
               else  ! deniz/oldCMRL
                  FACE_ELS_1(N_FACE_ELS)=((nNx-1)*(nNy-1)*(nNz-2)+(j-1)*(nNx-1)+(i-1))*6 + 3
                  FACE_ELS_2(N_FACE_ELS)=((nNx-1)*(nNy-1)*(nNz-2)+(j-1)*(nNx-1)+(i-1))*6 + 5
               endif            
            enddo
         enddo
      elseif (loadAxis==2) then
         do i=1,nNx-1
            do j=1,nNz-1
               N_FACE_ELS=N_FACE_ELS+1
               
               if (conventionABQ_CMRL) then
                  FACE_ELS_1(N_FACE_ELS)=((j-1)*(nNx-1)*(nNy-1)+(nNx-1)*(nNy-2)+(i-1))*6 + 4
                  FACE_ELS_2(N_FACE_ELS)=((j-1)*(nNx-1)*(nNy-1)+(nNx-1)*(nNy-2)+(i-1))*6 + 6
               else  ! deniz/oldCMRL
                  FACE_ELS_1(N_FACE_ELS)=((j-1)*(nNx-1)*(nNy-1)+(nNx-1)*(nNy-2)+(i-1))*6 + 4
                  FACE_ELS_2(N_FACE_ELS)=((j-1)*(nNx-1)*(nNy-1)+(nNx-1)*(nNy-2)+(i-1))*6 + 6
               endif            
            enddo
         enddo
      elseif (loadAxis==1) then
         do i=1,nNy-1
            do j=1,nNz-1
               N_FACE_ELS=N_FACE_ELS+1
               
               if (conventionABQ_CMRL) then
                  FACE_ELS_1(N_FACE_ELS)=((j-1)*(nNx-1)*(nNy-1)+(nNx-2)+(nNx-1)*(i-1))*6 + 1
                  FACE_ELS_2(N_FACE_ELS)=((j-1)*(nNx-1)*(nNy-1)+(nNx-2)+(nNx-1)*(i-1))*6 + 3
               else  ! deniz/oldCMRL
                  FACE_ELS_1(N_FACE_ELS)=((j-1)*(nNx-1)*(nNy-1)+(nNx-2)+(nNx-1)*(i-1))*6 + 1
                  FACE_ELS_2(N_FACE_ELS)=((j-1)*(nNx-1)*(nNy-1)+(nNx-2)+(nNx-1)*(i-1))*6 + 3
               endif            
            enddo
         enddo
               
      endif
      if(N_FACE_ELS.NE.(nNx-1)*(nNy-1)) write(*,*) 'assertion error ! ! !'
      write(104,*) '//ELSETS'      
      write(104,*) 2
      write(104,*) 1,0
      write(104,*) N_FACE_ELS,FACE_ELS_1(1:N_FACE_ELS)
      write(104,*) 2,0
      write(104,*) N_FACE_ELS,FACE_ELS_2(1:N_FACE_ELS)

! ***** load history
      write(104,*) '//TOTIME,IDT,MAXDT,MINDT,MAXNUMSTEPS'
      write(104,*) 40.0D0,0.1D0,3.0D0,0.05D0,15000000

! ***** MATERIAL PROPERTIES
      defaultMatProp(:,:)=0
      open(201,file='defaultMatProp.dat')
      read(201,*) nMatProp
      nLine=nMatProp/8
      nRemain=MOD(nMatProp,8)
      do iLine=1,nLine
         read(201,*) (defaultMatProp(iLine,II),II=1,8)
      enddo
      read(201,*) (defaultMatProp(nLine+1,II),II=1,nRemain)
      read(201,*) matVar
      close(201)

      write(104,*) '//MATERIAL PROPERTIES'
      nMat = 1
      do iMat=1,nMat
         write(104,*) iMat
         write(104,*) nMatProp
         do iLine=1,nLine
         write(104,*) (defaultMatProp(iLine,II),II=1,8)
         enddo
         write(104,*) (defaultMatProp(nLine+1,II),II=1,nRemain)         
         write(104,*) matVar
      enddo
      
! **** DISP BOUNDARY CONDITIONS
      if (creep_or_CSR.EQ.'S') then
         ! displacement controlled
         write(104,*) '//DISP BOUND CONDITIONS'
         write(104,*) 4
         write(104,*) 1,1,0,loadAxis,loadAxis  ! bottom face - constrain u_z
         write(104,*) 0
         write(104,*) 2,1,1,loadAxis,loadAxis  ! top face - program-controlled
         ! suppress lateral RGBs
         write(104,*) 1,0,0,1,3  ! origin: fix all
         write(104,*) 0,0,0
         if (loadAxis==1) write(104,*) nNx*(nNy-1)+1,0,0,3,3  ! node (0,L,0): fix u_z
         if (loadAxis==2) write(104,*) nNx,0,0,3,3  ! node (L,0,0): fix u_z
         if (loadAxis==3) write(104,*) nNx,0,0,2,2  ! node (L,0,0): fix u_y
         write(104,*) 0

      elseif (creep_or_CSR.EQ.'C') then
         ! traction controlled
         write(104,*) '//DISP BOUND CONDITIONS'
         write(104,*) 3
         write(104,*) 1,1,0,loadAxis,loadAxis  ! bottom face - constrain u_z
         write(104,*) 0
         ! suppress lateral RGBs
         write(104,*) 1,0,0,1,3  ! origin: fix all
         write(104,*) 0,0,0
         if (loadAxis==1) write(104,*) nNx*(nNy-1)+1,0,0,3,3  ! node (0,L,0): fix u_z
         if (loadAxis==2) write(104,*) nNx,0,0,3,3  ! node (L,0,0): fix u_z
         if (loadAxis==3) write(104,*) nNx,0,0,2,2  ! node (L,0,0): fix u_y
         write(104,*) 0
      elseif (creep_or_CSR.EQ.'M') then
         ! only suppress RGB
         write(104,*) '//DISP BOUND CONDITIONS: only RGMs suppressed'
         write(104,*) 3
         write(104,*) 1,0,0,1,3  ! node 1: fix u_x, u_y, u_z
         write(104,*) 0,0,0
         write(104,*) nNx,0,0,2,3  ! node nNx: fix u_y, u_z
         write(104,*) 0,0
         write(104,*) nNx*(nNy-1)+1,0,0,3,3  ! node at (0,Ly,0): fix u_z
         write(104,*) 0
      endif
! **** Print out BCs
      write(*,*) 'Applied default BCs:'
      write(*,*) 'bottom face: u_z fixed'
      write(*,*) 'top face: u_z controlled by femuel'
      write(*,*) 'node 1: u_x, u_y fixed'
      write(*,'(A,1I5,A)') ' node ',nNx,': u_y fixed'
      write(*,*) 

! **** TRACTION BCs
      write(104,*) '//FORCE BOUNDARY CONDITIONS'
      write(104,*) 0
      if(creep_or_CSR.EQ.'C')then   ! load controlled
         write(104,*) '//PRESSURE BOUNDARY CONDITIONS'
         write(104,*) 1
         if (conventionABQ_CMRL) then
            if (loadAxis==3) write(104,*) 2,1,1,2,3,3  ! pull both element sets 1 and 2 from face 3
            if (loadAxis==2) write(104,*) 2,1,1,2,3,3  ! pull both element sets 1 and 2 from face 3
            if (loadAxis==1) write(104,*) 2,1,1,2,3,4  ! pull element set 1 from face 3, elem set 2 from face 4
         else ! deniz/oldCMRL
            if (loadAxis==3) write(104,*) 2,1,1,2,3,3  ! pull both element sets 1 and 2 from face 3
            if (loadAxis==2) write(104,*) 2,1,1,2,3,3  ! pull both element sets 1 and 2 from face 3
            if (loadAxis==1) write(104,*) 2,1,1,2,3,4  ! pull element set 1 from face 3, elem set 2 from face 4
         endif
      else
         write(104,*) '//PRESSURE BOUNDARY CONDITIONS'
         write(104,*) 0
      endif
      write(104,*) '//PERIODIC BOUNDARY CONDITIONS'
      write(104,*) 0
      
      close(104)

      !set up element phases
      allocate(elemPhase(NELX))
      elemPhase(:)=1
      allocate(elemGrain(NELX))
      elemGrain(:)=0
      !allocate(elemOri(NELX,3))
      !elemOri(:,:)=0.D0
      ! generate orientation angles 
      ! according to user directives
      if (problemType.EQ.'S') then           ! SINGLE CRYSTAL
         ! # grains
         nGrains = 1
         allocate(grainTexture(3,nGrains))
         allocate(grainSizes(nGrains))
         ! grain crystal orientations
         grainTexture(:,1) = Ori_single(:)
         ! grain size
         grainSizes(1) = Lx
         ! assign elements
         elemGrain(1:NELX)=1
      elseif (problemType.EQ.'I') then       ! INCLUSION
         ! # grains
         nGrains=2
         allocate(grainTexture(3,nGrains))
         allocate(grainSizes(nGrains))
         ! grain crystal orientations
         grainTexture(:,1) = Ori_softHCP(:)  ! outer grain is soft, has id=1
         grainTexture(:,2) = Ori_hardHCP(:)  ! inner grain is hard, has id=2
         ! grain sizes
         grainSizes(1) = Lz/2.d0
         grainSizes(2) = Lz/2.d0
         ! assign elements
         elemGrain(1:NELX)=1
         ! for offsetting the hard grain from the center:
         xShift = 0  ! placement of the inclusions are s.t. their separation vector is perpendicular to the load axis
         yShift = 0
         zShift = 0
         if (loadAxis==1) yShift = offsetFromCenter
         if (loadAxis==2) zShift = offsetFromCenter
         if (loadAxis==3) xShift = offsetFromCenter
         ! assign elements to 'hard' inclusion (id=2)
         do kBlock=(nNz-1)/3+1+zShift,2*(nNz-1)/3+zShift
            do jBlock=(nNy-1)/3+1+yShift,2*(nNy-1)/3+yShift
               do iBlock=(nNx-1)/3+1+xShift,2*(nNx-1)/3+xShift
                  nElem_sta=((kBlock-1)*(nNx-1)*(nNy-1) &
                        +(jBlock-1)*(nNx-1)+iBlock-1)*nElemsPerBlock
                  do iElem=nElem_sta+1,nElem_sta+nElemsPerBlock
                     elemGrain(iElem)=2  ! inner/hard grain has grainID = 1
                  enddo
               enddo
            enddo
         enddo
      elseif (problemType.EQ.'B') then       !  BICRYSTAL
         bicrystalLayerDirection=' '
         stackX=1.0
         stackY=1.0
         stackZ=1.0
         do while (bicrystalLayerDirection.NE.'X'.AND. &
                   bicrystalLayerDirection.NE.'Y'.AND. &
                   bicrystalLayerDirection.NE.'Z')
            write(*,*) 'Indicate stacking direction for the bicrystal layer (X,Y,Z):'
            read(*,'(A1)') bicrystalLayerDirection
            CALL toUpperCase(bicrystalLayerDirection) !convert to uppercase
         enddo
         if (bicrystalLayerDirection.EQ.'X') stackX = 0.5D+00
         if (bicrystalLayerDirection.EQ.'Y') stackY = 0.5D+00
         if (bicrystalLayerDirection.EQ.'Z') stackZ = 0.5D+00
         ! # grains
         nGrains=2
         allocate(grainTexture(3,nGrains))
         allocate(grainSizes(nGrains))
         ! grain crystal orientations
         grainTexture(:,1) = Ori_LeftTilt(:)  ! top layer has id=1
         grainTexture(:,2) = Ori_RightTilt(:) ! bottom layer had id=2
         ! grain sizes
         grainSizes(1) = Lz/2.d0
         grainSizes(2) = Lz/2.d0
         ! assign elements to grains
         !top layer
         elemGrain(1:NELX) = 1
         !bottom layer
         do kBlock=1,NINT((nNz-1)*stackZ)
            do jBlock=1,NINT((nNy-1)*stackY)
               do iBlock=1,NINT((nNx-1)*stackX)
                  nElem_sta=((kBlock-1)*(nNx-1)*(nNy-1) &
                        +(jBlock-1)*(nNx-1)+iBlock-1)*nElemsPerBlock
                  do iElem=nElem_sta+1,nElem_sta+nElemsPerBlock
                     elemGrain(iElem)=2
                  enddo
               enddo
            enddo
         enddo
      elseif (problemType=='P' .or. problemType=='H') then
         if (problemType=='P') nInclusions = 2
         if (problemType=='H') nInclusions = 1
         ! # grains
         nGrainsMedium = nBlocks - nInclusions*(inclusionSizeBlocks**3)
         if (.not.breakHardGrainIntoSmallerGrains) then
            nGrains = nGrainsMedium + nInclusions
         else
            nGrains = nGrainsMedium + (nInclusions - 1) + (inclusionSizeBlocks**3)
         endif
         
         allocate(grainTexture(3,nGrains))
         allocate(grainSizes(nGrains))
         ! assign orientations to the surrounding medium
         do iGrain = 1,nGrainsMedium
            CALL RANDOM_NUMBER(randomNumber)                ! randomly draw orientations from list of euler angles
            randomNumber = randomNumber * nTextureEntries
            iTextureEntry = ceiling(randomNumber)
            grainTexture(:,iGrain) = textureEntries(:,iTextureEntry)
         enddo
         ! assign orientations to the soft/hard pair/inclusion
         grainTexture(:,nGrainsMedium+1) = Ori_softHCP(:)   ! soft grain of the pair
         do iGrain=nGrainsMedium+1,nGrains
            grainTexture(:,iGrain) = Ori_hardHCP(:) ! hard grain(s)/inclusion whether split into blocks or not
         enddo
         do iGrain = 1,nGrainsMedium
            CALL RANDOM_NUMBER(randomNumber)                ! randomly draw orientations from list of euler angles
            randomNumber = randomNumber * nTextureEntries
            iTextureEntry = ceiling(randomNumber)
            grainTexture(:,iGrain) = textureEntries(:,iTextureEntry)
         enddo
         ! grain sizes in the medium and of the inclusions
         grainSizes(1:nGrainsMedium) = ((Lx*Ly*Lz)/nBlocks)**(1.d0/3.d0)
         grainSizes(nGrainsMedium+1:nGrains) = (((Lx*Ly*Lz)/nBlocks)**(1.d0/3.d0))*inclusionSizeBlocks
         ! assign elements
         if (problemType=='P') then            

            elemGrain(1:NELX) = 0
            
            xShift = 0  ! placement of the inclusions are s.t. their separation vector is perpendicular to the load axis
            yShift = 0
            zShift = 0
            if (loadAxis==1) yShift = inclusionSizeBlocks / 2
            if (loadAxis==2) zShift = inclusionSizeBlocks / 2
            if (loadAxis==3) xShift = inclusionSizeBlocks / 2

            inclusionSoftGrainID = nGrainsMedium + 1
            inclusionHardGrainID = nGrainsMedium + 2
            iGrain = 0
            do kBlock=1,nNz-1
               do jBlock=1,nNy-1
                  do iBlock=1,nNx-1
                  
                  nElem_sta=((kBlock-1)*(nNx-1)*(nNy-1) &
                           +(jBlock-1)*(nNx-1)+iBlock-1)*nElemsPerBlock

            if (kBlock >= (nNz-1-inclusionSizeBlocks)/2+1 -zShift .and. &
                kBlock <= (nNz-1+inclusionSizeBlocks)/2   -zShift .and. &
                jBlock >= (nNy-1-inclusionSizeBlocks)/2+1 -yShift .and. &
                jBlock <= (nNy-1+inclusionSizeBlocks)/2   -yShift .and. &
                iBlock >= (nNx-1-inclusionSizeBlocks)/2+1 -xShift .and. &
                iBlock <= (nNx-1+inclusionSizeBlocks)/2   -xShift) then
                
               ! this block is within the soft grain
               do iElem=nElem_sta+1,nElem_sta+nElemsPerBlock
                  elemGrain(iElem) = inclusionSoftGrainID
               enddo
            elseif(kBlock >= (nNz-1-inclusionSizeBlocks)/2+1 +zShift .and. &
                   kBlock <= (nNz-1+inclusionSizeBlocks)/2   +zShift .and. &
                   jBlock >= (nNy-1-inclusionSizeBlocks)/2+1 +yShift .and. &
                   jBlock <= (nNy-1+inclusionSizeBlocks)/2   +yShift .and. &
                   iBlock >= (nNx-1-inclusionSizeBlocks)/2+1 +xShift .and. &
                   iBlock <= (nNx-1+inclusionSizeBlocks)/2   +xShift) then

               ! this block is within the hard grain   
               do iElem=nElem_sta+1,nElem_sta+nElemsPerBlock
                  elemGrain(iElem) = inclusionHardGrainID
               enddo
               if (breakHardGrainIntoSmallerGrains) inclusionHardGrainID = inclusionHardGrainID + 1
            else
            
               iGrain = iGrain + 1
               
               ! this block is within the surrounding medium
               do iElem=nElem_sta+1,nElem_sta+nElemsPerBlock
                  elemGrain(iElem) = iGrain
               enddo
               
            endif
            
                  enddo         
               enddo   
            enddo   
            
         elseif (problemType=='H') then

            elemGrain(1:NELX) = 0
            
            inclusionGrainID = nGrainsMedium + 1
            iGrain = 0
            do kBlock=1,nNz-1
               do jBlock=1,nNy-1
                  do iBlock=1,nNx-1
                  
                  nElem_sta=((kBlock-1)*(nNx-1)*(nNy-1) &
                           +(jBlock-1)*(nNx-1)+iBlock-1)*nElemsPerBlock
                  
            if (kBlock >= (nNz-1-inclusionSizeBlocks)/2+1 .and. &
                kBlock <= (nNz-1+inclusionSizeBlocks)/2   .and. &
                jBlock >= (nNy-1-inclusionSizeBlocks)/2+1 .and. &
                jBlock <= (nNy-1+inclusionSizeBlocks)/2   .and. &
                iBlock >= (nNx-1-inclusionSizeBlocks)/2+1 .and. &
                iBlock <= (nNx-1+inclusionSizeBlocks)/2) then
                
               ! this block is within the inclusion (hard)
               do iElem=nElem_sta+1,nElem_sta+nElemsPerBlock
                  elemGrain(iElem) = inclusionGrainID
               enddo
               if (breakHardGrainIntoSmallerGrains) inclusionGrainID = inclusionGrainID + 1
            else
            
               iGrain = iGrain + 1
               
               ! this block is within the surrounding medium
               do iElem=nElem_sta+1,nElem_sta+nElemsPerBlock
                  elemGrain(iElem) = iGrain
               enddo
               
            endif
            
                  enddo         
               enddo   
            enddo   

         endif        
         
      endif
      
      open(104,file=trim(meshName)//'/grainSizes.inp')
      write(104,*) nGrains
      do iGrain = 1,nGrains
         write(104,*) iGrain, grainSizes(iGrain)
      enddo
      close(104)

      open(104,file=trim(meshName)//'/grainPhases.inp')
      write(104,*) nGrains
      do iGrain = 1,nGrains
         write(104,*) iGrain, 1
      enddo
      close(104)

! ***** Generate grain texture file
      open(106,file=trim(meshName)//'/grainTexture.inp')
      write(106,*) nGrains
      do iGrain=1,nGrains
         write(106,*) iGrain,grainTexture(:,iGrain)
      enddo
      close(106)
      
! ***** Generate elementGrains.inp
      open(107,file=trim(meshName)//'/elementGrains.inp')
      write(107,*) NELX
      do iElem=1,NELX
         write(107,*) iElem, elemGrain(iElem)
      enddo
      close(107)
      
! ***** Generate Patches (each block of six elements is assigned to one patch)
      open(108,file=trim(meshName)//'/PATCH.inp')
      write(108,*) nBlocks, 6, 6
      do iBlock=1,nBlocks
         write(108,*) 6, (iElem,iElem=(iBlock-1)*6+1,(iBlock-1)*6+6)
      enddo
      close(108)
      
      deallocate(elemGrain)
      deallocate(tBlockIndices)
      deallocate(boundaryNodes)
      deallocate(nodeSetNodes)
      deallocate(nNodeSetNodes)
      deallocate(FACE_ELS_1,FACE_ELS_2)

      END PROGRAM
      
      function getNodeIDbyCoords(i,j,k,nNx,nNy,nNz) result(nodeID)
         implicit none
         integer, intent(in) :: i,j,k,nNx,nNy,nNz
         integer:: nodeID
         nodeID = (k-1)*nNx*nNy+(j-1)*nNx+i
      end function
      
      SUBROUTINE printMatrix(A,N,M,matName,printTo)
      implicit none
      integer, intent(in):: N,M
      real(8), intent(in):: A(N,M)
      character*(*),intent(in):: matName
      integer, intent(in):: printTo
      integer:: i,j
      if(printTo.GT.0) then
         write(printTo,*) matName,':'
      else
         write(*,*) matName,':'
      endif
      do i=1,N
         if(printTo.GT.0) then
            write(printTo,*) (A(i,j), j=1,M)
         else
            write(*,*) (A(i,j), j=1,M)
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
      if(printTo.GT.0) then
         write(printTo,*) matName,':'
      else
         write(*,*) matName,':'
      endif
      do i=1,N
         if(printTo.GT.0) then
            write(printTo,*) (A(i,j), j=1,M)
         else
            write(*,*) (A(i,j), j=1,M)
         endif
      enddo
      END SUBROUTINE
      
      SUBROUTINE toUpperCase(ch)
      implicit none
      character, intent(inout) :: ch
      integer :: ascii
      ascii = IACHAR(ch)
      if(ascii.GE.97.AND.ascii.LE.122) ch=ACHAR(ascii-32)
      END SUBROUTINE

      SUBROUTINE getNc_fromEuler(nx,ny,nz,phi1,theta,phi2)
      implicit none
      real(8), intent(out) :: nx,ny,nz
      real(8), intent(in) :: phi1,theta,phi2
      real(8), parameter :: PI=3.14159265359D0

      nx = sin(theta)*sin(phi1)
      ny = -sin(theta)*cos(phi1)
      nz = cos(theta)

      END SUBROUTINE

      SUBROUTINE getNas_fromEuler(na1,na2,na3,phi1,theta,phi2)
      implicit none
      real(8), intent(out) :: na1(3),na2(3),na3(3)
      real(8), intent(in) :: phi1,theta,phi2
      real(8), parameter :: PI=3.14159265359D0
      real(8) :: n0(3), rot(3,3)
      
      !CALL getRot(rot,phi1,theta,phi2)
      CALL getRotationMatrix(phi1,theta,phi2,rot)
      
      n0(1) = 1.D0
      n0(2) = 0.D0
      n0(3) = 0.D0      
      na1(:) = MATMUL(rot(:,:),n0(:))
      n0(1) = cos(2.D0/3.D0*PI)
      n0(2) = sin(2.D0/3.D0*PI)
      n0(3) = 0.D0
      na2(:) = MATMUL(rot(:,:),n0(:))
      n0(1) = cos(4.D0/3.D0*PI)
      n0(2) = sin(4.D0/3.D0*PI)
      n0(3) = 0.D0
      na3(:) = MATMUL(rot(:,:),n0(:))

      END SUBROUTINE
      
      subroutine getRotationMatrix(phi,theta,omega,tlgt)
      implicit double precision (a-h,o-z)  

!      include 'aba_param.inc'

      dimension tlg(3,3),tlgt(3,3)
      real(8) :: phi,theta,omega

      pi=4.d0*datan(1.d0)
  
      sp=dsin(phi)                      
      cp=dcos(phi)                     
      st=dsin(theta)                     
      ct=dcos(theta)                    
      so=dsin(omega)                    
      co=dcos(omega)   
      tlg(1,1)=co*cp-so*sp*ct
      tlg(1,2)=co*sp+so*ct*cp   
      tlg(1,3)=so*st   
      tlg(2,1)=-so*cp-sp*co*ct 
      tlg(2,2)=-so*sp+ct*co*cp
      tlg(2,3)=co*st
      tlg(3,1)=sp*st       
      tlg(3,2)=-st*cp       
      tlg(3,3)=ct

      tlgt = TRANSPOSE(tlg)

      return
      end   
      
      


!-------------------------------------------------------------------------

!    THIS SUBROUTINE INITIALISES THE CRYSTAL SLIP SYSTEMS, 
!    CALCULATES THE TRANSFORMATION MATRIX AND CONTAINS 
!    STATE VARIABLES

!----------------------------------------------------------
!     cmt - slip plane in crystal system
!     cst - slip directions in crystal system
!     tlg0 - transformation matrix (kalidindis paper)
!     islip - number of slip systems: HCP:30, bcc:48
!     nslips(ifam) - number of slip systems in each slip family. 
!           ifam: slip familty:
!  {0001}<11-20>  (3) basal, a
!  {10-10}<11-20> (3) prismatic, a
!  {10-1-1}<11-20>(6) 1st pyramidal, a
!  {10-11}<11-23>(12) 1st pyramidal, c+a
!  {11-22}<11-23> (6) 2nd pyramidal c+a
!------------------------------------------------------------------------

      subroutine inicrys_hcp(cst,cmt,nslips)
     
      implicit real*8 (a-h,o-z)

      parameter(nslip=30)
      parameter(nslip_fam=5)

      dimension cmt(3,nslip),cst(3,nslip),nslips(nslip_fam)
        
      r2=1.d0/2.d0
      r3=dsqrt(3.d0)
      ic=0

         !{0001}<11-20>(3) basal, a
         !{10-10}<11-20>(3) prismatic, a
         !{10-1-1}<11-20>(6)1st pyramidal, a
         !{10-11}<11-23>(12)1st pyramidal, c+a
         !{11-22}<11-23>(6) 2nd pyramidal c+a

      ! n slip systems for slip family
      nslips(1)=3 ! basal<a>
      nslips(2)=3 ! prismatic<a>
      nslips(3)=6 ! 1st pyramidal<a>
      nslips(4)=12! 1st pyramidal<c+a>
      nslips(5)=6 ! 2nd pyramidal<c+a>

!---basal slip system
         ! slip planes
         do is=1,3
            cmt(1,is)=  0.0
            cmt(2,is)=  0.0 
            cmt(3,is)=  1.0d0
         enddo

         ! slip directions
         cst(1,1)=  r2
         cst(2,1)= -r2*r3 
         cst(3,1)=  0.0
         
         cst(1,2)=  r2
         cst(2,2)=  r3*r2
         cst(3,2)=  0.0
      
         cst(1,3)= -1.0d0
         cst(2,3)=  0.0
         cst(3,3)=  0.0
         ic=ic+nslips(1)

!---    prismatic slip system
         cmt(1,1+ic)=  0.0
         cmt(2,1+ic)=  1.0d0 
         cmt(3,1+ic)=  0.0
         
         cmt(1,2+ic)= -r3*r2
         cmt(2,2+ic)=  r2
         cmt(3,2+ic)=  0.0
         
         cmt(1,3+ic)= -r3*r2
         cmt(2,3+ic)= -r2 
         cmt(3,3+ic)=  0.0
         
         
         cst(1,1+ic)=  1.0d0
         cst(2,1+ic)=  0.0
         cst(3,1+ic)=  0.0
         
         cst(1,2+ic)=  r2
         cst(2,2+ic)=  r3*r2
         cst(3,2+ic)=  0.0
         
         cst(1,3+ic)= -r2
         cst(2,3+ic)=  r3*r2
         cst(3,3+ic)=  0.0
         ic=ic+nslips(2)

!       pyramidal a slip
         bunbo=sqrt(4*(1.587)**2+3)
         bunshi=1.587
         cmt(1,1+ic)=0.0
         cmt(2,1+ic)=-2*bunshi/bunbo
         cmt(3,1+ic)=r3/bunbo
         
         cmt(1,2+ic)=r3*bunshi/bunbo
         cmt(2,2+ic)=-bunshi/bunbo
         cmt(3,2+ic)=r3/bunbo
         
         cmt(1,3+ic)= r3*bunshi/bunbo
         cmt(2,3+ic)= bunshi/bunbo
         cmt(3,3+ic)= r3/bunbo 
         
         cmt(1,4+ic)= 0.0
         cmt(2,4+ic)= 2*bunshi/bunbo
         cmt(3,4+ic)= r3/bunbo
         
         cmt(1,5+ic)=-r3*bunshi/bunbo
         cmt(2,5+ic)= bunshi/bunbo
         cmt(3,5+ic)= r3/bunbo
         
         cmt(1,6+ic)=-r3*bunshi/bunbo
         cmt(2,6+ic)=-bunshi/bunbo
         cmt(3,6+ic)= r3/bunbo
         

         cst(1,1+ic)= 1.0d0
         cst(2,1+ic)= 0.0
         cst(3,1+ic)= 0.0 
         
         cst(1,2+ic)= r2
         cst(2,2+ic)= r2*r3
         cst(3,2+ic)= 0.0        
         
         cst(1,3+ic)=-r2
         cst(2,3+ic)= r2*r3
         cst(3,3+ic)= 0.0 
         
         cst(1,4+ic)=-1.0d0
         cst(2,4+ic)= 0.0
         cst(3,4+ic)= 0.0
         
         cst(1,5+ic)=-r2
         cst(2,5+ic)=-r2*r3
         cst(3,5+ic)= 0.0  
         
         cst(1,6+ic)= r2
         cst(2,6+ic)=-r2*r3
         cst(3,6+ic)= 0.0
         ic=ic+nslips(3)

!       1st order <c+a> slip
       bunbo=sqrt(4.0*(1.587)**2+3.0)
         bunshi=1.587
          
       cmt(1,1+ic)=  0.0
       cmt(2,1+ic)=  -2*bunshi/bunbo 
       cmt(3,1+ic)=  r3/bunbo
       
       cmt(1,2+ic)=  r3*bunshi/bunbo
       cmt(2,2+ic)=  -bunshi/bunbo 
       cmt(3,2+ic)=  r3/bunbo
       
       cmt(1,3+ic)=  r3*bunshi/bunbo
       cmt(2,3+ic)=  bunshi/bunbo 
       cmt(3,3+ic)=  r3/bunbo
      
       cmt(1,4+ic)=  0.0
       cmt(2,4+ic)=  2*bunshi/bunbo 
       cmt(3,4+ic)=  r3/bunbo
       
       cmt(1,5+ic)=  -r3*bunshi/bunbo
       cmt(2,5+ic)=  bunshi/bunbo 
       cmt(3,5+ic)=  r3/bunbo
       
       cmt(1,6+ic)=  -r3*bunshi/bunbo
       cmt(2,6+ic)=  -bunshi/bunbo 
       cmt(3,6+ic)=  r3/bunbo
      
       cmt(1,7+ic)=  0.0
       cmt(2,7+ic)=  -2*bunshi/bunbo 
       cmt(3,7+ic)=  r3/bunbo
       
       cmt(1,8+ic)=  r3*bunshi/bunbo
       cmt(2,8+ic)=  -bunshi/bunbo 
       cmt(3,8+ic)=  r3/bunbo
       
       cmt(1,9+ic)=  r3*bunshi/bunbo
       cmt(2,9+ic)=  bunshi/bunbo 
       cmt(3,9+ic)=  r3/bunbo
       
       cmt(1,10+ic)=  0.0
       cmt(2,10+ic)=  2*bunshi/bunbo 
       cmt(3,10+ic)=  r3/bunbo
       
       cmt(1,11+ic)= - r3*bunshi/bunbo
       cmt(2,11+ic)=  bunshi/bunbo 
       cmt(3,11+ic)=  r3/bunbo
      
       cmt(1,12+ic)=  -r3*bunshi/bunbo
       cmt(2,12+ic)=  -bunshi/bunbo 
       cmt(3,12+ic)=  r3/bunbo
      
!-----------slip directions

       bunb=2.0*(sqrt((1.587)**2+1.0))
       bush=1.587
       cst(1,1+ic)=  1.0/bunb
       cst(2,1+ic)=  r3/bunb 
       cst(3,1+ic)=  2*bush/bunb

       cst(1,2+ic)=  -1.0/bunb
       cst(2,2+ic)=  r3/bunb 
       cst(3,2+ic)=  2*bush/bunb

       cst(1,3+ic)=  -2.0/bunb
       cst(2,3+ic)=  0.0 
       cst(3,3+ic)=  2*bush/bunb

       cst(1,4+ic)=  -1.0/bunb
       cst(2,4+ic)=  -r3/bunb 
       cst(3,4+ic)=  2*bush/bunb

       cst(1,5+ic)=  1.0/bunb
       cst(2,5+ic)=  -r3/bunb 
       cst(3,5+ic)=  2*bush/bunb

       cst(1,6+ic)=  2.0/bunb
       cst(2,6+ic)=  0.0 
       cst(3,6+ic)=  2*bush/bunb

       cst(1,7+ic)=  -1.0/bunb
       cst(2,7+ic)=  r3/bunb 
       cst(3,7+ic)=  2*bush/bunb

       cst(1,8+ic)=  -2.0/bunb
       cst(2,8+ic)=  0.0 
       cst(3,8+ic)=  2*bush/bunb

       cst(1,9+ic)=  -1.0/bunb
       cst(2,9+ic)=  -r3/bunb 
       cst(3,9+ic)=  2*bush/bunb

       cst(1,10+ic)=  1.0/bunb
       cst(2,10+ic)=  -r3/bunb 
       cst(3,10+ic)=  2*bush/bunb
      
       cst(1,11+ic)=  2.0/bunb
       cst(2,11+ic)=  0.0 
       cst(3,11+ic)=  2*bush/bunb

       cst(1,12+ic)=  1.0/bunb
       cst(2,12+ic)=  r3/bunb 
       cst(3,12+ic)=  2*bush/bunb

       ic=ic+nslips(4)

!       2nd order <c+a> slip
      bunb=2.0*(sqrt((1.587)**2+1.0))
      bush=1.587
      cmt(1,1+ic)=  bush/bunb
      cmt(2,1+ic)=  -r3*bush/bunb 
      cmt(3,1+ic)=  2.0/bunb
      
      cmt(1,2+ic)= 2*bush/bunb
      cmt(2,2+ic)=  0.0
      cmt(3,2+ic)=  2.0/bunb
      
      cmt(1,3+ic)=  bush/bunb
      cmt(2,3+ic)=  r3*bush/bunb 
      cmt(3,3+ic)=  2.0/bunb

      cmt(1,4+ic)=  -bush/bunb
      cmt(2,4+ic)=  r3*bush/bunb 
      cmt(3,4+ic)=  2.0/bunb
      
      cmt(1,5+ic)= -2.0*bush/bunb
      cmt(2,5+ic)= 0.0 
      cmt(3,5+ic)=  2.0/bunb
      
      cmt(1,6+ic)=  -bush/bunb
      cmt(2,6+ic)=  -r3*bush/bunb 
      cmt(3,6+ic)=  2.0/bunb
      
!----------------slip directions------------
      
      cst(1,1+ic)=  -1.0/bunb
      cst(2,1+ic)=  r3/bunb 
      cst(3,1+ic)=  2.0*bush/bunb
      
      cst(1,2+ic)=  -2.0/bunb
      cst(2,2+ic)=  0.0 
      cst(3,2+ic)=  2.0*bush/bunb
      
      cst(1,3+ic)=  -1.0/bunb
      cst(2,3+ic)=  -r3/bunb 
      cst(3,3+ic)=  2.0*bush/bunb
      
      cst(1,4+ic)=  1.0/bunb
      cst(2,4+ic)=  -r3/bunb 
      cst(3,4+ic)=  2.0*bush/bunb
      
      cst(1,5+ic)=  2.0/bunb
      cst(2,5+ic)=  0.0 
      cst(3,5+ic)=  2.0*bush/bunb
      
      cst(1,6+ic)=  1.0/bunb
      cst(2,6+ic)=  r3/bunb 
      cst(3,6+ic)=  2.0*bush/bunb
      ic=ic+nslips(5)

      return
      end

      subroutine create_schmidt(sch,cst,cmt,nslips)
      implicit none
      
      integer, parameter :: nslip=30, nslip_fam=5
      
      real(8),intent(out):: sch(3,3,nslip)
      real(8),intent(in) :: cmt(3,nslip),cst(3,nslip)
      integer,intent(in) :: nslips(nslip_fam)
      
      integer :: i,j,isys
      ! deniz - comment
      ! here the schmidt tensors are formed, using the slip direction and plane normal vectors
      sch=0.0d0      
      do isys=1,nslip
         do j=1,3
          do i=1,3
             sch(i,j,isys)=cst(i,isys)*cmt(j,isys)
          enddo
         enddo
      enddo

         
         
      end subroutine
      
      SUBROUTINE print_schmidt(sch,dir)
      implicit none

      integer, parameter :: nslip=30, nslip_fam=5
      real(8), intent(in) :: sch(3,3,nslip)
      integer, intent(in) :: dir

      character(len=3) :: strDir
      strDir = 'XYZ'

      write(120,*) 'Schmidt Factors for loading along direction ',strDir(dir:dir)
      write(120,*) '      <a>basal       <a>basal       <a>basal       <a>prismatic   <a>prismatic   <a>prismatic'
      write(120,'(6F15.3)') abs(sch(dir,dir,1)),abs(sch(dir,dir,2)),abs(sch(dir,dir,3)),  &  ! 3 <a>basals
                                 abs(sch(dir,dir,4)),abs(sch(dir,dir,5)),abs(sch(dir,dir,6))      ! 3 <a>prismatics
      write(120,*) 'MAX among <a>basals:', MAXVAL(abs(sch(dir,dir,1:3)))
      write(120,*) 'MAX among <a>prismatics:', MAXVAL(abs(sch(dir,dir,4:6)))
      write(120,*) 'MAX SCHMIDT FACTOR among all:', MAXVAL(abs(sch(dir,dir,1:6)))
      write(120,*) ''

      END SUBROUTINE