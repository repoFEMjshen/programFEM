      PROGRAM generateMesh
      
      use meshtools
      
      implicit none
      real(8), parameter :: PI = 3.14159265358979323846264338327d0
      real(8) :: Lx,Ly,Lz
      
      logical :: DEBUG
      parameter(DEBUG=.TRUE.)
      
      integer :: NX_max, NELX_max, MAXNODELEM, N_HISTORY_MAX
      integer :: NNODES, NFACES
      integer :: MAX_ELEM_SETS,MAX_ELEMS_PER_SET
      integer :: N_GRAINS_MAX, maxGrainIDInABQ, iGrain_ABQ, iGrains_ABQ, nGrains_ABQ, nGrains_D3D, nGrains_removedD3D
      parameter(NX_max = 500000,NELX_max = 1000000, MAXNODELEM=120, N_HISTORY_MAX=50)
      parameter(NNODES=4,NFACES=4)
      parameter(MAX_ELEM_SETS=24,MAX_ELEMS_PER_SET=NELX_max/MAX_ELEM_SETS)
      parameter(N_GRAINS_MAX=1000)
      

      real(8) :: G0XYZ(3*NX_max)
      integer :: IJK(4*NELX_max)
      integer:: IELEMNO(MAXNODELEM*NX_max),NINDX(NX_max)
      real(8) :: xElemNodes(3,4)
      integer :: dir, dirDomainFace
      
      integer :: NX,NELX
      integer :: nDim
      integer :: iNode,nodeID,iElem,iGrain
      integer :: elemID,iFace,grainID,regionID,iDOF,iOri,iVar
      integer :: grainID_D3D
      real(8) :: grainSize_D3D(N_GRAINS_MAX), grainOri_D3D(3,N_GRAINS_MAX)
      integer :: firstInd,lastInd,iDomainFace,indElemSet,indElem
      integer :: i,j,k
      integer :: transB_ratio
      real(8) :: transB_rand
      
      ! vars on BCs/loading types
      integer :: nBCS                    ! number of BCs
      character::BCtype                  ! type of boundary conditions
      character::loadDir                 ! CSR loading direction
      character::periodicBC              ! impose periodic BCs Y/N
      character::thermLoad               ! thermal loading conditions
      integer :: N_temp_history, iTemp
      real(8) :: temp_value, temp_time
      real(8) :: temp_value_history(N_HISTORY_MAX), &
                 temp_time_history(N_HISTORY_MAX)
      character::temp_cycling_type       ! A: cyclic around initial temp, S: cyclic above init. temp
      real(8) :: temp_stress_free        ! initial/stress-free temperature
      real(8) :: temp_cyclic_amplitude   ! in Kelvins
      real(8) :: temp_cyclic_period      ! in seconds
      
      character :: elemType !T: TET B:BRICK      
      real(8), allocatable :: elemOri(:,:)
      integer, allocatable :: elemPhase(:)
      integer :: elemGrain(NELX_max)
      real(8) :: grainOri(3,N_GRAINS_MAX), grainSize(N_GRAINS_MAX)
      integer :: grainMapABQ_to_D3D(N_GRAINS_MAX)
      logical :: grainDetectedInABQ(N_GRAINS_MAX)
      
!      integer :: nInitStateVars
!      real(8) :: initStateVars(100)
      integer :: iMat,nMat,nMatProp
      real(8) :: defaultMatProp(50,8)
      integer :: matVar,iLine,nLine,nRemain

      real(8) :: Ori0(3) = (/7.32999966D-02,0.717800021D+00,2.44009995D+00/)
      real(8) :: Ori1(3) = (/0.00000000D+00,1.57079632679D+00,0.52359878D+00/) 
      real(8) :: Ori2(3) = (/0.00000000D+00,0.00000000D+00,0.00000000D+00/)
      real(8) :: Ori_softHCP(3)   =(/0.0D+00,  PI/4.0,  PI/6.0/)       ! if HCP, soft for tension in the Z-dir
      real(8) :: Ori_hardHCP(3)   =(/0.0D+00, 0.0D+00, 0.0D+00/)       ! upright. if HCP, hard for tension in the Z direction
      real(8) :: Ori_LeftTilt(3)  =(/0.0D+00, -PI/4.0,  PI/6.0/)
      real(8) :: Ori_RightTilt(3) =(/0.0D+00,  PI/4.0,  PI/6.0/)
      real(8) :: Ori_TiltOnSide(3)=(/0.0D+00,  PI/4.0, 0.0D+00/)       ! HPC tilted 45 deg on its side
      real(8) :: Ori_OnSide(3)    =(/0.0D+00,  PI/2.0,  PI/6.0/)       ! HPC on its side
      real(8) :: Ori_OnSideFlat(3)=(/0.0D+00,  PI/2.0, 0.0D+00/)       ! HPC on its side, prismatic plane on ground
      
      integer :: cornerNodes(8)
      real(8) :: domainRange(3,2)
      integer, allocatable:: elemFaceNeighbors(:,:)
      integer, allocatable :: boundaryNodes(:,:) ! stores nodeIDs given (faceID,nodeIndx)
      integer :: nBoundaryNodes(6)               ! nodes on each face
      integer, allocatable :: nodeSetNodes(:,:)
      integer, allocatable :: nNodeSetNodes(:)
      integer:: elemSetElems(MAX_ELEM_SETS,MAX_ELEMS_PER_SET)
      integer:: nElemSetElems(MAX_ELEM_SETS)
      logical, allocatable :: elemFacesMarked(:)
      integer :: MAXNODESonEachFace
      integer :: iNodeSet,nNodeSets,iElemSet,nElemSets
      integer :: II, JJ
      integer :: getNodeIDbyCoords
      
      character(len=50) :: dirINP, fileABQ, fileDream3DCSV, fileANSYS, fileCPFEM
      integer :: unitABQ, unitDream3DCSV, unitANSYS
      parameter(unitABQ = 10, unitDream3DCSV=11, unitANSYS=12)
      integer :: error
      logical :: fileExists, found, foundElements, foundNodes
      character(len=150) :: fieldsCSV
      character :: readingWhat   ! E: element list, N: node list, X:neither
      character(len=300) :: strLine
      character(len=20)  :: strWord
      
      CALL init_random_seed()
      
      write(*,*) 'Specify the name of the directory to look for:'
      write(*,*) '1) the Dream3D CSV output storing grain size and orientation info: size_euler.csv'
      write(*,*) '2) ABAQUS mesh file: ABQ.inp'
      write(*,*) '3) ANSYS mesh file: ANSYS.inp'
      read(*,'(A)') dirINP
      write(*,*)
      
      fileABQ = trim(dirINP)//'/ABQ.inp'
      fileDream3DCSV = trim(dirINP)//'/size_euler.csv'
      fileANSYS = trim(dirINP)//'/ANSYS.inp'
	  

      INQUIRE(file=fileABQ,exist=fileExists)
      if(.NOT.fileExists) then
         write(*,*) 'ABQ mesh not found: ',fileABQ
         STOP
      endif
      INQUIRE(file=fileANSYS,exist=fileExists)
      if(.NOT.fileExists) then
         write(*,*) 'ANSYS mesh not found: ',fileANSYS
         write(*,*) 'ANSYS mesh is needed to extract the grain-region mapping in simmetrix'
         STOP
      endif
      INQUIRE(file=fileDream3DCSV,exist=fileExists)
      if(.NOT.fileExists) then
         write(*,*) 'Dream3D CSV file not found: ',fileDream3DCSV
         STOP
      endif
      
      write(*,*) 'Choose a name for the CPFE model:'
      read(*,*) fileCPFEM
      write(*,*)

      CALL system('mkdir '//trim(fileCPFEM))
      CALL system('cp options.inp '//trim(fileCPFEM)//'/options.inp')      
      CALL system('cp props_Ti6242.dat '//trim(fileCPFEM)//'/props_Ti6242.dat')      
      CALL system('cp runfem.job '//trim(fileCPFEM)//'/runfem.job')      
      
      open(101,file=trim(fileCPFEM)//'/loadtime.inp')
      write(101,*) '//specify whether cyclic (1) or creep/constant strain rate(2)'
      write(101,*) 2
      write(101,*) '//specify file name'
      write(101,*) fileCPFEM
      close(101)

      open(103,file=trim(fileCPFEM)//'/slu_proc_distrib.inp')
      write(103,*) 2, 2
      close(103)
      
      ! IMPORT ABQ MESH
      ! FORM NX, NELX, G0XYZ, CONNECTIVITY
      NX=0
      NELX=0
      IJK(:)=0
      IELEMNO(:)=0
      NINDX(:)=0
      G0XYZ(:)=0.D0
      elemGrain(:)=0
      open(unit=unitABQ,file=trim(fileABQ))
      ! ---------- NODES -------------!
      ! seek the beginning of node data
      found=.FALSE.
      foundElements=.FALSE.
      foundNodes=.FALSE.
      readingWhat='X'
      grainDetectedInABQ(:) = .FALSE.
      maxGrainIDInABQ = 0
      nGrains_ABQ = 0
      
      do while(.TRUE.)
         read(unitABQ,'(A)',iostat=error) strLine
         if (error.NE.0) then
            if(foundElements.AND.foundNodes) then
               ! EOF. Got some elements and nodes. proceed.
               exit
            else
               write(*,*) 'ABQ file: unexpected end of file. These headers were not found:'
               if(.NOT.foundElements) write(*,*) '*Element'
               if(.NOT.foundNodes) write(*,*) '*Node'
               STOP
            endif
         endif
         CALL toUpperCase(strLine)
         if (strLine(1:5).EQ.'*NODE') then
            found=.TRUE.
            foundNodes=.TRUE.
            readingWhat='N'
            cycle
         endif
         if (strLine(1:8).EQ.'*ELEMENT') then
            found=.TRUE.
            foundElements=.TRUE.
            readingWhat='E'                        ! now reading elements
            CALL extractRegionID(strLine,regionID) ! for this ABQ grainID/regionID
            cycle
         endif
         if (readingWhat.EQ.'N') then
            NX=NX+1
            read(strLine, *) nodeID, G0XYZ((nodeID-1)*3+1:(nodeID-1)*3+3)         
            if (NX.GT.NX_max) then
               write(*,*) 'Increase NX_max!'
               STOP
            endif         
         elseif(readingWhat.EQ.'E') then
         ! got a new element connectivity line
            NELX=NELX+1
            read(strLine, *) iElem, IJK((NELX-1)*4+1), &
                                    IJK((NELX-1)*4+2), &
                                    IJK((NELX-1)*4+4), & ! conventions differ.
                                    IJK((NELX-1)*4+3)    ! switch nodes 3 and 4.
            elemGrain(NELX) = regionID
            maxGrainIDInABQ = MAX(regionID,maxGrainIDInABQ)
            if (.NOT.grainDetectedInABQ(regionID)) then
               nGrains_ABQ = nGrains_ABQ + 1
            endif
            grainDetectedInABQ(regionID) = .TRUE.
            if (NELX.GT.NELX_max) then
               write(*,*) 'Increase NELX_max!'
               STOP
            endif         
         endif
      enddo
      close(unitABQ)
      
      if (DEBUG) then
         write(*,*) '# of distinct grain IDs (regions) in the ABQ file:', nGrains_ABQ
         do grainID = 1,maxGrainIDInABQ
            if (grainDetectedInABQ(grainID)) then
               !write(*,'(I0,A)',advance='no') grainID,' '
            else
               !write(*,'(A)',advance='no') 'XX '
            endif
         enddo
         write(*,*) '' 
      endif

      write(*,'(A)',advance='no') ' constructing node-to-element connectivity table... '
      CALL NODE_ELEM_CONNEC(IELEMNO,NINDX,IJK,NX,NELX,4*NELX_max,MAXNODELEM)

      write(*,*) 'done.'
      CALL writeNODE_ELEM_CONNEC(trim(fileCPFEM)//'/nodeElementConnect.dat', &
                                 IELEMNO,NINDX,NELX,NX,MAXNODELEM,NX_max)
      write(*,*) 'saved connectivity data file: nodeElementConnect.dat'
      
      write(*,*) 'Nodes imported: ', NX
      write(*,*) 'Elements imported: ', NELX
      write(*,*) 'Element-node connectivity imported.'
      write(*,*)
            
      write(*,'(A)',advance='no') 'generating element-element connectivity table...'
      allocate(elemFaceNeighbors(NELX,NFACES))
      ! CALCULATE ELEMENT-ELEMENT CONNECTIVITY
      CALL calcFaceNeighborsTET(elemFaceNeighbors(:,:),IJK,IELEMNO,NINDX, & 
                                   NELX,NNODES,MAXNODELEM,NX_max)
      write(*,*) 'done.'
      CALL writeFaceNeighborsTET(trim(fileCPFEM)//'/elementNeighbors.dat', &
                              elemFaceNeighbors,NELX)
      write(*,*) 'saved connectivity data file: elementNeighbors.dat'
      write(*,*)                      
      write(*,*) 'Identifying domain boundaries...'      
      MAXNODESonEachFace = (1.0*NX)**(2.0/3.0)*50
      allocate(boundaryNodes(6,MAXNODESonEachFace))
      ! identify the nodes on each face of the cubic domain
      ! determine nodes at the corners of the domain
      ! determine the spatial range of the domain (DX, DY, DZ)
      CALL identifyNodes(boundaryNodes,nBoundaryNodes,cornerNodes, &
                          domainRange,G0XYZ,elemFaceNeighbors, &
                          IJK,NELX,NX,MAXNODESonEachFace)
      
      elemType='T'   ! assume TET

      open(104,file=trim(fileCPFEM)//'/'//trim(fileCPFEM)//'.inp')
      write(104,*) '//DIMENSION,ADDITIONAL DOF PER NODE'
      write(104,*) 3, 0
                  
! ***** write node information
      write(104,*) '//NODE DATA'
      write(104,*) NX
      do iNode=1,NX
         write(104,*) iNode, G0XYZ((iNode-1)*3+1:(iNode-1)*3+3)
      enddo
      
! ***** write element-node connectivity information
      write(104,*) '//ELEM DATA'      
      write(104,*) NELX
      do iElem=1,NELX
         write(104,*) iElem, IJK((iElem-1)*4+1:(iElem-1)*4+4)
      enddo
      
! ***** NODESETS: 6 domain faces: loX, hiX, loY, hiY, loZ, hiZ

      ! Default BCs:
      ! 2 nodesets: TOP and BOTTOM faces
      ! bottom face: keep u_z fixed
      ! top face: u_z = assigned by function
      ! necessary and sufficient BCs to restrain rigid body motions:
      ! node 1: keep u_x, u_y fixed
      ! node nNx: keep u_y fixed
      nNodeSets = 6 + 8
      write(104,*) '//NODESETS'
      write(104,*) nNodeSets
      do iNodeSet=1, 6           ! nodes on the boundary faces
         write(104,*) iNodeSet,0
         write(104,*) nBoundaryNodes(iNodeSet),  &
          (boundaryNodes(iNodeSet,iNode),iNode=1,nBoundaryNodes(iNodeSet))
      enddo
      do iNodeSet=7, nNodeSets   ! corner nodes. 7=(0,0,0), 8=(L,0,0), 9=(0,L,0), 10=(L,L,0) etc..
         write(104,*) iNodeSet,0
         write(104,*) 1, cornerNodes(iNodeSet-6)
      enddo

! ***** element sets
      ! elements on face 1 of the domain, with local face number 1 on the domain:
      !  ---> elemSet 1
      ! elements on face 1 of the domain, with local face number 2 on the domain:
      !  ---> elemSet 2
      ! elements on face 1 of the domain, with local face number 3 on the domain:
      !  ---> elemSet 3
      ! elements on face 1 of the domain, with local face number 4 on the domain:
      !  ---> elemSet 4
      ! elements on face 2 of the domain, with local face number 1 on the domain:
      !  ---> elemSet 5
      ! etc...
      
      allocate(elemFacesMarked(NELX))
      
      nElemSetElems(:)=0
      elemFacesMarked(:)=.FALSE.
      elemSetElems(:,:)=0
      
      if(DEBUG) CALL printIntMatrix(nBoundaryNodes,1,6,'number of nodes on each face:',0)
      
      do iDomainFace=1,6
      
         elemFacesMarked(:) = .FALSE.
         dirDomainFace = ((iDomainFace-1)/2)+1
         
         do iNode=1,nBoundaryNodes(iDomainFace)
            nodeID = boundaryNodes(iDomainFace,iNode)
            lastInd = NINDX(nodeID)-1
            firstInd = 1
            if(nodeID.NE.1) firstInd = NINDX(nodeID-1)
            do iElem=firstInd,lastInd
               elemID = IELEMNO(iElem)
               if(elemFacesMarked(elemID)) CYCLE
               elemFacesMarked(elemID) = .TRUE.
               do iFace=1,4
                  if(elemFaceNeighbors(elemID,iFace).EQ.0) then   ! this face is on a boundary
                     CALL getElemXNodes(elemID,IJK,G0XYZ,NNODES,xElemNodes)
                     CALL faceNormal(iFace,xElemNodes,dir)  ! dir: 1/2/3
                     if (.NOT.dir.EQ.dirDomainFace) Cycle   ! this face is not on iDomainFace, skip

                     ! add element to the proper element set
                     indElemSet = (iDomainFace-1)*4+iFace
                     indElem = nElemSetElems(indElemSet)+1
                     elemSetElems(indElemSet,indElem) = elemID
                     nElemSetElems(indElemSet) = indElem
                     
                     exit

                  endif
               enddo
            enddo
         enddo
      enddo
      
      write(104,*) '//ELSETS'      
      write(104,*) 24
      do iDomainFace=1,6
         do iFace=1,4
            indElemSet = (iDomainFace-1)*4+iFace
            write(104,*) indElemSet,0
            write(104,*) nElemSetElems(indElemSet),   &
                         elemSetElems(indElemSet,1:nElemSetElems(indElemSet))
         enddo
      enddo

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
      if (nRemain > 0) read(201,*) (defaultMatProp(nLine+1,II),II=1,nRemain)
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
      
      nBCS=3
      BCtype=' '
      do while (BCtype.NE.'C'.AND.BCtype.NE.'S'.AND.BCtype.NE.'M')
         write(*,*) 'Specify the boundary conditions:'
         write(*,*) 'C: Creep or traction-controlled cyclic loading.'
         write(*,*) 'S: Constant Strain Rate, or displacement '// &
                    'controlled cyclic loading:'
         write(*,*) 'M: minimal BCs. only suppress RBMs.'
         read(*,'(A1)') BCtype
         CALL toUpperCase(BCtype) !convert to uppercase
      enddo
      write(*,*)
      
      loadDir = 'Z'
      
      if(BCtype=='S'.or.BCtype=='C') then
         loadDir=' '
         do while (loadDir.NE.'X'.AND.loadDir.NE.'Y'.AND.loadDir.NE.'Z')
            write(*,*) 'Specify the deformation/loading direction'
            write(*,*) 'X, Y or Z ?'
            read(*,'(A1)') loadDir
            CALL toUpperCase(loadDir) !convert to uppercase
         enddo
         write(*,*)         
      endif
      
      periodicBC=' '
      do while (periodicBC.NE.'Y'.AND.periodicBC.NE.'N')
         write(*,'(A)',advance='no') ' Impose periodic BCs ? (Y/N) '
         read(*,'(A1)') periodicBC
         CALL toUpperCase(periodicBC)
      enddo

      thermLoad=' '
      do while (thermLoad.NE.'N'.AND.thermLoad.NE.'C' &
                                .AND.thermLoad.NE.'H')
         write(*,*) 'Specify thermal loading'
         write(*,*) 'N: None'
         write(*,*) 'C: Thermal cycles'
         write(*,*) 'H: Provide temperature-time history'
         
         read(*,'(A1)') thermLoad
         CALL toUpperCase(thermLoad)
      enddo

      write(*,'(A)',advance='no') 'Initial temperature (Kelvins): '
      read(*,*) temp_stress_free
      
      if (thermLoad.EQ.'C') then
         write(*,'(A)',advance='no') 'Thermal cycle amplitude (Kelvins): '
         read(*,*) temp_cyclic_amplitude
         write(*,'(A)',advance='no') 'Thermal cycle period (seconds): '
         read(*,*) temp_cyclic_period
         temp_cycling_type=' '
         do while (temp_cycling_type.NE.'A'.AND.temp_cycling_type.NE.'S')
            write(*,'(A)',advance='no') 'Cycle (A)round the initial temp, or (S)trictly above? '
            read(*,'(A1)') temp_cycling_type
            CALL toUpperCase(temp_cycling_type)
         enddo
      elseif(thermLoad.EQ.'H') then
         ! temperature-time history given
         write(*,*) 'enter time (t), temperature (T) pairs: t, T [enter]'
         write(*,*) 'when done: 0, 0 [enter]'

         N_temp_history=0
         temp_value = 1.0
         temp_time = 1.0
         do while (temp_value.GT.0.0.OR.temp_time.GT.0.0)
            read(*,*) temp_time, temp_value
            N_temp_history=N_temp_history+1
            temp_time_history(N_temp_history)=temp_time
            temp_value_history(N_temp_history)=temp_value
         enddo
         N_temp_history=N_temp_history-1
      endif
      
! **** DISP BOUNDARY CONDITIONS

      if (BCtype.EQ.'S') then
         nBCS=5 ! 3 for suppressing Free Body Modes, +2 for constraining bottom and top faces
      elseif (BCtype.EQ.'C') then
         nBCS=4 ! 3 for suppressing Free Body Modes, +1 for constraining bottom face
      elseif (BCtype.EQ.'M') then
         nBCS=3 ! 3 for suppressing Free Body Modes
      endif
      
      write(104,*) '//DISP BND CONDS. Ln1:#BC,[periodicBC(1:Yes,0:No)]'
      if(periodicBC.EQ.'N') then 
         write(104,*) nBCS
      else
         write(104,*) nBCS, 1 ! mark the second field to signal femuel to impose periodic BCs.
      endif
      ! suppressing free body modes, in any case..
      if (loadDir=='Z') then
         write(104,*) 6+1,1,0,1,3  ! corner node at (0,0,0). Fix uX, uY, uZ.
         write(104,*) 0,0,0
         write(104,*) 6+2,1,0,2,3  ! corner node at (L,0,0). Fix uY, uZ.
         write(104,*) 0,0
         write(104,*) 6+3,1,0,3,3  ! corner node at (0,L,0). Fix uZ.
         write(104,*) 0
      elseif (loadDir=='X') then
         write(104,*) 6+1,1,0,1,3  ! corner node at (0,0,0). Fix uX, uY, uZ.
         write(104,*) 0,0,0
         write(104,*) 6+5,1,0,1,2  ! corner node at (0,0,L). Fix uX, uY.
         write(104,*) 0,0
         write(104,*) 6+3,1,0,1,1  ! corner node at (0,L,0). Fix uX.
         write(104,*) 0
      elseif (loadDir=='Y') then
         write(104,*) 6+1,1,0,1,3  ! corner node at (0,0,0). Fix uX, uY, uZ.
         write(104,*) 0,0,0
         write(104,*) 6+2,1,0,2,3  ! corner node at (L,0,0). Fix uY, uZ.
         write(104,*) 0,0
         write(104,*) 6+5,1,0,2,2  ! corner node at (0,0,L). Fix uY.
         write(104,*) 0
      endif
      
      if (BCtype.EQ.'S') then    ! if displacement-controlled test:
         if  (loadDir=='Z') then
            write(104,*) 5,1,0,3,3  ! Fix uZ at the bottom face
            write(104,*) 0
            write(104,*) 6,1,1,3,3  ! and mark the top face to be controlled by the code.
         elseif (loadDir=='X') then
            write(104,*) 1,1,0,1,1  ! Fix uX at the -X face
            write(104,*) 0
            write(104,*) 2,1,1,1,1  ! and mark the +X face to be controlled by the code.
         elseif (loadDir=='Y') then
            write(104,*) 3,1,0,2,2  ! Fix uY at the -Y face
            write(104,*) 0
            write(104,*) 4,1,1,2,2  ! and mark the +Y face to be controlled by the code.
         endif
      elseif (BCtype.EQ.'C') then! if traction-controlled test:
         if  (loadDir=='Z') then
            write(104,*) 5,1,0,3,3  ! Fix uZ at the bottom face
            write(104,*) 0
         elseif (loadDir=='X') then
            write(104,*) 1,1,0,1,1  ! Fix uX at the -X face
            write(104,*) 0
         elseif (loadDir=='Y') then
            write(104,*) 3,1,0,2,2  ! Fix uY at the -Y face
            write(104,*) 0
         endif
      endif

! **** TRACTION BCs
      write(104,*) '//FORCE BOUNDARY CONDITIONS'
      write(104,*) 0
      
      if(BCtype.EQ.'C')then   ! load controlled
         write(104,*) '//PRESSURE BOUNDARY CONDITIONS'
         write(104,*) 1
         if(loadDir=='Z') then
            write(104,*) 4,1,21,22,23,24,3,4,2,1
         elseif(loadDir=='X') then
            write(104,*) 4,1,5,6,7,8,3,4,2,1
         elseif(loadDir=='Y') then
            write(104,*) 4,1,13,14,15,16,3,4,2,1
         endif
      else
         write(104,*) '//PRESSURE BOUNDARY CONDITIONS'
         write(104,*) 0
      endif
      write(104,*) '//PERIODIC BOUNDARY CONDITIONS'
      write(104,*) 0
      
      close(104)
      
! **** loading conditions
      open(102,file=trim(fileCPFEM)//'/loadtime_creep.inp')
      if(BCtype.EQ.'C') then
         write(102,*) '//tramp'
         write(102,*) 10.00
         write(102,*) '//creep (1), CSR (2), tract cyclic/history(3/5), displ. cyclic/history conrolled(4/6)'      
         write(102,*) 1
         write(102,*) '//screep/strain rate'
         write(102,*) 750
      else
         write(102,*) '//tramp'
         write(102,*) 0.d0
         write(102,*) '//creep (1), CSR (2), tract cyclic/history(3/5), displ. cyclic/history conrolled(4/6)'      
         write(102,*) 2
         write(102,*) '//screep/strain rate'
         write(102,*) 0.001
      endif
      close(102)      
      
! ***** Generate thermal_loading.inp
      open(109,file=trim(fileCPFEM)//'/thermal_loading.inp')
      write(109,*) '// thermal loading type. 1: provide temperature-time history &
                    2: thermal cycling around init. temp. 3: thermal cycling above init. temp'
      if(thermLoad.EQ.'H'.OR.thermLoad.EQ.'N') then
         write(109,*) 1
      elseif(thermLoad.EQ.'C') then
         if(temp_cycling_type.EQ.'A') write(109,*) 2
         if(temp_cycling_type.EQ.'S') write(109,*) 3
      endif
      write(109,*) '// initial/stress-free temperature (in Kelvins)'
      write(109,*) temp_stress_free
      write(109,*) '// for cyclic loading, enter: Thermal_cycle_amplitude Thermal_cycle_period'
      write(109,*) '// for temp-time history. line1: number of temp-time points provided.'// &
                   'lines below, list of pairs: time, temperature'
      if(thermLoad.EQ.'C') then
         write(109,*) temp_cyclic_amplitude, temp_cyclic_period
      elseif(thermLoad.EQ.'H') then
         write(109,*) N_temp_history
         do iTemp = 1,N_temp_history
            write(109,*) temp_time_history(iTemp),temp_value_history(iTemp)
         enddo
         write(109,*)
      elseif(thermLoad.EQ.'N') then
         ! no thermal loading
         write(109,*) 1
         write(109,*) 0, temp_stress_free
      endif
      close(109)
      
      write(*,*) 'reading grain information...'
      
      ! read SIMMETRIX - to - Dream3D region ID mapping (from Kourosh)
      grainMapABQ_to_D3D = 0
      iGrains_ABQ = 0
      error=0
      open(unit=unitANSYS,file=trim(fileANSYS))
      do while(error.eq.0)
         read(unitANSYS,'(A)',iostat=error) strLine
         if (error==0.and.strLine(1:3)=='MAT') then
            read(strLine,*) strWord,grainID_D3D
            iGrains_ABQ=iGrains_ABQ+1
            grainMapABQ_to_D3D(iGrains_ABQ) = grainID_D3D
         endif
      enddo	
      if(iGrains_ABQ.NE.nGrains_ABQ) then
         write(*,*) 'Number of grains in ABQ and ANSYS files do not match:', nGrains_ABQ,iGrains_ABQ
         stop
      else
         write(*,*) '# of grains in the mesh:', nGrains_ABQ
      endif
      close(unitANSYS)            

      if(nGrains_ABQ>N_GRAINS_MAX) then
         write(*,*) 'increase N_GRAINS_MAX to ',nGrains_ABQ,', the number of grains in the ABQ file'
         stop
      endif

      ! read grain orientation angles 
      open(unit=unitDream3DCSV,file=fileDream3DCSV)
      open(301,file='debug_convABQ.txt')
      read(unitDream3DCSV,*) nGrains_D3D
      ! check consistency of the CSV file and the ABQ-D3D map:
      do iGrains_ABQ=1,nGrains_ABQ
         grainID_D3D = grainMapABQ_to_D3D(iGrains_ABQ)
         if (grainID_D3D == 0) then
            write(*,*) 'D3D grain ID=0 remains. Make sure no bad cells remain in the Simmetrix output.'
            stop
         elseif (grainID_D3D > nGrains_D3D) then
            write(*,*) 'D3D grain ID',grainID_D3D,' referenced in the ANSYS file cannot be found in the CSV file'
            stop
         endif
      enddo

      ! check the data fields in the CSV file. assumed ordering of the fields: Grain_ID,EquivalentDiameters,EulerAngles_0,EulerAngles_1,EulerAngles_2
      read(unitDream3DCSV,'(A)') fieldsCSV
      if (fieldsCSV(1:70) /= 'Grain_ID,EquivalentDiameters,EulerAngles_0,EulerAngles_1,EulerAngles_2' .and. &
          fieldsCSV(1:72) /= 'Feature_ID,EquivalentDiameters,EulerAngles_0,EulerAngles_1,EulerAngles_2') then
         write(*,*) 'Unexpected CSV file format. Fields expected in order:'
         write(*,*) 'Grain_ID,EquivalentDiameters,EulerAngles_0,EulerAngles_1,EulerAngles_2 (Dream3D 5)'
         write(*,*) 'Feature_ID,EquivalentDiameters,EulerAngles_0,EulerAngles_1,EulerAngles_2 (Dream3D 6)'
         stop
      endif
      
      ! read in the size and orientation of D3D grains
      nGrains_removedD3D = 0
      do iLine=1,nGrains_D3D
         read(unitDream3DCSV,*) grainID_D3D,grainSize_D3D(grainID_D3D),grainOri_D3D(1:3,grainID_D3D)
         !CALL getABQgrainID_from_D3DgrainID(grainID,grainID_D3D,grainMapABQ_to_D3D,N_GRAINS_MAX)
         !if(grainID==0) then
         !   nGrains_removedD3D = nGrains_removedD3D + 1
         !   write(301,*) 'Simmetrix-removed grain. D3D_ID:',grainID_D3D,'Diameter',grainSize_D3D
         !   cycle
         !endif
         if(grainID_D3D>N_GRAINS_MAX) then
            write(*,*) 'increase N_GRAINS_MAX to',grainID_D3D,', the number of grains in the D3D file'
            stop
         endif
      enddo
      !close(301)
      close(unitDream3DCSV)
      
      !write(*,*) 'number of grains removed by Simmetrix: ', nGrains_removedD3D

      ! assign sizes&texture to ANSYS/FEM grains from D3D grains
      ! this assignment method is compatible with periodic microstructures:
      !     for multiple ANSYS/FEM grains have the same D3D grainID
      do iGrain_ABQ = 1,nGrains_ABQ
         grainID_D3D = grainMapABQ_to_D3D(iGrain_ABQ)
         grainSize(iGrain_ABQ) = grainSize_D3D(grainID_D3D)
         grainOri(:,iGrain_ABQ) = grainOri_D3D(:,grainID_D3D)
      enddo
      
      !set up element orientations and phases
      allocate(elemPhase(NELX))
      elemPhase(:)=0
      allocate(elemOri(NELX,3))
      elemOri(:,:)=0.D0
      do iElem=1,NELX
         elemOri(iElem,:)=grainOri(:,elemGrain(iElem))
      enddo
      
      write(*,*) 'generating input files...'

!==========================NOT USED ANYMORE==============================
! ***** Generate statevariables.dat
!      initStateVars(:)=0
!      open(200,file='initStateVars.dat')
!      read(200,*) nInitStateVars
!      read(200,*) (initStateVars(iVar),iVar=1,nInitStateVars)
!      close(200)
!      open(105,file=trim(fileCPFEM)//'/statevariable.dat')
!      do iElem=1,NELX
!         write(105,*) iElem,((elemOri(iElem,iOri)),iOri=1,3),1, &
!                     (initStateVars(iVar),iVar=1,nInitStateVars),elemPhase(iElem)
!      enddo
!      close(105)
! ***** Generate texture.dat
!         open(106,file=trim(fileCPFEM)//'/texture.dat')
!      do iElem=1,NELX
!         write(106,*) (elemOri(iElem,iOri),iOri=1,3),elemPhase(iElem)
!      enddo
!      close(106)
!==========================NOT USED ANYMORE=============================

      write(*,*) 'Enter the percent ratio of the transformed-Beta grains:'
      write(*,*) 'e.g.: '
      write(*,*) '0   : no transformed-Beta grains'
      write(*,*) '30  : 30% of grains are transformed-Beta (a typical value for bimodal Ti6242)'
      read(*,*) transB_ratio
      
! ***** Generate grainPhases.inp - orientations for each grain
      open(105,file=trim(fileCPFEM)//'/grainPhases.inp')
      write(105,*) nGrains_ABQ
      do iGrain=1,nGrains_ABQ
         CALL RANDOM_NUMBER(transB_rand)
         if (transB_rand*100.d0 > transB_ratio) then
            write(105,*) iGrain, 1                    ! primary Alpha
         else
            write(105,*) iGrain, 2                    ! transformed Beta
         endif
      enddo
      close(105)
      

! ***** Generate grainTexture.inp - orientations for each grain
      open(105,file=trim(fileCPFEM)//'/grainTexture.inp')
      write(105,*) nGrains_ABQ
      do iGrain=1,nGrains_ABQ
         write(105,*) iGrain, grainOri(1:3,iGrain)
      enddo
      close(105)
      
! ***** Generate elementGrains.inp
      open(107,file=trim(fileCPFEM)//'/elementGrains.inp')
      write(107,*) NELX
      do iElem=1,NELX
         write(107,*) iElem, elemGrain(iElem)
      enddo
      close(107)
      
! ***** Generate grainSizes.inp
      open(108,file=trim(fileCPFEM)//'/grainSizes.inp')
      write(108,*) nGrains_ABQ
      do grainID=1,nGrains_ABQ
         write(108,*) grainID, grainSize(grainID)
      enddo
      close(108)
      
      write(*,*) 'done.'
      
      deallocate(elemFacesMarked)
      deallocate(elemFaceNeighbors)
      deallocate(boundaryNodes)
      deallocate(elemPhase)
      deallocate(elemOri)

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
      
      SUBROUTINE toUpperCase(string)
      implicit none
      character(len=*), intent(inout) :: string
      integer :: ascii, iPos
      do iPos=1,len(string)
         ascii = IACHAR(string(iPos:iPos))
         if(ascii.GE.97.AND.ascii.LE.122) string(iPos:iPos)=ACHAR(ascii-32)
      enddo
      END SUBROUTINE
      
      ! identifies the nodes on each face of the cubic domain
      SUBROUTINE identifyNodes(boundaryNodes,nBoundaryNodes,cornerNodes, &
                                domainRange,G0XYZ,elemFaceNeighbors,      &
                                IJK,NELX,NX,MAXNODESonEachFace)
                              
      use meshtools

      implicit none
      
      integer, intent(out) :: boundaryNodes(6,MAXNODESonEachFace)
      integer, intent(out) :: nBoundaryNodes(6)
      integer, intent(out) :: cornerNodes(8)
      real(8), intent(out) :: domainRange(3,2)
      real(8), intent(in)  :: G0XYZ(NX*3)
      integer, intent(in)  :: elemFaceNeighbors(NELX,4)
      integer, intent(in)  :: IJK(NELX*4)
      integer, intent(in)  :: NELX,NX,MAXNODESonEachFace

      integer :: iElem, iNode, iFace, nodeID, iLocNode, iDIM, I, J, coeff
      integer :: MDIM
      real(8) :: domainLen(3)
      logical :: nodeAlreadyAdded(NX,3)
      real(8) :: pos, tol(3), xNodes(3,4), xNode(3)
      real(8) :: nodeDiagonalPos(8),cornerNodeDiagonalPos(8), diagTransform(8,3)
      
      nodeAlreadyAdded(:,:)=.FALSE.
      nBoundaryNodes(:)=0
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
      
      do iElem=1,NELX
         do iFace=1,4
            if (elemFaceNeighbors(iElem,iFace).EQ.0) then
               ! this is a face at the boundary
               ! determine the direction it faces
               CALL getElemXNodes(iElem,IJK,G0XYZ,4,xNodes)
               CALL faceNormal(iFace,xNodes,iDIM)
               do iLocNode=1,4
                  if(iLocNode.EQ.iFace) cycle
                  ! this is a node on domain boundary that faces direction iDIM
                  ! add to the list of boundary nodes
                  nodeID = IJK((iElem-1)*4+iLocNode)
                  if (nodeAlreadyAdded(nodeID,iDIM)) cycle
                  pos = xNodes(iDIM,iLocNode)
                  if (pos+tol(iDIM).GT.domainRange(iDIM,2).AND.   &
                      pos-tol(iDIM).LT.domainRange(iDIM,2)) then
                     ! node is on the lower bound edge in direction iDIM
                     nBoundaryNodes(iDIM*2)=nBoundaryNodes(iDIM*2)+1
                     boundaryNodes(iDIM*2,nBoundaryNodes(iDIM*2)) = nodeID
                     nodeAlreadyAdded(nodeID,iDIM)=.TRUE.
                  elseif (pos+tol(iDIM).GT.domainRange(iDIM,1).AND.  &
                          pos-tol(iDIM).LT.domainRange(iDIM,1)) then
                     ! node is on the upper bound edge in direction iDIM
                     nBoundaryNodes(iDIM*2-1)=nBoundaryNodes(iDIM*2-1)+1     
                     boundaryNodes(iDIM*2-1,nBoundaryNodes(iDIM*2-1)) = nodeID
                     nodeAlreadyAdded(nodeID,iDIM)=.TRUE.
                  endif
               enddo
            endif
         enddo
      enddo
      
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
      
      END SUBROUTINE
      
      SUBROUTINE faceNormal(iFace, xElemNodes, dir)
      implicit none
      integer, intent(in) :: iFace
      real(8), intent(in) :: xElemNodes(3,4)
      integer, intent(out):: dir
      
      !Local
      integer :: edgeV(4,3)
      real(8):: edge1(3),edge2(3), norm(3)
      
      edgeV(1,1) = 2
      edgeV(1,2) = 3
      edgeV(1,3) = 4
      
      edgeV(2,1) = 1
      edgeV(2,2) = 3
      edgeV(2,3) = 4
      
      edgeV(3,1) = 1
      edgeV(3,2) = 2
      edgeV(3,3) = 4
      
      edgeV(4,1) = 1
      edgeV(4,2) = 2
      edgeV(4,3) = 3

      edge1(:) = xElemNodes(:,edgeV(iFace,3)) - xElemNodes(:,edgeV(iFace,1))
      edge2(:) = xElemNodes(:,edgeV(iFace,2)) - xElemNodes(:,edgeV(iFace,1))
      
      CALL cross(norm,edge1,edge2)
      
      norm(:)=abs(norm(:))
      dir = MAXLOC(norm(:),1)
      END SUBROUTINE
      
      SUBROUTINE cross(vc, v1, v2)
      implicit none
      real(8), intent(in) :: v1(3),v2(3)
      real(8), intent(out):: vc(3)
      
      vc(1)=+v1(2)*v2(3)-v1(3)*v2(2)
      vc(2)=-v1(1)*v2(3)+v1(3)*v2(1)
      vc(3)=+v1(1)*v2(2)-v1(2)*v2(1)
      
      END SUBROUTINE
      
      SUBROUTINE extractRegionID(strLine, regionID)
      character(len=*),intent(in):: strLine
      integer, intent(out) :: regionID
      
      integer :: iPos,ABQregionIDstartsAt
      parameter(ABQregionIDstartsAt=23)
      
      do iPos=ABQregionIDstartsAt,len(strLine)
         if(strLine(iPos:iPos).EQ.',') then
            read(strLine(ABQregionIDstartsAt:iPos-1),*) regionID
         endif
      enddo
      
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
      
      
!***********************************************************************
!C     WRITES LINEAR ARRAY OF ELEMENT NUMBERS FOR EACH NODE - deniz
!***********************************************************************
      SUBROUTINE writeNODE_ELEM_CONNEC(strConnectFile,IELEMNO,NINDX, &
                                       NELX,NX,MAXNODELEM,MAXNODE)
                                       
      implicit none
      character(len=*), intent(in) :: strConnectFile      

      integer, intent(in) :: NX,NELX,MAXNODELEM,MAXNODE
      integer, intent(in) :: IELEMNO(MAXNODELEM*MAXNODE),NINDX(MAXNODE)
!***********************************************************************
!C     IELEMNO: STORES ELEMENT NUMBERS
!C     NINDX: STORES START AND END INDEX IN IELEMNO ARRAY FOR
!C            EACH NODE
!***********************************************************************      
      !Local
      logical :: fileExists
      integer :: error, iElem, nElem, iNode, I, lastIndex
      integer :: elemStart, elemEnd

      
      open(716,file=trim(strConnectFile))
      write(716,*,iostat=error) NX
      lastIndex = 1
      do iNode=1,NX
         elemStart=1
         if(iNode.ne.1) elemStart=NINDX(iNode-1)
         elemEnd = NINDX(iNode)-1
         nElem = elemEnd - elemStart + 1
         
         write(716,*,iostat=error) nElem
         write(716,*,iostat=error) (IELEMNO(iElem),iElem=elemStart,elemEnd)
         lastIndex = lastIndex + nElem
      enddo
      
      close(716)
      
      END SUBROUTINE

            
      SUBROUTINE getABQgrainID_from_D3DgrainID(grainID_ABQ,grainID_D3D,grainMapABQ_to_D3D,N_GRAINS_MAX)
      implicit none
      integer, intent(out):: grainID_ABQ
      integer, intent(in) :: grainID_D3D
      integer, intent(in) :: grainMapABQ_to_D3D(N_GRAINS_MAX),N_GRAINS_MAX
      
      integer :: iGrainID_ABQ
      
      grainID_ABQ=0
      
      do iGrainID_ABQ=1,N_GRAINS_MAX
         if (grainMapABQ_to_D3D(iGrainID_ABQ)==grainID_D3D) then
            grainID_ABQ = iGrainID_ABQ
            return
         endif
      enddo
      END SUBROUTINE
      
      