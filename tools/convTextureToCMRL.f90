      PROGRAM convTextureToCMRL
      
      use meshtools
      use crystal_Titanium
      !use options
      use grains
      
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
      CALL READGM(G0XYZ,IJK,NX,NELX,meshName)    
      CALL READNSET(NSET,NNSET,NLIST)
      CALL READELSET(NESET,NELSET,NELLIST)
      CALL READSTR(LR,lineStr,error)
      CALL READSTR(LR,lineStr,error)   ! simulation time
      CALL SKIP_READMT()

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
         write(*,*) 'this mesh already has the newCMRL/ABAQUS convention for node numbering'
         write(*,*) 'connectivity and faceIDs will not be swapped'
      elseif (nElems_oldConvention == NELX) then
         write(*,*) 'this mesh has the oldCMRL convention for node numbering'
         write(*,*) 'connectivity and faceIDs will be swapped'
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
      
      INQUIRE(file='MS_nodes.inp',exist=fileExists)
      if(.NOT.fileExists) then

         ! create list of nodes
         write(*,*) 'creating MS_nodes.inp ... from the main input file'
         open(105,file='MS_nodes.inp')
         write(105,*) NX
         do iNode=1,NX
            nodePositions(:,iNode) = G0XYZ((iNode-1)*3+1:(iNode-1)*3+3)
            write(105,*) iNode, nodePositions(:,iNode)
         enddo
         close(105)
         write(*,*) 'created MS_nodes.inp'
      endif
      
      INQUIRE(file='MS_connectivity.inp',exist=fileExists)
      if(.NOT.fileExists) then

         ! create list of nodes
         write(*,*) 'creating MS_connectivity.inp ... from the main input file'
         open(105,file='MS_connectivity.inp')
         write(105,*) NELX, 'TET4F'
         do iElem=1,NELX
            write(105,*) iElem, IJK((iElem-1)*4+1:(iElem-1)*4+4)
         enddo
         close(105)
         write(*,*) 'created MS_connectivity.inp'
      endif
      
      CALL readGrainsFromElementList('elementGrains.inp',NX,NELX)
      if(.NOT.grainsImported) CALL createSingleGrain(NX,NELX)
      if (grainsImported) write(*,'(I0,A)') nGrains,' grains imported.'
      
      INQUIRE(file='MS_featureElements.inp',exist=fileExists)
      if(.NOT.fileExists) then
         
         ! create list of grain <- element set assignments
         write(*,*) 'creating MS_featureElements.inp ... from elementGrains.inp'
         open(105,file='MS_featureElements.inp')
         write(105,*) nGrains
         do iGrain=1,nGrains
            write(105,*) 'ElementSet',iGrain,iGrain
         enddo
         close(105)
         write(*,*) 'created MS_featureElements.inp'
         
         ! create element sets, listing elements in grains
         write(*,*) 'creating MS_elemSets.inp ...'
         open(105,file='MS_elemSets.inp')
         nElementSetsCMRL = nGrains + NESET + 1    ! nGrains + top face elements + all elements
         write(105,*) nElementSetsCMRL
         nTotalElements = 0
         do iGrain=1,nGrains
            staIndex = 1
            if (iGrain.NE.1) staIndex = grainElementIndx(iGrain-1)
            endIndex = grainElementIndx(iGrain)-1
            nElem = endIndex - staIndex + 1
            write(strGrain,'(A,I0)') 'grain',iGrain
            write(105,*) iGrain, trim(strGrain), '  Explicit'
            write(105,*) nElem, (grainElements(iIndex),iIndex=staIndex,endIndex)         
            nTotalElements = nTotalElements + nElem
         enddo
         do iElemSet=1,NESET
            write(105,*) nGrains+iElemSet, 'nameElemSet Explicit'
            write(105,*) NELSET(iElemSet), NELLIST(1:NELSET(iElemSet),iElemSet)
         enddo
         write(105,*) nGrains + NESET + 1, 'entireMesh Incremental'
         write(105,*) 1, nTotalElements, 1
         close(105)
         write(*,*) 'created MS_elemSets.inp'
         
         ! create element sets, listing elements in grains
         write(*,*) 'creating MS_nodeSets.inp'
         open(105,file='MS_nodeSets.inp')
         write(105,*) NSET
         do iNodeSet = 1,NSET
            write(105,*) iNodeSet,'nameNodeSet Explicit'
            write(105,*) NNSET(iNodeSet), NLIST(1:NNSET(iNodeSet),iNodeSet)
         enddo
         close(105)
         
      endif
      
      INQUIRE(file='MS_featureProps.inp',exist=fileExists)
      if(.NOT.fileExists) then
      
         write(*,*) 'creating MS_featureProps.inp ... from grainTexture.inp, grainSize.inp, grainPhases.inp'
            
         write(*,*) 'Convert between active/passive conventions? (transpose the rotations). T/F' 
         read(*,*) convertConvention
         
         INQUIRE(file='grainTexture.inp',exist=fileExists)
         if (.not. fileExists) then
            write(*,*) 'grainTexture.inp not found'
            stop
         endif
         CALL readTexture('grainTexture.inp')
         if (convertConvention) then
            do iGrain=1,nGrains
               euler(:) = grainTexture(:,iGrain)
               grainTexture(1,iGrain) = - euler(3)
               grainTexture(2,iGrain) = - euler(2)
               grainTexture(3,iGrain) = - euler(1)
            enddo
         endif

         ! import grain sizes
         CALL readGrainSizes('grainSizes.inp')
         if(.not.grainSizesImported) then
            write(*,*) 'grainSizes.inp not found'
            stop
         endif
         
         ! import grain phases
         CALL readGrainPhases('grainPhases.inp')
         if(.not.grainPhasesImported) then
            write(*,*) 'grainPhases.inp not found -- assuming phase = 1'
         endif
         
         open(210,file='MS_featureProps.inp')
         write(210,*) nGrains, 5
         do iGrain=1,nGrains
            write(210,*) iGrain, 5, grainPhase(iGrain), grainTexture(1:3,iGrain), grainSize(iGrain)
         enddo
         close(210)
         
         write(*,*) 'created MS_featureProps.inp'

      endif
      
      ! create load case folder
      CALL execute_command_line('mkdir loadCase')
      
      ! read and write displacement and pressure boundary conditions
      write(*,*) 'load axis? 1=X, 2=Y, 3=Z'
      read(*,*) loadAxis
      nShiftElemSetIDs = nGrains
      nShiftNodeSetIDs = 0
      sampleLength = abs(maxval(nodePositions(loadAxis,1:NX))-minval(nodePositions(loadAxis,1:NX)))
      call READ_WRITE_BC('loadCase/BC_disp.inp','loadCase/BC_pressure.inp', &
                         swapIDs,nShiftElemSetIDs,nShiftNodeSetIDs,sampleLength,loadAxis)

      write(*,*) 'created BC_disp.inp, BC_pressure.inp'
      
      ! copy the PATCH file
      CALL execute_command_line('cp PATCH.inp loadCase/')
      
      ! create no-PATCH file
      open(210,file='loadCase/noPATCH.inp')
      write(210,*) NELX,1,1
      do i=1,NELX
         write(210,*) 1, i
      enddo
      close(210)
      
      ! create a default options file
      open(210,file='loadCase/options.inp')
      write(210,*) 'typeInitialGuessForDisplacement StiffnessBased'
      write(210,*) 'nCommunicationSplits 10'
      write(210,*) 'considerGNDs .false.'
      write(210,*) 'shouldIgnoreDisplacementsNegligible .true.'
      write(210,*) 'shouldIgnoreForcesNegligible .true.'
      write(210,*) 'toleranceRelativeStress 1.0e-5'
      write(210,*) 'toleranceRelativeStressMinimum 1.0e-10'
      write(210,*) 'frequencyOutputNodal 4'
      write(210,*) 'frequencyOutputElemental 4'
      close(210)
      
      ! finally create inputMain.inp
      open(210,file='loadCase/inputMain.inp')
      write(210,*) '// # dimensions, # DOFs'
      write(210,*) '3 3'
      write(210,*) '// nodes'
      write(210,*) 'File=../MS_nodes.inp'
      write(210,*) '// connectivity'
      write(210,*) 'File=../MS_connectivity.inp'
      write(210,*) '// node sets'
      write(210,*) 'File=../MS_nodeSets.inp'
      write(210,*) '// element sets'
      write(210,*) 'File=../MS_elemSets.inp'
      write(210,*) '// MATERIAL PROPERTIES'
      write(210,*) 'File=../material_Ti7_FBar.inp'
      write(210,*) '// Outputting Keywords'
      write(210,*) '1'
      write(210,*) '1 2 stressCauchy deformationGradientPlastic'
      write(210,*) '// Grain (feature) properties: # props, Euler angles, grain size, phase id, ...'
      write(210,*) 'File=../MS_featureProps.inp'
      write(210,*) '//Assign Feature(grain) IDs to element(sets)'
      write(210,*) 'File=../MS_featureElements.inp'
      write(210,*) '//Loading Details'
      write(210,*) 'File=XXX.inp'
      write(210,*) '//DISP BOUND CONDITIONS'
      write(210,*) 'File=BC_disp.inp'
      write(210,*) '//PRESSURE BOUNDARY CONDITIONS'
      write(210,*) 'File=BC_pressure.inp'
      write(210,*) '//Thermal loading'
      write(210,*) '1'
      write(210,*) '293.00   '
      write(210,*) '//PERIODIC BOUNDARY CONDITIONS'
      write(210,*) '0'

      write(*,*) 'all input files converted'
      
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

      return
      end
      
      SUBROUTINE READNSET(NSET,NNSET,NLIST)
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'PARDIS.H'
      INTEGER, intent(out):: NNSET(MAXNSET),NLIST(MAXNS,MAXNSET)
      INTEGER:: N
      integer :: error
      character(len=150) lineStr

      call READSTR(LR,lineStr,error)
      READ(LR,*)NSET

      IF(NSET.GT.MAXNSET)THEN
         WRITE(*,*)'INSUFFICIENT MEMORY-MAXNSET'
         STOP
      ENDIF

      DO I=1,NSET
         READ(LR,*)NSNO,NGEN
         IF(NGEN.EQ.0)THEN
            READ(LR,*)NUMN,(NLIST(J,NSNO),J=1,NUMN)
         ENDIF
         ! error check added - deniz
         IF(NUMN.GT.MAXNS) THEN
            WRITE(*,*)'Increase MAXNS to ', NUMN
            STOP
         END IF
         IF(NSET.GT.MAXNSET) THEN
            WRITE(*,*)'Increase MAXNSET to ', NSET
            STOP
         END IF
         IF(NGEN.EQ.1)THEN
            READ(LR,*)ISTART,IEND,INC
            NUMN=0
            DO N=ISTART,IEND,INC
               NUMN=NUMN+1
               NLIST(NUMN,NSNO)=N
            ENDDO
         ENDIF
         NNSET(NSNO)=NUMN
      ENDDO
      RETURN
      END SUBROUTINE

!*******************************************************
      SUBROUTINE READELSET(NESET,NELSET,NELLIST)
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'PARDIS.H'
      integer, intent(out) :: NESET,NELSET(MAXELSET),NELLIST(MAXELS,MAXELSET)
      integer :: error
      character(len=150) lineStr
      call READSTR(LR,lineStr,error)
      READ(LR,*)NESET

      IF(NESET.GT.MAXELSET)THEN
         WRITE(*,*)'INSUFFICIENT MEMORY-MAXELSET'
         WRITE(*,*)'INSUFFICIENT MEMORY-MAXELSET'
         STOP
      ENDIF

      DO I=1,NESET
         READ(LR,*)NSNO,NGEN
         IF(NGEN.EQ.0)THEN
            READ(LR,*)NUMN,(NELLIST(J,NSNO),J=1,NUMN)
         ENDIF
         if(NUMN.GT.MAXELS) then
            WRITE(*,*)'INSUFFICIENT MEMORY-MAXELS'
            WRITE(*,*)'INSUFFICIENT MEMORY-MAXELS'
            STOP
         endif
         IF(NGEN.EQ.1)THEN
            READ(LR,*)ISTART,IEND,INC
            NUMN=0
            DO N=ISTART,IEND,INC
               NUMN=NUMN+1
               NELLIST(NUMN,NSNO)=N
            ENDDO
         ENDIF
         NELSET(NSNO)=NUMN
      ENDDO
      
      END SUBROUTINE
      
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

      SUBROUTINE SKIP_READMT()
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'PARDIS.H'
      real(8) :: PROMAT2D(MAXPROPS,MAXGRP)
      integer :: NPROPSG(MAXGRP)
      character(len=150) :: lineStr
      integer :: error
      
      CALL READSTR(LR,lineStr,error)
      READ(LR,*)NGROUP
      IF(NGROUP.GT.MAXGRP)THEN
         WRITE(*,*)'INSUFFICIENT MEM -MAXGRP'
         STOP
      ENDIF

      READ(LR,*)(NPROPSG(N),N=1,NGROUP)
      DO I=1,NGROUP
         IF(NPROPSG(I).GT.MAXPROPS)THEN
            WRITE(*,*)'INSUFFICIENT MEM - MAXPROPS'
            STOP
         ENDIF
      ENDDO

      DO N=1,NGROUP
         MATGRP=NPROPSG(N)
         NCOUNT=MATGRP/8
         NREMAIN=MOD(MATGRP,8)
         NSHIFT=0
            WRITE(*,*)'NUM OF PROPS: ',MATGRP
         DO J=1,NCOUNT
            READ(LR,*)(PROMAT2D(NSHIFT+I,N),I=1,8)
            NSHIFT=NSHIFT+8
         ENDDO
         READ(LR,*)(PROMAT2D(NSHIFT+I,N),I=1,NREMAIN)
      ENDDO
      READ(LR,*)NSTATE
      IF(NSTATE.GT.MSV)THEN
         WRITE(*,*)'INSUFFICIENT MEMORY -MSV'
         STOP
      ENDIF
      RETURN
      END SUBROUTINE
      

      SUBROUTINE READ_WRITE_BC(fileNameBCdisp,fileNameBCpress,swapIDs, &
                               nShiftElemSetIDs,nShiftNodeSetIDs,sampleLength,loadAxis)
      IMPLICIT NONE

      INCLUDE 'PARDIS.H'
      
      integer, intent(in) :: nShiftElemSetIDs,nShiftNodeSetIDs
      real(8), intent(in) :: sampleLength
      integer, intent(in) :: loadAxis
      logical, intent(in) :: swapIDs
      character(len=*), intent(in) :: fileNameBCdisp,fileNameBCpress
      
      real(8) :: loadAxisVector(3)

      REAL(8):: dispBC(MAXDIM)
      INTEGER:: NBC3(MBOUND3),IK(MAXDIM)
      INTEGER:: NBC1(MBOUND1),NBC2(MBOUND2)
      INTEGER:: NNSET(MAXNSET),NLIST(MAXNS,MAXNSET)
      INTEGER:: IDF(MBOUND1), IBC(MDOFX),IBCT(MDOFX)
      INTEGER:: IPRESS_SET(MAXPRESS,MBOUND3),IPRESS_SF(MAXPRESS,MBOUND3)
      INTEGER:: ILF(MAXEL),IFACE_TET(12),IFACE_BRICK(24)
      INTEGER:: IBELEM(MAXEL)
      INTEGER:: NDF,I,IFUNC,IEND,II
      INTEGER:: IPF
      INTEGER:: ISTART,ISET,NCONSTR
      INTEGER:: N,NBN
      INTEGER:: NBOUD1,NBOUD2,NBOUD3
      INTEGER:: NCONST, NPLANE, NODE1,NN,NESET
      INTEGER:: NDF1,NEQ,NELX,NODE,NDFF,NGAUSS
      INTEGER:: NNOD
      INTEGER:: NPRESS(MBOUND3), iPress
      
      integer :: nTotalPressBCs
      
      character(len=50) :: strFunc
      character(len=150) :: lineStr
      integer :: error
      
      CALL READSTR(LR,lineStr,error)

      !    imposed displacements
      !----------------------------
      open(250,file=trim(fileNameBCdisp))
      READ(LR,*) NBOUD1  ! # of entries to read
      write(250,*) NBOUD1  ! # of entries to read
      DO N=1,NBOUD1
         READ(LR,*) NBN,ISET,IFUNC,ISTART,IEND
         if (IFUNC==1) then     ! program should pull the nodes with internal function
            strFunc = 'ConstantStrainRate'
         elseif (IFUNC==0) then ! displacements given as inputs
            strFunc = 'Explicit'
         else
            write(*,*) 'unknown IFUNC controller in input file'
            stop
         endif
         if (ISET==1) then ! nodeSet
            write(250,*) 'NodeSet ',NBN+nShiftNodeSetIDs,trim(strFunc),ISTART,IEND
         else              ! node
            write(250,*) 'Node ',NBN+nShiftNodeSetIDs,trim(strFunc),ISTART,IEND
         endif
         if (IFUNC==1) then               ! program should pull the nodes with internal function
            write(250,*) sampleLength
         elseif(IFUNC==0) THEN            ! displacements given as inputs
            dispBC(:)=0.D0
            read(LR,*)(dispBC(I),I=ISTART,IEND)
            write(250,*) (dispBC(I),I=ISTART,IEND)  
         ENDIF
      ENDDO
      close(250)
      

!     SKIP THE 2-ND TYPE BOUNDARY CONDITION : POINT FORCES
      CALL READSTR(LR,lineStr,error)
      READ(LR,*) NBOUD2
      DO N=1,NBOUD2
         READ(LR,*) NBN,ISET,IFUNC,ISTART,IEND
         IF(IFUNC.EQ.0)THEN
            CALL READSTR(LR,lineStr,error)
         ENDIF
      ENDDO
      
!      READ THE THIRD TYPE BOUNDARY CONDITION (PRESSURE)
      loadAxisVector = 0.d0
      loadAxisVector(loadAxis) = 1.d0
      open(250,file=trim(fileNameBCpress))
      CALL READSTR(LR,lineStr,error)
      READ(LR,*) NBOUD3
      if (NBOUD3 > MBOUND3) then
         write(*,*) 'increase MBOUND3-max # of pressure BCs'
         stop
      endif
      nTotalPressBCs = 0
      DO N=1,NBOUD3
         READ(LR,*)NPRESS(N),IPF,(IPRESS_SET(I,N),I=1,NPRESS(N)),(IPRESS_SF(I,N),I=1,NPRESS(N))
         if (swapIDs) then
            do I=1,NPRESS(N)
               if(IPRESS_SF(I,N)==1) then ! CMRL/ABAQUS convention has the nodes 3/4 swapped wrt oldCMRL convention.
                  IPRESS_SF(I,N)=2        ! this requires swapping faceIDs 1/2 for consistency
               elseif(IPRESS_SF(I,N)==2)  then
                  IPRESS_SF(I,N)=1
               endif
            enddo
         endif
         if (IPF /= 1) then
            WRITE(*,*)'unsupported option IPF=',IPF,' for the pressure BC'
            STOP
         endif
         nTotalPressBCs = nTotalPressBCs + NPRESS(N)
         if (NPRESS(N).GT.MAXPRESS)THEN
            WRITE(*,*)'Increase MAXNSET to ', NPRESS(N)
            STOP
         ENDIF
      ENDDO
      write(250,*) nTotalPressBCs
      do N = 1,NBOUD3
         do iPress=1,NPRESS(N)
            !write(250,*) 'ElementSet',iPress_SET(iPress,N)+nShiftElemSetIDs,iPress_SF(iPress,N),'NormalToFace'
            write(250,*) 'ElementSet',iPress_SET(iPress,N)+nShiftElemSetIDs,iPress_SF(iPress,N),'Explicit'
            write(250,*) loadAxisVector(1:3)
         enddo
      enddo
      close(250)

      CALL READSTR(LR,lineStr,error)
!     READ CONSTRAINTS
      READ(LR,*)NCONST
      NCONSTR=0
      
      ! finally, close LR
      close(LR)
      
!-----------------------------------
!      ! face 1 (CODE) ==== face 4 (deniz)
!      IFACE_TET(1)=1
!      IFACE_TET(2)=2
!      IFACE_TET(3)=3

!      ! face 2 (CODE) ==== face 3 (deniz)
!      IFACE_TET(4)=1
!      IFACE_TET(5)=2
!      IFACE_TET(6)=4

!      ! face 3 (CODE) ==== face 1 (deniz)
!      IFACE_TET(7)=2
!      IFACE_TET(8)=3
!      IFACE_TET(9)=4

      ! face 4 (CODE) ==== face 2 (deniz)
!      IFACE_TET(10)=1
!      IFACE_TET(11)=3
!      IFACE_TET(12)=4
!------------------------------------

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