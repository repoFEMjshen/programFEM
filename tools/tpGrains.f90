! development log
!----------------
!
! June 2015 - grain IDs, phases and Euler angles are plotted for inspecting consistency between D3D file and the input files (needs debug)

      program PostProcessing
      use meshtools
      use grains
      use crystal_Titanium
      implicit none
      include 'PARDIS.H'
      
      integer :: NX, NELX
      integer, parameter :: NODE = 4
      integer :: istep, nstep, nstepInp, nstepTP,nstep_max, lastProgress
      integer :: iTime
      real(8) :: tTime
      integer :: ENTRIESPERLINE      
      parameter(ENTRIESPERLINE = 500)

      real(8), allocatable :: nodalPositions(:,:)    ! stores the node positions (read from disp.out)
      real(8), allocatable :: AveStress(:,:)
      real(8) :: xNodes(3,4), elemJacMat(3,3), elemJac
      real(8) :: dispFactor
      real(8) :: G0XYZ(MAXCRD)
      integer :: IJK(MNELX), IJKNEL(MAXELN)
      INTEGER:: IELEMNO(MAXNODELEM*MAXNODE),NINDX(MAXNODE)
      real(8) :: xc(3)
      
   ! variables for domain geometry
      real(8) :: domainRange0(3,2), domainLen0(3)
      real(8) :: domainRange(3,2), domainLen(3)
      
      ! RVE-averaged quantities
      real(8), allocatable :: elemVol(:)
      real(8) :: totVol
      real(8), allocatable :: gaussRVE(:)

      ! element groupings
      integer, allocatable:: elemFaceNeighbors(:,:)

      ! grains
      character(len=50) :: grainFileName
      character(len=150) :: lineStr
      logical :: fileExists, successReadNeighbors
      integer :: iGrain
      integer :: lastIndex, iNodeIndex, iGloNode

      real(8) :: xset(4),yset(4),zset(4), xpos, ypos, zpos
      real(8) :: distance, weight, totWeight
      integer :: elemNodes(4)
      integer :: val_ind, outputTimeStepSkip
      integer :: timeStepForDeformationField
      character(len=30) :: meshName, fileName
      character :: GaussOrNodal, plotConfig, commentsInData, commentsInDISP, timeawareData, hasGrainInfo
      character :: plotGrainInfo, eulerDirection
      integer :: iElem,iFace,iNode,iNodeBlock,iElemBlock
      integer :: iEuler,iDOF,isys,iDir,nQuantDOF,nOutQuantDOF
      real(8) :: euler(3), rot(3,3), nFace(3)
      real(8) :: thisSchmidt
      real(8) :: gst(3,nslip_max),gmt(3,nslip_max)
      real(8), allocatable :: grainSchmidtBasalMax(:,:), grainSchmidtPrismMax(:,:)
      real(8), allocatable :: grainSchmidtBasalAll(:,:,:), grainSchmidtPrismAll(:,:,:)
      real(8), allocatable :: grainBasalPole(:,:)
      real(8), allocatable :: elementSlipTrace(:,:)
      real(8) :: rveSchmidtBasal(3), rveSchmidtPrism(3)
      integer, allocatable :: outQuantDOF(:)
      logical, allocatable :: export_QuantDOF(:)
      logical :: onlyPlotSchmidFactors
      logical :: errorL
      integer :: len,i,j,k,kk,ITEMP,error
      
!------------------------------------------------------------------------------      

      ! READ IN THE MESH, CONNECTIVITY
      CALL READGM(G0XYZ,IJK,NX,NELX,meshName)      

      write(*,*) 'Mesh imported:'
      write(*,*) '# of elements: ', NELX
      write(*,*) '# of nodes: ', NX
      
      ! construct the node-elem connectivity array
      CALL NODE_ELEM_CONNEC(IELEMNO,NINDX,IJK,NX,NELX,MNELX,MAXNODELEM)
      
      ! -------- read node displacements --------!
      allocate(nodalPositions(NX,3))
      do iNode=1,NX
         nodalPositions(iNode,1:3)=G0XYZ((iNode-1)*3+1:iNode*3)   !plot on reference configuration
      enddo

      dispFactor = 1.D0
      write(*,*) 'Displacement multiplier for the plot?'
      write(*,*) '(enter 0 to plot on undeformed shape)'
      read(*,*) dispFactor
      
      if (dispFactor /= 0.d0) then
         write(*,*) 'Enter the line number in time.out to indicate the deformed shape'
         read(*,*) timeStepForDeformationField
         open(101,file='disp.out')
         do istep = 1,timeStepForDeformationField
            read(101,*) iTime,tTime
            do iNode=1,NX
               read(101,*,iostat=error)(nodalPositions(iNode,iDOF), iDOF=1,3)           !only last step
               if (error /= 0) exit
               nodalPositions(iNode,1:3) = G0XYZ((iNode-1)*3+1:iNode*3) & !plot on deformed configuration
                                          +dispFactor*nodalPositions(iNode,1:3)
            enddo
            if (error /= 0) then
               write(*,*) 'displacement file cannot be opened, or is empty'
               stop
            endif   
         enddo
         close(101)
      endif

      write(*,*) 'Nodal positions imported'
      write(*,*)
      
      write(*,*) 'Plot only Schmid Factors (T), or full crystallographic information? (F)'
      read(*,*) onlyPlotSchmidFactors


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
                                   NELX,4,MAXNODELEM,MAXNODE)
                                   
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
      CALL identifyBoundaryNodes(G0XYZ,NELX,NX,IJK,elemFaceNeighbors,domainRange0)

      ! domainRange0: initial domain range
      ! domainRange: may evolve in time
      domainRange = domainRange0
      
      write(*,'(A,I0,A)') ' ', nTotGrainNodes,' # of grain-nodes'
            
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
      do iElem=1,NELX
         CALL getElemXNodes(iElem,IJK,G0XYZ,NODE,xNodes)
         CALL calcJacTET4(xNodes, elemJacMat, elemJac, elemVol(iElem))   
      enddo
      
      plotGrainInfo='Y'
      
      eulerDirection=' '
      do while (eulerDirection.NE.'A'.AND.eulerDirection.NE.'P')
         write(*,*) 'Grain orientations are provided in form of Euler angles'
         write(*,*) 'What is the direction of the active transformation associated with the Euler angles?'
         write(*,*) '(A) Euler angles rotate the specimen coordinate axes onto the crystal coordinate axes '&
                   &'(OIM standard, Kourosh, Deniz)'
         write(*,*) '(P) Euler angles rotate the crystal coordinate axes onto the specimen coordinate axes '&
                   &'(was used by Deka)'
         write(*,*) 'What is the direction of the active transformation associated with the Euler angles?'
         read(*,'(A1)') eulerDirection
         if (eulerDirection.EQ.'a') eulerDirection='A'
         if (eulerDirection.EQ.'p') eulerDirection='P'
      enddo
      write(*,*)
      
      ! calculate schmidt factors for each grain
      ! first initialize the crystal to get cst_HCP,cmt_HCP
      call init_crys_HCP()
      allocate(grainBasalPole(3,nGrains))
      allocate(elementSlipTrace(3,NELX))
      allocate(grainSchmidtBasalMax(3,nGrains))
      allocate(grainSchmidtPrismMax(3,nGrains))
      allocate(grainSchmidtBasalAll(3,3,nGrains))
      allocate(grainSchmidtPrismAll(3,6,nGrains))
      grainSchmidtBasalMax(:,:) = 0.d0
      grainSchmidtPrismMax(:,:) = 0.d0
      grainSchmidtBasalAll(:,:,:) = 0.d0
      grainSchmidtPrismAll(:,:,:) = 0.d0
      
      open(unit=505, file='grainSchmidtFactors.out')      
      write(505,*) 'grainID, max basal Schmidt factor along X,Y,Z directions, max prism Schmidt factor along X,Y,Z directions'
      do iGrain=1,nGrains

         euler(:) = grainTexture(:,iGrain)
         call euler_slip(euler,rot)
         if (eulerDirection=='P') then
            rot = TRANSPOSE(rot)
         endif
         
         do isys=1,nslip_HCP
            gst(:,isys) = MATMUL(rot,cst_HCP(:,isys))
            gmt(:,isys) = MATMUL(rot,cmt_HCP(:,isys))
         enddo
         
         grainBasalPole(:,iGrain) = gmt(:,1) ! basal pole
         if (grainBasalPole(3,iGrain) < 0.d0) then
            grainBasalPole(:,iGrain) = - grainBasalPole(:,iGrain)
         endif
         
         do isys=1,3   
            do iDir=1,3
               thisSchmidt = dabs(gst(iDir,isys)*gmt(iDir,isys))
               grainSchmidtBasalAll(iDir,isys,iGrain) = thisSchmidt
               if(thisSchmidt > grainSchmidtBasalMax(iDir,iGrain)) &
                  grainSchmidtBasalMax(iDir,iGrain) = thisSchmidt
            enddo
         enddo

         do isys=4,6
            do iDir=1,3
               thisSchmidt = dabs(gst(iDir,isys)*gmt(iDir,isys))
               grainSchmidtPrismAll(iDir,isys,iGrain) = thisSchmidt
               if(thisSchmidt > grainSchmidtPrismMax(iDir,iGrain)) &
                  grainSchmidtPrismMax(iDir,iGrain) = thisSchmidt
            enddo
         enddo
         
         
         
         write(505,*) iGrain,grainSchmidtBasalMax(:,iGrain),grainSchmidtPrismMax(:,iGrain)
      enddo
      
      ! calculate slip traces in each element face
      do iElem = 1, NELX
      
         iGrain = grainIDelem(iElem)
         
         elementSlipTrace(:,iElem) = 0.d0
         do iFace = 1,4
            if (elemFaceNeighbors(iElem,iFace) /= 0) cycle
            
            call getElemXNodes(iElem,IJK,G0XYZ,NODE,xNodes)
            call faceNormalVector(iFace, xNodes, nFace)
            call cross(elementSlipTrace(:,iElem),grainBasalPole(:,iGrain),nFace)
            
         enddo
         
      enddo
      
      write(505,*) ''
      write(505,*) 'list of Schmidt factors for individual slip directions:'
      do iGrain=1,nGrains
         write(505,*)
         write(505,*) iGrain
         do isys=1,3   
               write(505,*) 'basal, isys=',isys,' along X,Y,Z'
               write(505,*) grainSchmidtBasalAll(:,isys,iGrain)
         enddo
         write(505,*) iGrain
         do isys=4,6
               write(505,*) 'prism, isys=',isys,' along X,Y,Z'
               write(505,*) grainSchmidtPrismAll(:,isys,iGrain)
         enddo
      enddo
      
      close(505)
      
      
      ! OUTPUT RVE-averaged schmidt factors
      totVol = 0.d0
      rveSchmidtBasal(:) = 0.d0
      rveSchmidtPrism(:) = 0.d0
      do iElem=1,NELX
         iGrain = grainIDelem(iElem)
         rveSchmidtBasal(:) = rveSchmidtBasal(:) + elemVol(iElem)*grainSchmidtBasalMax(:,iGrain)
         rveSchmidtPrism(:) = rveSchmidtPrism(:) + elemVol(iElem)*grainSchmidtPrismMax(:,iGrain)
         totVol = totVol + elemVol(iElem)
      enddo
      rveSchmidtBasal = rveSchmidtBasal / totVol
      rveSchmidtPrism = rveSchmidtPrism / totVol
      open(504, file='rveSchmidtFactors.out')
      write(504,*) 'RVE-averaged: max basal Schmidt factor along X,Y,Z directions, max prism Schmidt factor along X,Y,Z directions'
      write(504,*) '   ',rveSchmidtBasal(:),rveSchmidtPrism(:)
      write(504,*) '   '
      write(504,*) 'total RVE volume:',totVol
      close(504)
      
      
      ! OUTPUT AS TECPLOT FILE
      ! ------------------------------!
      
      ! first plot phases / grain IDs !
      ! ------------------------------!
      if (plotGrainInfo=='Y'.and.grainsImported) then
         write(*,*) 'plotting grain phases and grain IDs'
         open(602, file='grains_TPG.plt')
         ! ---------------- HEADER ----------------- !
         ! ----------------------------------------- !
         write(602,'(A)',advance='no') (' VARIABLES= "X","Y","Z"')
         if (onlyPlotSchmidFactors) then
            write(602,'(A)') ',"grainID","schBasalX","schBasalY","schBasalZ",&
                            &"schPrismX","schPrismY","schPrismZ"'
            write(602,'(A,I0,A,I0,A)') &
               'ZONE NODES=',nTotGrainNodes,', ELEMENTS=',NELX,', ZONETYPE=FETETRAHEDRON, DATAPACKING=BLOCK, &
              &VARLOCATION=([4-10]=CELLCENTERED)'
         else  ! print only Schmid factors
            write(602,'(A)') ',"grainID","grainPhase","euler1","euler2","euler3","schBasalX","schBasalY","schBasalZ",&
                            &"schPrismX","schPrismY","schPrismZ","basalSlipTrace1","basalSlipTrace2","basalSlipTrace3",&
                            &"cAxisX","cAxisY","cAxisZ","elemID"'
            write(602,'(A,I0,A,I0,A)') &
               'ZONE NODES=',nTotGrainNodes,', ELEMENTS=',NELX,', ZONETYPE=FETETRAHEDRON, DATAPACKING=BLOCK, &
               VARLOCATION=([4-21]=CELLCENTERED)'
         endif
         
         ! there are grains. export with duplicate boundary nodes
         do i=1,3
            do iNodeBlock=1,INT(nTotGrainNodes/ENTRIESPERLINE)
               write(602,*) (nodalPositions(grainNodes(iNodeIndex),i),   &
                              iNodeIndex=(iNodeBlock-1)*ENTRIESPERLINE+1,iNodeBlock*ENTRIESPERLINE)
            enddo
            !if(iNodeBlock*ENTRIESPERLINE.NE.NX) then  ! print the remaining
               write(602,*) (nodalPositions(grainNodes(iNodeIndex),i),   &
                              iNodeIndex=INT(nTotGrainNodes/ENTRIESPERLINE)*ENTRIESPERLINE+1,nTotGrainNodes)
            !end if
         enddo
         
         ! grain IDs
         do iElemBlock=1,INT(NELX/ENTRIESPERLINE)
            write(602,*) (grainIDelem(iElem), &
                          iElem=(iElemBlock-1)*ENTRIESPERLINE+1,iElemBlock*ENTRIESPERLINE)
         enddo
         write(602,*) (grainIDelem(iElem),iElem=INT(NELX/ENTRIESPERLINE)*ENTRIESPERLINE+1,NELX)
         
         if (.not.onlyPlotSchmidFactors) then
         
            ! grain Phases
            do iElemBlock=1,INT(NELX/ENTRIESPERLINE)
               write(602,*) (grainPhase(grainIDelem(iElem)), &
                             iElem=(iElemBlock-1)*ENTRIESPERLINE+1,iElemBlock*ENTRIESPERLINE)
            enddo
            write(602,*) (grainPhase(grainIDelem(iElem)),iElem=INT(NELX/ENTRIESPERLINE)*ENTRIESPERLINE+1,NELX)
            
            ! euler angles
            do iEuler=1,3
               do iElemBlock=1,INT(NELX/ENTRIESPERLINE)
                  write(602,*) (grainTexture(iEuler,grainIDelem(iElem)), &
                                iElem=(iElemBlock-1)*ENTRIESPERLINE+1,iElemBlock*ENTRIESPERLINE)
               enddo
               write(602,*) (grainTexture(iEuler,grainIDelem(iElem)),iElem=INT(NELX/ENTRIESPERLINE)*ENTRIESPERLINE+1,NELX)
            enddo
         endif
         
         ! basal schmidt factors
         do iDir=1,3
            do iElemBlock=1,INT(NELX/ENTRIESPERLINE)
               write(602,*) (grainSchmidtBasalMax(iDir,grainIDelem(iElem)), &
                             iElem=(iElemBlock-1)*ENTRIESPERLINE+1,iElemBlock*ENTRIESPERLINE)
            enddo
            write(602,*) (grainSchmidtBasalMax(iDir,grainIDelem(iElem)),iElem=INT(NELX/ENTRIESPERLINE)*ENTRIESPERLINE+1,NELX)
         enddo
         
         ! prismatic schmidt factors
         do iDir=1,3
            do iElemBlock=1,INT(NELX/ENTRIESPERLINE)
               write(602,*) (grainSchmidtPrismMax(iDir,grainIDelem(iElem)), &
                             iElem=(iElemBlock-1)*ENTRIESPERLINE+1,iElemBlock*ENTRIESPERLINE)
            enddo
            write(602,*) (grainSchmidtPrismMax(iDir,grainIDelem(iElem)),iElem=INT(NELX/ENTRIESPERLINE)*ENTRIESPERLINE+1,NELX)
         enddo
         
         if (.not.onlyPlotSchmidFactors) then

            ! slip trace vectors on surface elements
            do iDir=1,3
               do iElemBlock=1,INT(NELX/ENTRIESPERLINE)
                  write(602,*) (elementSlipTrace(iDir,iElem), &
                                iElem=(iElemBlock-1)*ENTRIESPERLINE+1,iElemBlock*ENTRIESPERLINE)
               enddo
               write(602,*) (elementSlipTrace(iDir,iElem),iElem=INT(NELX/ENTRIESPERLINE)*ENTRIESPERLINE+1,NELX)
            enddo
            ! c-axis
            
            do iDir=1,3
               do iElemBlock=1,INT(NELX/ENTRIESPERLINE)
                  write(602,*) (grainBasalPole(iDir,grainIDelem(iElem)), &
                                iElem=(iElemBlock-1)*ENTRIESPERLINE+1,iElemBlock*ENTRIESPERLINE)
               enddo
               write(602,*) (grainBasalPole(iDir,grainIDelem(iElem)),iElem=INT(NELX/ENTRIESPERLINE)*ENTRIESPERLINE+1,NELX)
            enddo
            
            ! element IDs
            do iElemBlock=1,INT(NELX/ENTRIESPERLINE)
               write(602,*) (iElem, &
                             iElem=(iElemBlock-1)*ENTRIESPERLINE+1,iElemBlock*ENTRIESPERLINE)
            enddo
            write(602,*) (iElem,iElem=INT(NELX/ENTRIESPERLINE)*ENTRIESPERLINE+1,NELX)
         endif
         ! ---------------- CONNECTIVITY ----------------- !
         ! ----------------------------------------------- !
         ! -------- connectivity information is ---------- !
         ! ------------imported from zone 1 -------------- !
         ! ----------- see CONNECTIVITYSHAREZONE --------- !
         ! ----------------------------------------------- !
         ! --------------only written in first zone------- !
         ! has grains. export with duplicate boundary nodes
         do iElem=1,NELX
            write(602,*) IJK_grainNodes(:,iElem)
         enddo         
         close(602)
         
      endif


      
      CALL grains_Destruct()
      ! element groupings
      deallocate(elemFaceNeighbors)
      deallocate(grainSchmidtPrismMax,grainSchmidtBasalMax)
      deallocate(grainSchmidtPrismAll,grainSchmidtBasalAll)
      deallocate(grainBasalPole)

      deallocate(elemVol)
      
      deallocate(nodalPositions)
      
      end
        
      
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
      

      SUBROUTINE findDispIndex(this_N,disp_ind,disp_Ind_last,nodalPositions_N,nStep,error)
      implicit none
      integer, intent(in) :: this_N
      integer, intent(out) :: disp_ind
      integer, intent(inout) :: disp_Ind_last
      integer, intent(in) :: nodalPositions_N(nStep),nStep
      logical, intent(out) :: error

      error = .FALSE.
      
      disp_ind = disp_Ind_last
      do while(nodalPositions_N(disp_ind).NE.this_N)
         disp_ind = disp_ind + 1
         if(disp_ind.GT.nStep) then
            error=.TRUE.
            return
         endif
      enddo
      
      disp_Ind_last = disp_ind
      
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
      
      ! tlgt : euler angles describe an active rotation that rotates the specimen axes onto the crystal axes
      ! tlg  : euler angles describe an active rotation that rotates the crystal axes onto the specimen axes
      subroutine euler_slip(euler,tlgt)  
      implicit double precision (a-h,o-z)  

!      include 'aba_param.inc'

      dimension euler(3),tlg(3,3),tlgt(3,3)

      phi=euler(1)
      theta =euler(2)
      omega  =euler(3)

  
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

      SUBROUTINE cross(vc, v1, v2)
      implicit none
      real(8), intent(in) :: v1(3),v2(3)
      real(8), intent(out):: vc(3)
      
      vc(1)=+v1(2)*v2(3)-v1(3)*v2(2)
      vc(2)=-v1(1)*v2(3)+v1(3)*v2(1)
      vc(3)=+v1(1)*v2(2)-v1(2)*v2(1)
      
      END SUBROUTINE
      