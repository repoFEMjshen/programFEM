! General purpose converter: FE solution field outputs --> TecPlot files
! supports:
!  element/node based fields
!  solution time stamps
!
! author: Deniz Ozturk, 2014
!
      program PostProcessing
      use meshtools
      use grains
      implicit none
      include 'PARDIS.H'
      
      integer :: NX, NELX, NODE
      integer :: nstepTime, nstepInp, nstepTP,nstep_max, lastProgress, beginStepTP,nstepsSpecific
      integer, allocatable :: plotTimeSteps(:)
      real(8) :: plotLastSeconds
      integer :: ENTRIESPERLINE      
      parameter(ENTRIESPERLINE = 500)
      
      integer :: valuesAtNodes
      integer, parameter :: valuesUnique = 1    ! enum for valuesAtNodes
      integer, parameter :: valuesPerGrain = 2

      real(8), allocatable :: nodalPositions(:,:)    ! stores the node positions (read from disp.out)
      real(8) :: nodalPositions_t      ! time of the imported positions
      integer :: nodalPositions_N      ! step number of the imported positions
      real(8), allocatable :: gaussValues(:,:)       ! stores the gauss values of the solution
      real(8), allocatable :: nodalValues(:,:)       ! stores the nodal values of the solution
      real(8) :: times(100000), plotBeginTime
      real(8) :: values_t              ! time of the imported solution
      integer :: values_N              ! step number of the imported solution
      integer :: connZone
      logical :: connExported, found
      real(8), allocatable :: AveStress(:,:)
      real(8) :: xNodes(3,4), elemJacMat(3,3), elemJac
      real(8) :: dispFactor
      real(8) :: G0XYZ(MAXCRD)
      integer :: IJK(MNELX), IJKNEL(MAXELN)
      INTEGER:: IELEMNO(MAXNODELEM*MAXNODE),NINDX(MAXNODE)
      real(8) :: xc(3)
      
   ! variables for domain geometry
      real(8) :: domainRange(3,2), domainLen(3)
      
      ! RVE-averaged quantities
      real(8), allocatable :: elemVol(:)
      real(8) :: totVol

      ! element groupings
      integer, allocatable:: elemFaceNeighbors(:,:)

      ! grains
      character(len=50) :: grainFileName
      character(len=150) :: lineStr
      logical :: fileExists, successReadNeighbors,successReadNodeConnect
      logical :: foundDispBlock, plotThisTime
      integer :: iGrain
      integer :: lastIndex, iNodeIndex, iGloNode
      integer :: iGrainNode, nNodes_file

      real(8) :: xset(4),yset(4),zset(4), xpos, ypos, zpos
      real(8) :: distance, weight, totWeight
      integer :: elemNodes(4)
      integer :: val_ind, outputTimeStepSkip
      character(len=30) :: meshName, fileName, solnFileName
      character :: GaussOrNodal, plotConfig, commentsInData, commentsInDISP, timeawareData, hasGrainInfo
      character :: plotGrainInfo, rangeOrSpecify
      integer :: iElem,iNode,iNodeBlock,iElemBlock
      integer :: iEuler,iDOF,nQuantDOF,nOutQuantDOF
      integer, allocatable :: outQuantDOF(:)
      logical, allocatable :: export_QuantDOF(:)
      logical :: errorL
      integer :: len,i,j,k,kk,ITEMP,error
      integer :: step_N, this_N, disp_ind
      
!------------------------------------------------------------------------------      

      ! READ IN THE MESH, CONNECTIVITY
      CALL READGM(G0XYZ,IJK,NX,NELX,meshName)      

      write(*,*) 'Mesh imported:'
      write(*,*) '# of elements: ', NELX
      write(*,*) '# of nodes: ', NX
      
      NODE = 4
      
      write(*,'(A)',advance='no') ' identifying node-to-element connectivity... '
      CALL readNODE_ELEM_CONNEC('nodeElementConnect.dat',IELEMNO,NINDX,NX,successReadNodeConnect)
      if (.not.successReadNodeConnect) then
         write(*,'(A)',advance='no') 'calculating... '
         CALL NODE_ELEM_CONNEC(IELEMNO,NINDX,IJK,NX,NELX,MNELX,MAXNODELEM)
         write(*,*) 'done.'
      else
         write(*,*) 'imported.'
      endif
      
      ! READ TIME STEPPING INFORMATION
      open(102,file='time.out')
      nstepTime=0
      times(:) = 0.d0
      do
         read(102,*,iostat=error) times(nstepTime+1)
         if(error.EQ.-1) EXIT
         nstepTime=nstepTime+1
      enddo
      close(102)
      if (nstepTime.EQ.0) then
         write(*,*) 'Cannot find file: time.out'
         STOP
      endif
      
      write(*,*) '# of time steps:', nstepTime

   !------------------------------------------------------------------------------

      write(*,*) 'Specify the file that contains the field data:'
      read(*,*) solnFileName
      write(*,*)

      write(*,*) 'Enter # of components of the field variable.'
      write(*,*) '(Temperature: 1, Cauchy stress: 6, Fp: 9, etc.)'
      read(*,*) nQuantDOF
      write(*,*)
      
      allocate(outQuantDOF(nQuantDOF))
      allocate(export_QuantDOF(nQuantDOF))
      
      write(*,*) 'specify the components of the field variable you''d like to export:'
      write(*,*) 'example: F F T F T F --> 3rd and 5th components. T: export, F: do not export'
      read(*,*) export_QuantDOF(1:nQuantDOF)
      write(*,*)
      nOutQuantDOF=0
      do iDOF=1,nQuantDOF
         if(export_QuantDOF(iDOF))then
            nOutQuantDOF=nOutQuantDOF+1
            outQuantDOF(nOutQuantDOF) = iDOF
         endif
      enddo
      
      plotConfig=' '
      do while (plotConfig.NE.'D'.AND.plotConfig.NE.'R')
         write(*,*) 'Plot on (R)eference configuration or (D)eformed configuration?'
         read(*,'(A1)') plotConfig
         if (plotConfig.EQ.'d') plotConfig='D'
         if (plotConfig.EQ.'r') plotConfig='R'
      enddo
      write(*,*)
      
      dispFactor = 1.D0
      if(plotConfig.EQ.'D') then
         write(*,*) 'Displacement multiplier for the plot:'
         read(*,*) dispFactor
      endif
      
      !-------------------------------------
      
      GaussOrNodal=' '
      do while (GaussOrNodal.NE.'G'.AND.GaussOrNodal.NE.'N')
         write(*,*) 'Does the file contain gauss pt. values (G), or nodal values (N)?'
         read(*,'(A1)') GaussOrNodal
         if (GaussOrNodal.EQ.'g') GaussOrNodal='G'
         if (GaussOrNodal.EQ.'n') GaussOrNodal='N'
      enddo
      write(*,*)
            
      commentsInData=' '
      do while (commentsInData.NE.'Y'.AND.commentsInData.NE.'N')
         write(*,*) 'Are there comment/info lines in at the beginning of each time step in data files? (Y/N)'
         read(*,'(A1)') commentsInData
         if (commentsInData.EQ.'y') commentsInData='Y'
         if (commentsInData.EQ.'n') commentsInData='N'
      enddo
      write(*,*)

      commentsInDISP=' '
      do while (commentsInDISP.NE.'Y'.AND.commentsInDISP.NE.'N')
         write(*,*) 'Are there comment/info lines at the beginning of each time step in DISP.OUT? (Y/N)'
         read(*,'(A1)') commentsInDISP
         if (commentsInDISP.EQ.'y') commentsInDISP='Y'
         if (commentsInDISP.EQ.'n') commentsInDISP='N'
      enddo
      write(*,*)

      timeawareData='Y'
      do while (timeawareData.NE.'Y'.AND.timeawareData.NE.'N')
         write(*,*) 'include time information in techPlot output?'
         write(*,*) 'Requirements:'
         write(*,*) '  1) solution field has time information (SOLUTIONTIME field in ZONE header)'
         write(*,*) '  2) disp.out contains the displacements at the corresponding times'
         read(*,'(A1)') timeawareData
         if (timeawareData.EQ.'y') timeawareData='Y'
         if (timeawareData.EQ.'n') timeawareData='N'
      enddo
      write(*,*)

      rangeOrSpecify=' '
      do while (rangeOrSpecify.NE.'R'.AND.rangeOrSpecify.NE.'S')
         write(*,*) 'Export a range of time steps, or specify time steps? R / S'
         read(*,*) rangeOrSpecify
         write(*,*)
         if (rangeOrSpecify=='r') rangeOrSpecify='R'
         if (rangeOrSpecify=='s') rangeOrSpecify='S'
      enddo

      if(rangeOrSpecify=='R') then
         write(*,*) 'Enter max # of time steps to limit the export from the solution file'
         write(*,*) 'Enter 0 to export all time steps'
         write(*,*) 'Enter (-)(#) to export last (#) number of seconds'
         read(*,*) nstep_max
         write(*,*)
      else
         write(*,*) 'Enter the total # of time steps you would like to export:'
         read(*,*) nstepsSpecific
         write(*,*)
         allocate(plotTimeSteps(nstepsSpecific))
         write(*,*) 'Enter the line numbers in time.out that you would like to export:'
         read(*,*) plotTimeSteps(1:nstepsSpecific)
         write(*,*)
         nstep_max = 0 
      endif
      
      if (nstep_max==0) then
         beginStepTP = 1
         nstepTP = nstepTime ! # of time steps to export to TP file.
      elseif(nstep_max < 0) then
         ! find that time
         plotBeginTime = times(nstepTime) + nstep_max
         do i=1,nstepTime
            if (times(i) >= plotBeginTime) then
               beginStepTP = i
               exit
            endif
         enddo
         nstepTP = nstepTime
      else
         beginStepTP = 1
         nstepTP = nstep_max
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
      endif

      ! IDENTIFY THE NODES IN INDIVIDUAL GRAINS
      CALL identifyGrainNodes(NX,NELX,IJK)

      ! IDENTIFY THE NODES AT 
      ! 1) DOMAIN BOUNDARY, DOMAIN CORNERS
      ! 2) GRAIN BOUNDARIES
      ! 3) GRAIN-GRAIN BOUNDARIES
      CALL identifyBoundaryNodes(G0XYZ,NELX,NX,IJK,elemFaceNeighbors,domainRange)
      domainLen(:)=domainRange(:,2)-domainRange(:,1)

      
      write(*,*) 'Domain range detected - X', domainRange(1,:)
      write(*,*) 'Domain range detected - Y', domainRange(2,:)
      write(*,*) 'Domain range detected - Z', domainRange(3,:)
      
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
      
      plotGrainInfo='N'
      if(.not.grainsImported) plotGrainInfo='N'
      do while (plotGrainInfo.NE.'Y'.AND.plotGrainInfo.NE.'N')
         write(*,*) 'Plot grain phases, grain IDs and Euler angles too? (Y/N)'
         read(*,'(A1)') plotGrainInfo
         if (plotGrainInfo.EQ.'y') plotGrainInfo='Y'
         if (plotGrainInfo.EQ.'n') plotGrainInfo='N'
      enddo
      write(*,*)

      ! calculate the nodal values, if not provided
      ! removed in this version.
      ! if nodal values are given, field is plotted through the nodal points
      ! if gauss values are given, field is plotted through the gauss points (TECPLOT: CELL-centered variables, Block data packing)
   
   
      ! allocate buffers, open files for in/out
   
      allocate(nodalPositions(3,NX))
      if(plotConfig.EQ.'D') then
         inquire(file='disp.out',exist=fileExists)
         if(.NOT.fileExists) then
            write(*,*) 'disp.out not found (required for plotting the deformed shape)'
            stop
         endif
         open(101,file='disp.out')
      endif
      
      open(103,file=trim(solnFileName)//'.out')
      
      if(GaussOrNodal.EQ.'G') then 
         allocate(gaussValues(nQuantDOF,NELX))
         open(603, file=trim(solnFileName)//'_TPG.plt')
      endif
      if(GaussOrNodal.EQ.'N') then
         read(103,*) nNodes_file
         if (nNodes_file == NX) then
            valuesAtNodes = valuesUnique
         elseif (nNodes_file == nTotGrainNodes) then
            valuesAtNodes = valuesPerGrain
         else
            write(*,*) 'data file inconsistent with the mesh. total # of nodes in file:'
            write(*,*) 'mesh (unique nodes):',NX
            write(*,*) 'mesh (grain-nodes):',nTotGrainNodes
            write(*,*) 'file:',nNodes_file
            goto 700
         endif
         
         allocate(nodalValues(nQuantDOF,nNodes_file))
         open(603, file=trim(solnFileName)//'_TPN.plt')
      endif
      
      ! OUTPUT AS TECPLOT FILE
      ! ------------------------------!
      
      ! now plot the field variables  !
      ! ------------------------------!
      write(*,*) 'Processing...'
      lastProgress=0
      do i=1,MIN(60,nstepTP)
         write(*,'(A)',advance='no') '_'
      enddo
      write(*,*) ''
         
      nstepInp=0
      val_ind = 0
      disp_ind = 0
      connExported = .false.
      
      do while (disp_ind < MIN(nstepTime,nstepTP))

         ! -------- read the solution at -------!
         ! --------- nodal/gauss points --------!

         ! ---------- IMPORT THE SOLUTION ----------!
         
         val_ind = val_ind + 1

         if(GaussOrNodal.EQ.'G') then
            error = 0
         
            if(commentsInData.EQ.'Y') then
               read(103,*,iostat=error) values_N, values_t
               if(error.NE.0) then
                  write(*,*) 'End of data file. Number of data blocks read:',val_ind
                  exit
               endif               
            else
               values_N = val_ind
               values_t = val_ind
            endif
            do iElem = 1, NELX      
               read(103,*,iostat=error) (gaussValues(iDOF,iElem),iDOF=1,nQuantDOF)
               if(error.NE.0) then
                  write(*,*) 'Solution file cannot be opened, or is corrupt @element',iElem
                  stop               
               endif  
            enddo
            
            if(error.NE.0) exit

         endif

         if(GaussOrNodal.EQ.'N') then
            
            error = 0

            if(commentsInData.EQ.'Y') then
               read(103,*,iostat=error) values_N, values_t
               if(error.NE.0) then
                  write(*,*) 'End of data file. Number of data blocks read:',val_ind
                  exit
               endif 
            else
               values_N = val_ind
               values_t = val_ind
            endif
            if(error.NE.0) then
               write(*,*) 'Solution file cannot be opened, or is empty'
               stop
            end if
            do iNode = 1, nNodes_file
               read(103,*,iostat=error) nodalValues(:,iNode)
               if(error.NE.0) then
                  write(*,*) 'Solution file cannot be opened, or is corrupt.'
                  stop               
               endif  
            enddo

            if(error.NE.0) exit
            
         endif
         
         ! -------- read node displacements --------!


         ! if plotting in deformed configuration, find the displacement field block corresponding to the solution
         if (plotConfig=='D') then
            foundDispBlock = .false.
            do while (.not.foundDispBlock)   ! search for the displacement data block that matches
                                             ! the time stamp of the solution
               disp_ind = disp_ind + 1

               if(plotConfig.EQ.'D' .and. commentsInDISP.EQ.'Y') then
                  read(101,*,iostat=error) nodalPositions_N, nodalPositions_t
                  if(error.EQ.-1) then
                     write(*,*) 'unexpected format/corrupt file: disp.out'
                     write(*,*) 'NSTEP time -- pair is expected before each data block'
                  endif
                  if (commentsInData=='Y') then    ! loop disp blocks until increment numbers match
                     if (nodalPositions_N == values_N) foundDispBlock = .true.
                  else                             ! assume data blocks in sync with disp blocks. accept disp block.
                     foundDispBlock = .true.
                  endif
               elseif(plotConfig.EQ.'D' .and. commentsInDISP.EQ.'N') then
                  ! no time stamp on the disp.out file. assume a one-to-one sequential correspondence between disp.out and time.out
                  if (commentsInData=='Y') then ! loop disp blocks until the corresp time in time.out matches time of the field block
                     if (times(disp_ind) == values_t) then
                        foundDispBlock = .true.
                     endif
                  else ! no time stamp in either disp or data block assume one-to-one sequential correspondence betw disp block and time 
                     foundDispBlock = .true.
                  endif
               endif
               
               do iNode=1,NX
                  read(101,*,iostat=error)(nodalPositions(iDOF,iNode), iDOF=1,3)           !only last step
                  if (error /= 0) exit
                  nodalPositions(1:3,iNode) = G0XYZ((iNode-1)*3+1:iNode*3) & !plot on deformed configuration
                                             +dispFactor*nodalPositions(1:3,iNode)
               enddo
               if (error /= 0) then
                  write(*,*) 'displacement file cannot be opened, or is empty'
                  stop
               endif

            enddo
         else ! plotConfig.EQ.'U' -- plot in undeformed configuration -- we're NOT reading disp.out file
            ! if not plotting in deformed configuration, plot using the initial positions
            disp_ind = disp_ind + 1
            do iNode=1,NX
               nodalPositions(1:3,iNode)=G0XYZ((iNode-1)*3+1:iNode*3)   !plot on reference configuration
            enddo    
         endif

         ! OUTPUT
         plotThisTime = .false.
         
         if (rangeOrSpecify=='R' .and. disp_ind >= beginStepTP) then
            plotThisTime = .true.
         elseif (rangeOrSpecify=='S') then
            found = .false.
            do i=1,nstepsSpecific
               if (plotTimeSteps(i)==disp_ind) found = .true.
            enddo
            if (found) plotThisTime = .true.
         endif
         
         if (rangeOrSpecify=='S') then
            if (disp_ind > MAXVAL(plotTimeSteps(:))) exit  ! no more time steps to plot
         endif

         ! now plot the field variables  !
         ! ------------------------------!

         if(plotThisTime.and.GaussOrNodal.EQ.'N') then
            ! EXPORT NODAL VALUES

            write(603,'(A)',advance='no') (' VARIABLES= "X","Y","Z"')
            do iDOF=1,nOutQuantDOF
               write(603,'(A,A,I0,A)',advance='no') ',"', solnFileName(1:3),outQuantDOF(iDOF),'"'
            enddo
            write(603,*)
             
            write(603,'(A,I0,A,I0,A)',advance='no') 'ZONE NODES=',nNodes_file,', ELEMENTS=',NELX, &
                                                    ', DATAPACKING=POINT ,ZONETYPE=FETETRAHEDRON, '

            write(603,'(A,F20.7)') 'STRANDID=1, SOLUTIONTIME=',values_t            
            
            do iNode=1,nNodes_file
               write(603,*) nodalPositions(:,grainNodes(iNode)), &
                           (nodalValues(outQuantDOF(j),iNode),j=1,nOutQuantDOF)
            enddo
            
            ! connectivity
            connZone = 1
            connExported = .true.
            if (valuesAtNodes==valuesUnique) then
               do iElem=1,NELX
                  write(603,*) IJK((iElem-1)*4+1:iElem*4)
               enddo
            elseif (valuesAtNodes==valuesPerGrain) then
               do iElem=1,NELX
                  write(603,*) IJK_grainNodes(:,iElem)
               enddo
            endif
            
         elseif(plotThisTime.and.GaussOrNodal.EQ.'G') then
         
            ! EXPORT GAUSS PT VALUES / CELL-CENTERED VALUES

            ! ---------------- HEADER ----------------- !
            ! ----------------------------------------- !
            write(603,'(A)',advance='no') (' VARIABLES= "X","Y","Z"')
            do iDOF=1,nOutQuantDOF
               write(603,'(A,A,I0,A)',advance='no') ',"', solnFileName(1:3),outQuantDOF(iDOF),'"'
            enddo
            write(603,*)
            
            write(603,'(A,I0,A,I0,A)',advance='no') &
            'ZONE NODES=',nTotGrainNodes,', ELEMENTS=',NELX,', ZONETYPE=FETETRAHEDRON, DATAPACKING=BLOCK, VARLOCATION=(['
            
            do i=4,nOutQuantDOF+2
               write(603,'(I0,A)',advance='no') i,','
            enddo
            write(603,'(I0,A)',advance='no') nOutQuantDOF+3,']=CELLCENTERED), '
            if(connExported) write(603,'(A,I0,A)',advance='no') 'CONNECTIVITYSHAREZONE=',connZone,', '
            write(603,'(A,F20.7)') 'STRANDID=1, SOLUTIONTIME=',values_t
                        
            ! ---------------- NODE POSITIONS ----------------- !
            ! ------------------------------------------------- !
            ! find the array index of the displacement data corresponding to this time step of the solution
            ! current time step of the solution: 
            do i=1,3
               do iNodeBlock=1,INT(nTotGrainNodes/ENTRIESPERLINE)
                  write(603,*) (nodalPositions(i,grainNodes(iNodeIndex)),   &
                                 iNodeIndex=(iNodeBlock-1)*ENTRIESPERLINE+1,iNodeBlock*ENTRIESPERLINE)
               enddo
               !if(iNodeBlock*ENTRIESPERLINE.NE.NX) then  ! print the remaining
                  write(603,*) (nodalPositions(i,grainNodes(iNodeIndex)),   &
                                 iNodeIndex=INT(nTotGrainNodes/ENTRIESPERLINE)*ENTRIESPERLINE+1,nTotGrainNodes)
               !end if
            enddo
            ! ---------------- FIELD VARIABLES ----------------- !
            ! -------------------------------------------------- !
            do j=1,nOutQuantDOF
               do iElemBlock=1,INT(NELX/ENTRIESPERLINE)
                  write(603,*) (gaussValues(outQuantDOF(j),iElem), &
                                iElem=(iElemBlock-1)*ENTRIESPERLINE+1,iElemBlock*ENTRIESPERLINE)
               enddo
               write(603,*) (gaussValues(outQuantDOF(j),iElem),iElem=INT(NELX/ENTRIESPERLINE)*ENTRIESPERLINE+1,NELX)
            enddo  
            ! ---------------- CONNECTIVITY ----------------- !
            ! ----------------------------------------------- !
            ! -------- connectivity information is ---------- !
            ! ------------imported from zone 1 -------------- !
            ! ----------- see CONNECTIVITYSHAREZONE --------- !
            ! ----------------------------------------------- !
            ! --------------only written in first zone------- !
            if(.not.connExported) then
               connZone = 1
               connExported = .true.
               do iElem=1,NELX
                  write(603,*) IJK_grainNodes(:,iElem)
               enddo
            endif            

         end if        
      
         !progress bar fancyness
         if((60*val_ind)/nstepTP.GT.lastProgress) then
            write(*,'(A)',advance='no') '*'
            call flush(6)
            lastProgress=(60*val_ind)/nstepTP
         endif
         
      enddo
      write(*,*)
      if (GaussOrNodal=='G') write(*,*) 'Exported tecPlot output file: ', trim(solnFileName)//'_TPG.dat'
      if (GaussOrNodal=='N') write(*,*) 'Exported tecPlot output file: ', trim(solnFileName)//'_TPN.dat'

      ! close files
      close(603)
700   if(plotConfig.EQ.'D') close(101)
      close(103)
      
      
      CALL grains_Destruct()
      ! element groupings
      deallocate(elemFaceNeighbors)
      
      deallocate(elemVol)
      
      deallocate(outQuantDOF,export_QuantDOF)
      deallocate(nodalPositions)
      if(GaussOrNodal.EQ.'N') deallocate(nodalValues)
      if(GaussOrNodal.EQ.'G') deallocate(gaussValues)
      
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