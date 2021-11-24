      program convertProcessorBasedOutputs
      implicit none
      integer:: nProcessors,nProcessors_OutputFiles, &
                I,J,JJ,I1,I2,IA,IPROC,KOUNT,KINC_MAX,&
                writeOption,timeStep,istep,kinc,ierror,iElement,iElementDummy,&
                MaxElemInPatch,NPATCH,nNodes,nElements,iNode,iNodeDummy
      integer, allocatable :: nElementsInProc(:)
      real(8) :: times(100000)
      integer, allocatable :: plotTimeSteps(:)
       integer :: nstepsSpecific
      logical :: plotThisTime, found, foundDispBlock, fileExists
      real(8) :: dummy
      real(8),allocatable :: STRESS(:,:),FP(:,:), NYE(:,:)
      real(8), allocatable :: disp(:,:)
      real(8) :: timeNodal
      integer, allocatable :: NElementsInPatch(:),ElementsInPatch(:,:)
      integer :: nComponents
      logical :: putCommentsInData,processNyeToo,writeTimeDispOut
      character(len=10) :: strProcessor
      character(len=100) :: fileName
      character(len=30) :: meshName
      
      write(*,*) 'detecting # of processors used for CPFE simulation...'
      nProcessors_OutputFiles = 0
      fileExists=.true.
      do while (fileExists)
         CALL GET_SUBSCRIPT(nProcessors_OutputFiles,strProcessor)
         fileName='stressCauchy'//strProcessor
         inquire(file=trim(fileName),exist=fileExists)
         if (fileExists) nProcessors_OutputFiles = nProcessors_OutputFiles + 1
      enddo
   
      if (nProcessors_OutputFiles == 0) then
         write(*,*) 'stressCauchy - output files are missing'
         stop
      endif
      
      write(*,*) 'No of processors used for CPFE simulation:', nProcessors_OutputFiles
      
      open(100,FILE='../postprocessing/elementDecompositionMetis.out')
      read(100,*) nProcessors
      if (nProcessors /= nProcessors_OutputFiles) then
         write(*,*) 'number of output files and the partitioning info in elementDecompositionMetis are inconsistent'
         stop         
      endif 
      allocate(nElementsInProc(0:nProcessors-1))
      do iProc=0,nProcessors-1
         read(100,*) nElementsInProc(iProc)
      enddo
      close(100)
      
      write(*,*) 'Enter 0 ---> writes info for ALL time steps into a single file'
      write(*,*) 'Enter 1 ---> writes info for LAST time step into a single file'
      write(*,*) 'Enter 2 ---> writes info for PREFERED time step(s) into a single file'
      read(*,*) writeOption
      
      if (writeOption == 2) then
         write(*,*) 'Enter the total # of time steps you would like to export:'
         read(*,*) nstepsSpecific
         write(*,*)
         allocate(plotTimeSteps(nstepsSpecific))
         write(*,*) 'Enter the line numbers in timeElemental.out that you would like to export:'
         read(*,*) plotTimeSteps(1:nstepsSpecific)
         write(*,*)
      
!         write(*,*) 'Info for which time step do you prefer to be written?'
!         read(*,*) timeStep
      endif

      write(*,*) 'put time information before each data block in the data file? (T/F)'
      read(*,*) putCommentsInData
      
      write(*,*) 'write time.out and disp.out too?'
      read(*,*) writeTimeDispOut
      
      !write(*,*) 'process the Nye field too?'
      !read(*,*) processNyeToo
      processNyeToo = .false.

      write(*,*) '*** reading timeElemental.out to determine # time steps ***'
      KINC_MAX=0
      OPEN(1,FILE='timeElemental.out')
   
      istep = 0
      do
         istep = istep + 1
         read(1,*,iostat=ierror) times(istep)
         if (ierror.ne.0) exit 
         KINC_MAX=KINC_MAX+1
      enddo
      close(1)
      write(*,*) '# time steps:',KINC_MAX
      write(*,*) '*** finished reading timeElemental.out  ***'      

      call readNumNodesElements(nNodes,nElements)  ! from the mesh file
      allocate(disp(3,nNodes))
      
      write(*,*) '*** reading PATCH.inp ***'
      OPEN(201,FILE='../input/PATCH.inp')
      READ(201,*) NPATCH,MaxElemInPatch
      
      ALLOCATE(NElementsInPatch(NPATCH))
      ALLOCATE(ElementsInPatch(NPATCH,MaxElemInPatch))
      
      NElementsInPatch=0
      ElementsInPatch=0
      
      DO I=1,NPATCH
         READ(201,*) NElementsInPatch(I),(ElementsInPatch(I,J),J=1,NElementsInPatch(I))      
      ENDDO
      
      CLOSE(201)          
      write(*,*) '*** finished reading PATCH.inp ***'      
      
      
      if (nElements /= sum(NElementsInPatch)) then
         write(*,*) '# of elements in the mesh is not equal to the total # elements in patches'
         stop
      endif
      
      allocate(STRESS(nElements,6),FP(nElements,9))
      if (processNyeToo) allocate(NYE(nElements,9))
      STRESS = 0.d0
      FP = 0.d0
      if (processNyeToo) NYE = 0.d0
      

      DO IPROC=0,(nProcessors-1)
         CALL GET_SUBSCRIPT(IPROC,strProcessor)
         fileName='stressCauchy'//strProcessor
         OPEN(100+IPROC,FILE=fileName)
         fileName='defGradPlastic'//strProcessor
         inquire(file=trim(fileName),exist=fileExists)
         if (.not.fileExists) fileName='deformationGradientPlastic'//strProcessor
         inquire(file=trim(fileName),exist=fileExists)
         if (.not.fileExists) then
            write(*,*) 'file not found: defGradPlastic_#.out or deformationGradientPlastic_#.out'
            stop
         endif
         OPEN(7000+IPROC,FILE=fileName) 
         if (processNyeToo) then
            fileName='curlfp'//strProcessor
            OPEN(8000+IPROC,FILE=fileName) 
         endif
      ENDDO
      
      open(15,file='cauchy.out')
      open(16,file='fp.out')
      if (processNyeToo) open(17,file='nye.out')
   
      open(9001,file='timeNodal.out')
      open(9002,file='displacements.out')
      if (writeTimeDispOut) then
         open(21,file='time.out')
         open(22,file='disp.out')
      endif

      do kinc = 1 , KINC_MAX
      
         write(*,'(A31,I6,A5)') '*** reading info for time step:', kinc, '  ***'
     
         DO IPROC=0,(nProcessors-1)
         
            DO iElementDummy=1,nElementsInProc(iProc)
               READ(100+IPROC,*) iElement, nComponents, STRESS(iElement,1:nComponents)
               if (nComponents /= 6) then
                  write(*,*) 'unexpected output file format. Cauchy Stress nComponents=',nComponents
                  stop
               endif
            ENDDO
         
         ENDDO
         
         
         DO IPROC=0,(nProcessors-1)
         
            DO iElementDummy=1,nElementsInProc(iProc)
               READ(7000+IPROC,*) iElement, nComponents, FP(iElement,1:nComponents)
               if (nComponents /= 9) then
                  write(*,*) 'unexpected output file format. FP nComponents=',nComponents
                  stop
               endif
            ENDDO
            
         ENDDO 
         

         
         if (processNyeToo) then
            DO IPROC=0,(nProcessors-1)

               DO iElementDummy=1,nElementsInProc(iProc)
                  READ(8000+IPROC,*) iElement, nComponents, NYE(iElement,1:nComponents)
                  if (nComponents /= 9) then
                     write(*,*) 'unexpected output file format. Cauchy Stress nComponents=',nComponents
                     stop
                  endif
               ENDDO
               
            ENDDO 
         endif

         plotThisTime = .false.
         if ( writeOption == 0 .or.                         &
             (writeOption == 1 .and. kinc == kinc_max)) then
         
            plotThisTime = .true.
         elseif (writeOption ==2 ) then
            found = .false.
            do i=1,nstepsSpecific
               if (plotTimeSteps(i)==kinc) found = .true.
            enddo
            if (found) plotThisTime = .true.
         endif
         
         if (plotThisTime) then
         
            if (putCommentsInData) then
               write(15,*) kinc, times(kinc)
               write(16,*) kinc, times(kinc)
               if (processNyeToo) write(17,*) kinc, times(kinc)
            endif
            do iElement=1,nElements
               write(15,*) STRESS(iElement,1:6)
            enddo
            
            do iElement=1,nElements
               write(16,*) FP(iElement,1:9)
            enddo

            if (processNyeToo) then
               do iElement=1,nElements
                  write(17,*) NYE(iElement,1:9)
               enddo
            endif
            
            ! find the displacements that corresp. to this solution time
            foundDispBlock = .false.
            do while (.not.foundDispBlock)   ! search for the displacement data block that matches
                                             ! the time stamp of the solution
               read(9001,*) timeNodal
               do iNodeDummy=1,nNodes
                  read(9002,*) iNode, disp(:,iNode)
               enddo
               if (timeNodal == times(kinc)) then
                  foundDispBlock = .true.
                  exit
               endif
            enddo
            if (.not. foundDispBlock) then
               write(*,*) 'cannot find displacement block corresp. to the solution time:',times(kinc)
               stop
            endif

            if (writeTimeDispOut) then
               write(21,*) times(kinc)
               write(22,*) kinc, times(kinc)
               do iNode=1,nNodes
                  write(22,*) disp(:,iNode)
               enddo
            endif
            
         endif
         
         if (writeOption == 2) then
            if (plotTimeSteps(nstepsSpecific)==kinc) exit
         endif

      enddo
            
      close(15)
      close(16)
      if (processNyeToo) close(17)
   
      close(9001)
      close(9002)
      if (writeTimeDispOut) then
         close(21)
         close(22)
      endif


      DO IPROC=0,(nProcessors-1)
         close(100+IPROC)
         close(7000+IPROC)
         if (processNyeToo) close(8000+IPROC)
      ENDDO      
      
      deallocate(STRESS,FP)
      if (processNyeToo) deallocate(NYE)

	write(*,*) 'exported cauchy.out, fp.out, nye.out, time.out'
            
            
      end
      
      
      
!*************************************************************
      SUBROUTINE GET_SUBSCRIPT(IPROC,FNAME)
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*10 FNAME
      DIMENSION IN(10),INTMP(10)           
      
      IPROC1=IPROC
      I=0

      DO WHILE(IPROC1.GT.9)
         I1=IPROC1/10
         I2=IPROC1-I1*10
         I=I+1
         IPROC1=I1
         IN(I)=I2
      ENDDO

      I=I+1
      IN(I)=IPROC1

      DO J=1,I
         INTMP(I-J+1)=IN(J)
      ENDDO

      DO J=1,I
         IN(J)=INTMP(J)
      ENDDO

      IF(I.EQ.1)THEN
         FNAME='_'//ACHAR(48+IN(1))//'.out'
      ELSEIF(I.EQ.2)THEN
         FNAME='_'//ACHAR(48+IN(1))//ACHAR(48+IN(2))//'.out'
      ELSE
         FNAME='_'//ACHAR(48+IN(1))//ACHAR(48+IN(2))//ACHAR(48+IN(3))//'.out'
      ENDIF
      
      
      RETURN
      END 
      
      SUBROUTINE readNumNodesElements(nNodes,nElements)
      implicit none
      
      integer, intent(out) :: nNodes, nElements
      
      character(len=150) :: lineStr
      integer :: error, icyc, MDIM,NDFADD, iNode, lenTrimmed
      real(8) :: xNode(3)
      character(len=30) :: meshName
      character(len=50) :: fileName
      
      open(201,file='loadtime.inp')
      call READSTR(201,lineStr,error)
      read(201,*) icyc
      call READSTR(201,lineStr,error)
      read(201,*) meshName
      close(201)
      
      lenTrimmed=len_trim(meshName)
      fileName=meshName(1:lenTrimmed)//'.inp'
      open(unit=202,file=fileName)

      call READSTR(202,lineStr,error)
      read(202,*) MDIM,NDFADD
      call READSTR(202,lineStr,error)
      read(202,*) nNodes
      do iNode=1,nNodes
         read(202,*) xNode(1:3)
      enddo
      call READSTR(202,lineStr,error)
      read(202,*) nElements
      close(202)
      
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
            
      subroutine doPartitioningLinear        &
      (                                      &
         nItems,                             &
         idProcessor,                        &
         nProcessors,                        &
         idStart,                            &
         idEnd                               &                     
      )

         implicit none

      !  Inputs ------------------------------------------------------------------------------------------

         integer, intent(in) ::     &
            nItems,                             &
            idProcessor,                        &
            nProcessors

      !  Outputs -----------------------------------------------------------------------------------------

         integer, intent(out) ::    &
            idStart,                            &
            idEnd                            

      ! --------------------------------------------------------------------------------------------------
      !  Local Variables ---------------------------------------------------------------------------------
      ! --------------------------------------------------------------------------------------------------

         integer ::     &
            nItemsInProcessor,      &
            nProcessorsCritical

      ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         nProcessorsCritical = nItems - (nItems / nProcessors) * nProcessors

         if ( (idProcessor + 1) <= nProcessorsCritical) then
            nItemsInProcessor = nItems / nProcessors + 1
            idStart = idProcessor * nItemsInProcessor + 1

         else 
            nItemsInProcessor = nItems / nProcessors 
            idStart = nProcessorsCritical * (nItemsInProcessor + 1) +                                    &
                                             (idProcessor - nProcessorsCritical) * nItemsInProcessor + 1
            
         end if

         idEnd = idStart + nItemsInProcessor - 1

      end subroutine doPartitioningLinear
