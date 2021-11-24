      program convertProcessorBasedOutputs
      implicit none
      integer:: nProcessors,I,J,JJ,I1,I2,IA,IPROC,NELEND,NELNOX,NELST,KOUNT,KINC_MAX,&
                writeOption,timeStep,istep,kinc,nElements,ierror,iElement,MaxElemInPatch,&
                NPATCH
      real(8) :: times(100000)
      integer, allocatable :: plotTimeSteps(:)
       integer :: nstepsSpecific
      logical :: plotThisTime, found
      real(8) :: dummy
      real(8),allocatable :: STRESS(:,:),FP(:,:), NYE(:,:)
      integer, allocatable :: NElementsInPatch(:),ElementsInPatch(:,:)
      logical :: putCommentsInData,processNyeToo
      character(len=10) :: strProcessor
      character(len=100) :: fileName
      logical :: oneComponentPerLine, fileExists

      write(*,*) 'detecting # of processors used for CPFE simulation...'
      nProcessors = 0
      fileExists=.true.
      do while (fileExists)
         CALL GET_SUBSCRIPT(nProcessors,strProcessor)
         fileName='stress'//strProcessor
         inquire(file=trim(fileName),exist=fileExists)
         if (fileExists) nProcessors = nProcessors + 1
      enddo
      
      if (nProcessors == 0) then
         write(*,*) 'output files are missing'
         stop
      endif
      
      write(*,*) 'No of processors used for CPFE simulation:', nProcessors
      
      write(*,*) 'Enter 0 ---> writes info for ALL time steps into a single file'
      write(*,*) 'Enter 1 ---> writes info for LAST time step into a single file'
      write(*,*) 'Enter 2 ---> writes info for PREFERED time step(s) into a single file'
      read(*,*) writeOption
      
      if (writeOption == 2) then
         write(*,*) 'Enter the total # of time steps you would like to export:'
         read(*,*) nstepsSpecific
         write(*,*)
         allocate(plotTimeSteps(nstepsSpecific))
         write(*,*) 'Enter the line numbers in time.out that you would like to export:'
         read(*,*) plotTimeSteps(1:nstepsSpecific)
         write(*,*)
      
!         write(*,*) 'Info for which time step do you prefer to be written?'
!         read(*,*) timeStep
      endif

      write(*,*) 'put time information before each data block in the data file? (T/F)'
      read(*,*) putCommentsInData

      write(*,*) 'process the Nye field too?'
      read(*,*) processNyeToo

      write(*,*) '*** reading time.out to determine # time steps ***'
      KINC_MAX=0
      OPEN(1,FILE='time.out')
   
      istep = 0
      do
         istep = istep + 1
         read(1,*,iostat=ierror) times(istep)
         if (ierror.ne.0) exit 
         KINC_MAX=KINC_MAX+1
      enddo
      close(1)
      write(*,*) '# time steps:',KINC_MAX
      write(*,*) '*** finished reading time.out  ***'      




      
      write(*,*) '*** reading PATCH.inp ***'
      OPEN(201,FILE='PATCH.inp')
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
      

      nElements=sum(NElementsInPatch)
      
      allocate(STRESS(nElements,6),FP(nElements,9))
      if (processNyeToo) allocate(NYE(nElements,9))
      STRESS = 0.d0
      FP = 0.d0
      if (processNyeToo) NYE = 0.d0
      

      DO IPROC=0,(nProcessors-1)
         CALL GET_SUBSCRIPT(IPROC,strProcessor)
         fileName='stress'//strProcessor
         OPEN(100+IPROC,FILE=fileName)
         fileName='fp'//strProcessor
         OPEN(7000+IPROC,FILE=fileName) 
         if (processNyeToo) then
            fileName='curlfp'//strProcessor
            OPEN(8000+IPROC,FILE=fileName) 
         endif
      ENDDO
      
      open(15,file='cauchy.out')
      open(16,file='fp.out')
      open(18,file='time_extracted.out')
      if (processNyeToo) open(17,file='nye.out')

      do kinc = 1 , KINC_MAX
      
         write(*,'(A31,I6,A5)') '*** reading info for time step:', kinc, '  ***'
     
         DO IPROC=0,(nProcessors-1)
                
            I1=NPATCH/nProcessors+1
            I2=NPATCH/nProcessors
            
            IA=NPATCH-I2*nProcessors
    
            
            IF((IPROC+1).LE.IA)THEN
               NELST=IPROC*I1+1
               NELNOX=I1
            ELSE
               NELST=IA*I1+(IPROC-IA)*I2+1
               NELNOX=I2
            ENDIF
            
            NELEND=NELST+NELNOX-1

               
            DO I=NELST,NELEND
            DO JJ=1,NElementsInPatch(I)
            
               DO J=1,6
                  READ(100+IPROC,*) STRESS(ElementsInPatch(I,JJ),J)
               ENDDO
               
            ENDDO
            ENDDO
         
         ENDDO
         
         
         DO IPROC=0,(nProcessors-1)
                
            I1=NPATCH/nProcessors+1
            I2=NPATCH/nProcessors
            
            IA=NPATCH-I2*nProcessors
    
            
            IF((IPROC+1).LE.IA)THEN
               NELST=IPROC*I1+1
               NELNOX=I1
            ELSE
               NELST=IA*I1+(IPROC-IA)*I2+1
               NELNOX=I2
            ENDIF
            
            NELEND=NELST+NELNOX-1
         
            
            DO I=NELST,NELEND
            DO JJ=1,NElementsInPatch(I)
               KOUNT=0
               DO I1=1,3
                  DO I2=1,3
                     KOUNT=KOUNT+1
                     READ(7000+IPROC,*) FP(ElementsInPatch(I,JJ),KOUNT)
                  ENDDO
               ENDDO
            ENDDO
            ENDDO
            
         ENDDO 
         

         
         DO IPROC=0,(nProcessors-1)
                
            I1=NPATCH/nProcessors+1
            I2=NPATCH/nProcessors
            
            IA=NPATCH-I2*nProcessors
    
            
            IF((IPROC+1).LE.IA)THEN
               NELST=IPROC*I1+1
               NELNOX=I1
            ELSE
               NELST=IA*I1+(IPROC-IA)*I2+1
               NELNOX=I2
            ENDIF
            
            NELEND=NELST+NELNOX-1
         
            if (processNyeToo) then
               DO I=NELST,NELEND
               DO JJ=1,NElementsInPatch(I)
                  KOUNT=0
                  DO I1=1,3
                     DO I2=1,3
                        KOUNT=KOUNT+1
                        READ(8000+IPROC,*) NYE(ElementsInPatch(I,JJ),KOUNT)
                     ENDDO
                  ENDDO
               ENDDO
               ENDDO
            endif
            
         ENDDO 

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
         
            write(18,*) times(kinc)
            if (putCommentsInData) then
               write(15,*) kinc, times(kinc)
               write(16,*) kinc, times(kinc)
               write(17,*) kinc, times(kinc)
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
         endif
         
         if (writeOption == 2 .and. plotTimeSteps(nstepsSpecific)==kinc) exit

      enddo
            
      close(15)
      close(16)
      close(17)
      close(18)

      DO IPROC=0,(nProcessors-1)
         close(100+IPROC)
         close(7000+IPROC)
         if (processNyeToo) close(8000+IPROC)
      ENDDO      
      
      deallocate(STRESS,FP)
      if (processNyeToo) deallocate(NYE)
            
            
      end
      
      
      
!*************************************************************
      SUBROUTINE GET_SUBSCRIPT(IPROC,strProcessor)
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*10 strProcessor
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
         strProcessor='_'//ACHAR(48+IN(1))//'.OUT'
      ELSEIF(I.EQ.2)THEN
         strProcessor='_'//ACHAR(48+IN(1))//ACHAR(48+IN(2))//'.OUT'
      ELSE
         strProcessor='_'//ACHAR(48+IN(1))//ACHAR(48+IN(2))//ACHAR(48+IN(3))//'.OUT'
      ENDIF
      
      
      RETURN
      END 
      