!  This program is a PARALLEL code which reads field data from the user-defined files and then calculates
!  the average value of the field over grains.
!
!  Author: Ahmad Shahba (2015)
!
! The following files are required for this program
!
!  1) loadtime.inp
!  2) Mesh file whose name is given in "loadtime.inp"
!  3) PATCH.inp (if multiProcessorInput)
!  4) time.out
!  5) optionsCI.inp (optional)
!  6) grainTexture.inp
!  7) elementGrains.inp or grains.inp (optional)

      program calculateAverageFieldOverGrains
      use optionsCI
      use meshtools
      use grains
      implicit none

      interface
         subroutine READGM(G0XYZ,IJK,NX,NELX,iproc_mpi)
            implicit none
            real(8), allocatable :: G0XYZ(:)
            integer, allocatable :: IJK(:)
            integer, intent(out) :: NX, NELX
            !locals
            integer :: IJKNEL(4)
            real(8) :: xc(3)
            integer :: MDIM,NDFADD,NDF,NODE,NGAUSS
            character(len=50) :: fileName
            !dummies
            integer :: icyc,len,i,ii,j,jj,nnum,iElem,iNode,NEL,NEQ,iproc_mpi
         end subroutine READGM
      
      end interface





      include 'mpif.h'       
      integer:: ierror, iproc_mpi, nprocs_mpi, NPROCS_SIMULATION, nstepMax, nstepTimeOUT, istep , &
                error, dummyInteger, iproc, proc_FIRST, proc_LAST, proc_IN_PROCESSOR, NELST, &
                NELEND, NELNOX, NPATCH, i, ii, j, jj, MaxElemInPatch, iElem, NX, NELX, MDIM, NDF, &
                NGAUSS, NODE, NFACES, NEDGES, NEQ, nQuantDOF, grainID, staIdx, endIdx, &
                staGrainID, endGrainID, nGrainsInProcessor, nGrainsToExport
      real(8):: dummyReal, timeStart, elemVol, grainVolume
      integer :: grainIDsToExport(3000)
      real(8):: times(100000), FP(3,3), eye(3,3), dummy(3,3), cauchy(3,3), cAxis(3), traction(3)
      integer, allocatable:: NElementsInPatch(:), ElementsInPatch(:,:), IJK(:)
      real(8), allocatable:: FieldValuePerElement(:,:) , G0XYZ(:), volumeElements(:), &
                             fieldAverageOverGrain(:,:), plasticStrainElement(:), &
                             plasticStrainOverGrain(:), stressBasalElement(:), &
                             stressBasalOverGrain(:), grain_rot(:,:,:)
      character(len=30) :: meshName, fileName, solnFileName
      character(len=10) :: strProcNumber
      logical :: componentsAppearOn1Line, computePlasticStrain, computeStressBasal, &
                 exportForCertainGrains 

      
 


 
      call mpi_init(ierror)
      call mpi_comm_rank(mpi_comm_world,iproc_mpi,ierror)
      call mpi_comm_size(mpi_comm_world,nprocs_mpi,ierror)   
      
      if (iproc_mpi ==0) write(*,*) '# processors used to run this program: ',nprocs_mpi
   ! import CI program options
      call readOptionsCI('optionsCI.inp')
      if (IPROC_MPI==0) CALL summarizeOptionsCI()   

   ! checking the number of processors vs the number of input files
      if (.not.multiProcessorInput .and. nprocs_mpi /= 1 ) then
         if (iproc_mpi ==0) &
            write(*,*) 'multiProcessorInput =.false. ; therefore, run this code with only one processor'
         call mpi_finalize(ierror)
         stop
      endif
      

   ! receiving input from the terminal   
      if (iproc_mpi == 0) then
         write(*,*) 'Specify the file that contains the field data over elements?'
         
         if (multiProcessorInput) then
            write(*,*) 'NOTE: just provide the main part of file name '
            write(*,*) '      For example, if your files are stress_0.out, stress_1.out and etc, here insert only stress'
         else
            write(*,*) 'NOTE: Insert the file name with extension, e.g. stress.out'
         endif
         
         read(*,*) solnFileName
         write(*,*)
         
         
         write(*,*) 'Enter # of components of the field variable.'
         write(*,*) '(Temperature: 1, Cauchy stress: 6, Fp: 9, etc.)'
         read(*,*) nQuantDOF
         write(*,*)
         
         write(*,*) 'All components of the field per element appear in 1 line line? (T/F)'
         read(*,*) componentsAppearOn1Line
         write(*,*)   
         
         write(*,*) 'If input files are Fp, would you like the code to compute average plastic' 
         write(*,*) '    strain as well? (T/F)'
         read(*,*) computePlasticStrain
         write(*,*)
         
         
         write(*,*) 'If input files are stress, would you like the code to compute average stress'
         write(*,*) '   normal to basal plane, as well? (T/F)'
         read(*,*) computeStressBasal
         write(*,*)
         

         write(*,*) 'Do you want to have average values for certain grains in separate files? (T/F)'
         read(*,*) exportForCertainGrains
         if (exportForCertainGrains) then
            write(*,*) 'How many grains do you want to track?'
            read(*,*) nGrainsToExport
            write(*,*) 'Which grains do you like to track?'
            grainIDsToExport = 0
            do i = 1, nGrainsToExport
               read(*,*) grainIDsToExport(i)
            enddo
         endif
         write(*,*)
         
                  
      

         if(multiProcessorInput) then   
            
            ! find the number of CPUs used in simulation/data export
            CALL detectNumProcessors(NPROCS_SIMULATION,solnFileName)

            if (NPROCS_SIMULATION ==0) then
               if (IPROC_MPI==0) write(*,*) 'cannot locate multiple-CPU solution data:' &
                                             //trim(adjustl(solnFileName))//'_#.OUT'
               CALL MPI_FINALIZE(IERROR)
               stop            
            else
               if (IPROC_MPI==0) write(*,*) '# of CPUs/export files detected:',NPROCS_SIMULATION
            endif     
            
         endif
         
      endif
      
      call mpi_barrier(mpi_comm_world,ierror)
      call mpi_bcast(nQuantDOF,1,mpi_integer,0,mpi_comm_world,ierror)
      call mpi_bcast(NPROCS_SIMULATION,1,mpi_integer,0,mpi_comm_world,ierror)
      call mpi_bcast(componentsAppearOn1Line,1,mpi_logical,0,mpi_comm_world,ierror)
      call mpi_bcast(computePlasticStrain,1,mpi_logical,0,mpi_comm_world,ierror)
      call mpi_bcast(computeStressBasal,1,mpi_logical,0,mpi_comm_world,ierror)
      call mpi_bcast(exportForCertainGrains,1,mpi_logical,0,mpi_comm_world,ierror)
      if (exportForCertainGrains) then
         call mpi_bcast(nGrainsToExport,1,mpi_integer,0,mpi_comm_world,ierror)
         call mpi_bcast(grainIDsToExport,3000,mpi_integer,0,mpi_comm_world,ierror)
      endif
      call mpi_bcast(solnFileName,30,mpi_character,0,mpi_comm_world,ierror)
      
   ! checking the number of processors vs the number of input files
      if (multiProcessorInput .and. nprocs_mpi > NPROCS_SIMULATION ) then
         if (iproc_mpi ==0) &
            write(*,*) '# processors used for this code shall be less than the number of input files'
         call mpi_finalize(ierror)
         stop
      endif

      
      ! opening input files
      if (multiProcessorInput) then
         do iproc= 0, nprocs_simulation-1
            call get_subscript(iproc,strprocnumber)
            fileName = trim(adjustl(solnFileName))//strprocnumber
            call openFile(500+iproc,fileName,iproc_mpi)
         enddo      
      else
         call openFile(500,solnFileName,iproc_mpi)
      endif
      
      

      ! load F-Bar patch information
      if (multiProcessorInput) then
         fileName = 'PATCH.inp'
         call openFile(201,fileName,iproc_mpi)
         READ(201,*) NPATCH,MaxElemInPatch
         
         ALLOCATE(NElementsInPatch(NPATCH))
         ALLOCATE(ElementsInPatch(NPATCH,MaxElemInPatch))
         
         NElementsInPatch=0
         ElementsInPatch=0
         
         DO I=1,NPATCH
            READ(201,*) NElementsInPatch(I),(ElementsInPatch(I,J),J=1,NElementsInPatch(I))      
         ENDDO
         
         CLOSE(201)        
      endif

      ! READ IN THE MESH, CONNECTIVITY
      CALL READGM(G0XYZ,IJK,NX,NELX,iproc_mpi)      
      ! set element defaults
      MDIM=3
      NDF=3
      NGAUSS=1
      NODE=4
      NFACES=4
      NEDGES=6
      NEQ=NX*NDF

      timeStart = MPI_WTIME()
      
      if(IPROC_MPI==0) then
         write(*,*) 'Mesh imported:'
         write(*,*) '# of elements: ', NELX
         write(*,*) '# of nodes: ', NX
      endif

      
      ! calculating volumes of elements
      allocate (volumeElements(NELX))
      volumeElements = 0.d0
      do iElem = 1 , NELX
         call calcElemVol(iElem, elemVol, IJK, G0XYZ, NODE)
         volumeElements(iElem) = dabs(elemVol)
      enddo
      
      ! ---------------------- IMPORT GRAINS ------------------------!
      ! first try grains.inp, if fails, elementGrains.inp, if fails create a trivial grain.
      ! import grains
      ! first try grains.inp, if fails, elementGrains.inp, if fails create a trivial grain.
      CALL readGrains('grains.inp',NX,NELX)
      if(grainsImported .and. iproc_mpi==0) &
         write(*,'(A,I0,A)') ' ', nGrains,' grains imported from grains.inp'   

      if(.NOT.grainsImported) then
         CALL readGrainsFromElementList('elementGrains.inp',NX,NELX)
         if(grainsImported .and. iproc_mpi==0) &
            write(*,'(A,I0,A)') ' ', nGrains,' grains imported from elementGrains.inp'
      endif
      
      if(.NOT.grainsImported) then
         CALL createSingleGrain(NX,NELX)
         if (iproc_mpi==0) write(*,'(A,I0,A)') 'single grain assumed'
      endif

      !read texture: grainTexture(:)
      CALL readTexture('grainTexture.inp')
      if(textureImported .and. iproc_mpi==0 ) write(*,*) 'texture imported'

      ! calculate grain rotation matrices -  Deniz - Sept 22
      allocate(grain_rot(3,3,nGrains))
      do grainID=1,nGrains
         call eulerRot(grainTexture(:,grainID),grain_rot(:,:,grainID))
         if (EulerDekaConvention) then
            grain_rot(:,:,grainID) = TRANSPOSE(grain_rot(:,:,grainID))
         endif
      ENDDO      
      
      ! READ TIME STEPPING INFORMATION
      fileName = 'time.out'
      call openFile(102,fileName,iproc_mpi)
      nstepTimeOUT=0
      times(:) = 0.d0
      do
         read(102,*,iostat=error) times(nstepTimeOUT+1)
         if(error.EQ.-1) EXIT
         nstepTimeOUT=nstepTimeOUT+1
      enddo
      close(102)
      if (nstepTimeOUT.EQ.0) then
         if(IPROC_MPI==0) write(*,*) 'Cannot find file: time.out'
         CALL MPI_FINALIZE(IERROR)
         STOP
      endif
      
      if(IPROC_MPI==0) write(*,*) '# of time steps in time.out:', nstepTimeOUT
      
      ! determine the # of time steps to analyse
      nstepMax = nStepsToAnalyse
      if (nStepsToAnalyse==0) nstepMax = nstepTimeOUT
      if (nStepsToAnalyse > nstepTimeOUT) nstepMax = nstepTimeOUT

      ! checking the number of processors
      if (nprocs_mpi > NELX .or. nprocs_mpi > nGrains ) then
         if (iproc_mpi ==0) &
            write(*,*) '# processors used for this code shall be less than # Elements and # Grains'
         call mpi_finalize(ierror)
         stop
      endif

       ! identity matrix
       EYE = 0.D0
       EYE(1,1)= 1.D0
       EYE(2,2)= 1.D0
       EYE(3,3)= 1.D0      
      
  
      ! opening output files
      allocate (FieldValuePerElement(nQuantDOF,NELX))
      allocate (fieldAverageOverGrain(nQuantDOF,nGrains))
      fileName='AverageGrainField.out'
      if (iproc_mpi == 0) call openFile(100,fileName,iproc_mpi)
      
      if (computePlasticStrain) then
         fileName='grainAvgPlasticStrain.out'
         call openFile(101,fileName,iproc_mpi)
         allocate(plasticStrainElement(NELX),plasticStrainOverGrain(nGrains))
      endif


      if (computeStressBasal) then
         fileName='AverageGrainBasalStress.out'
         call openFile(101,fileName,iproc_mpi)
         allocate(stressBasalElement(NELX),stressBasalOverGrain(nGrains))
         if (iproc_mpi == 0) open(105,file='stressBasalElement.out')
         
      endif
      
     
      if (iproc_mpi == 0 .and. exportForCertainGrains) then
         do i = 1, nGrainsToExport
            write(fileName,'(I4)') grainIDsToExport(i)
            fileName='AverageGrainField_'//trim(adjustl(fileName))//'.out'
            call openFile(1000+i,fileName,iproc_mpi)
         enddo
         
         if (computePlasticStrain) then
            do i = 1, nGrainsToExport
               write(fileName,'(I4)') grainIDsToExport(i)
               fileName='grainAvgPlasticStrain_'//trim(adjustl(fileName))//'.out'
               call openFile(4000+i,fileName,iproc_mpi)
            enddo
         endif

         if (computeStressBasal) then
            do i = 1, nGrainsToExport
               write(fileName,'(I4)') grainIDsToExport(i)
               fileName='AverageGrainBasalStress_'//trim(adjustl(fileName))//'.out'
               call openFile(4000+i,fileName,iproc_mpi)
            enddo
         endif
         
      endif
      
   ! reading files 
      do istep = 1, nstepMax
         
         if (iproc_mpi == 0) write(*,*) 'time:', istep, times(istep)
         
         if(multiProcessorInput) then
         
            call PARTITION_CALCULATOR(proc_FIRST,proc_LAST,proc_IN_PROCESSOR,&
                                      NPROCS_SIMULATION,NPROCS_MPI,IPROC_MPI) 
            proc_FIRST=proc_FIRST-1
            proc_LAST=proc_LAST-1

            FieldValuePerElement=0.d0

            
            do iproc=proc_first,proc_last

               error = 0
               if (commentsInDataFile) read(500+iproc,*,iostat=error) dummyInteger, dummyReal
               if(error.NE.0) then   !end of input file      
                  write(*,*) 'End of input file'
                  CALL MPI_FINALIZE(IERROR)
                  stop
               endif

               CALL PARTITION_CALCULATOR(NELST,NELEND,NELNOX,NPATCH,NPROCS_SIMULATION,IPROC)        
            
               DO I=NELST,NELEND
               DO JJ=1,NElementsInPatch(I)
               
                  if (componentsAppearOn1Line) then
                     read(500+iproc,*,iostat=error) &
                        FieldValuePerElement(1:nQuantDOF,ElementsInPatch(I,JJ))
                     if(error.NE.0) then   !end of input file      
                        write(*,*) 'End of input file'
                        CALL MPI_FINALIZE(IERROR)
                        stop
                     endif  
 
                  else
                     DO J=1,nQuantDOF
                        READ(500+IPROC,*,iostat=error) &
                           FieldValuePerElement(J,ElementsInPatch(I,JJ))
                        if(error.NE.0) then   !end of input file      
                           write(*,*) 'End of input file'
                           CALL MPI_FINALIZE(IERROR)
                           stop
                        endif                

                     ENDDO                   
                  endif
                  
               ENDDO
               ENDDO
               
            enddo

            
            CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
            call MPI_ALLREDUCE(MPI_IN_PLACE,FieldValuePerElement,NELX*nQuantDOF,&
                               MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR)  
        
         else 
         
            FieldValuePerElement = 0.d0
   
            error = 0
            if (commentsInDataFile) read(500,*,iostat=error) dummyInteger, dummyReal
            if(error.NE.0) then   !end of input file      
               write(*,*) 'End of input file'
               CALL MPI_FINALIZE(IERROR)
               stop
            endif                

            
            do iElem = 1, NELX      
               if (componentsAppearOn1Line) then
                  read(500,*,iostat=error) FieldValuePerElement(1:nQuantDOF,iElem)
                  if(error.NE.0) then   !end of input file      
                     write(*,*) 'End of input file'
                     CALL MPI_FINALIZE(IERROR)
                     stop
                  endif                

               else
                  DO J=1,nQuantDOF
                     READ(500,*,iostat=error) FieldValuePerElement(J,iElem)
                     if(error.NE.0) then   !end of input file      
                        write(*,*) 'End of input file'
                        CALL MPI_FINALIZE(IERROR)
                        stop
                     endif                

                  ENDDO               
               endif
            enddo

            
         endif
 
         ! compute effective plastic strain
         if (computePlasticStrain) then
            plasticStrainElement = 0.d0
            CALL PARTITION_CALCULATOR(NELST,NELEND,NELNOX,NELX,NPROCS_MPI,IPROC_MPI) 
            
            do iElem = NELST, NELEND
               FP(1,1)=FieldValuePerElement(1,iELem)
               FP(1,2)=FieldValuePerElement(2,iELem)
               FP(1,3)=FieldValuePerElement(3,iELem)
               FP(2,1)=FieldValuePerElement(4,iELem)
               FP(2,2)=FieldValuePerElement(5,iELem)
               FP(2,3)=FieldValuePerElement(6,iELem)
               FP(3,1)=FieldValuePerElement(7,iELem)
               FP(3,2)=FieldValuePerElement(8,iELem)
               FP(3,3)=FieldValuePerElement(9,iELem)               
               DUMMY=1./2.*(MATMUL(TRANSPOSE(FP),FP)-EYE)
               DUMMY=MATMUL(TRANSPOSE(DUMMY),DUMMY)
               plasticStrainElement(iElem)= DSQRT((DUMMY(1,1)+DUMMY(2,2)+DUMMY(3,3))*2./3.)               
            enddo
            
            
            CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)      
            call MPI_ALLREDUCE(MPI_IN_PLACE,plasticStrainElement,NELX,MPI_DOUBLE_PRECISION, &
                               MPI_SUM,MPI_COMM_WORLD,IERROR)  
            
         endif
         
         
         ! compute Stress normal to Basal plane
         if (computeStressBasal) then
            stressBasalElement = 0.d0
            CALL PARTITION_CALCULATOR(NELST,NELEND,NELNOX,NELX,NPROCS_MPI,IPROC_MPI) 
            
            do iElem = NELST, NELEND
               cauchy(1,1)=FieldValuePerElement(1,iELem)
               cauchy(2,2)=FieldValuePerElement(2,iELem)
               cauchy(3,3)=FieldValuePerElement(3,iELem)
               cauchy(1,2)=FieldValuePerElement(4,iELem)
               cauchy(1,3)=FieldValuePerElement(5,iELem)
               cauchy(2,3)=FieldValuePerElement(6,iELem)
               cauchy(2,1)=cauchy(1,2)
               cauchy(3,1)=cauchy(1,3)
               cauchy(3,2)=cauchy(2,3)
               grainID = grainIDelem(iElem)
               cAxis= grain_rot(1:3,3,grainID)
               traction= matmul(cauchy,cAxis)
               stressBasalElement(iElem)= dot_product(traction,cAxis)             
            enddo
            
            
            CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)      
            call MPI_ALLREDUCE(MPI_IN_PLACE,stressBasalElement,NELX,MPI_DOUBLE_PRECISION, &
                               MPI_SUM,MPI_COMM_WORLD,IERROR)  
                               
            if (iproc_mpi == 0) then
            
               do iElem=1,NELX
                  write(105,*) stressBasalElement(iElem)
               enddo
            
            endif
            
         endif
         
         ! calculating average over grains
         fieldAverageOverGrain = 0.d0
         
         if (computePlasticStrain) plasticStrainOverGrain = 0.d0
         if (computeStressBasal)   stressBasalOverGrain = 0.d0

         CALL PARTITION_CALCULATOR(staGrainID,endGrainID,nGrainsInProcessor,nGrains,&
                                   NPROCS_MPI,IPROC_MPI) 
                                   
         do grainID = staGrainID , endGrainID
     
            grainVolume = 0.d0
            staIdx = 1
            if ( grainID /= 1) staIdx =  grainElementIndx( grainID - 1 )
            endIdx = grainElementIndx( grainID ) - 1

            do ii = staIdx , endIdx
               fieldAverageOverGrain(1:nQuantDOF, grainID) = &
                  fieldAverageOverGrain(1:nQuantDOF, grainID) + &
                  FieldValuePerElement(1:nQuantDOF,grainElements(ii)) * &
                  volumeElements(grainElements(ii))
                  
                  
                  
               if (computePlasticStrain)  plasticStrainOverGrain(grainID) = &
                                          plasticStrainOverGrain(grainID) + &
                                          plasticStrainElement(grainElements(ii)) * &
                                          volumeElements(grainElements(ii))
                                          
               if (computeStressBasal)  stressBasalOverGrain(grainID) = &
                                          stressBasalOverGrain(grainID) + &
                                          stressBasalElement(grainElements(ii)) * &
                                          volumeElements(grainElements(ii))                                          
               
               grainVolume = grainVolume + volumeElements(grainElements(ii))
               
            enddo
              
            fieldAverageOverGrain(1:nQuantDOF, grainID) = &
               fieldAverageOverGrain(1:nQuantDOF, grainID) / grainVolume
              

            if (computePlasticStrain)  plasticStrainOverGrain(grainID) = &
                                       plasticStrainOverGrain(grainID) / grainVolume
                                       
            if (computeStressBasal)  stressBasalOverGrain(grainID) = &
                                       stressBasalOverGrain(grainID) / grainVolume                                       
         enddo
         
         
         CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)      
         call MPI_ALLREDUCE(MPI_IN_PLACE,fieldAverageOverGrain,nGrains*nQuantDOF, &
                            MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR) 
                            
         if (computePlasticStrain) then
            CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR) 
            call MPI_ALLREDUCE(MPI_IN_PLACE,plasticStrainOverGrain,nGrains,MPI_DOUBLE_PRECISION, &
                               MPI_SUM,MPI_COMM_WORLD,IERROR)                                
         endif

         if (computeStressBasal) then
            CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR) 
            call MPI_ALLREDUCE(MPI_IN_PLACE,stressBasalOverGrain,nGrains,MPI_DOUBLE_PRECISION, &
                               MPI_SUM,MPI_COMM_WORLD,IERROR)                                
         endif         
 

         ! output results for the time step
         if (iproc_mpi == 0) then
            write(100,'(I6,F18.8)') istep, times(istep)
            if (computePlasticStrain) write(101,'(I6,F18.8)') istep, times(istep)
            if (computeStressBasal)   write(101,'(I6,F18.8)') istep, times(istep)
            
            do grainID = 1, nGrains
               write(100,*) grainID, fieldAverageOverGrain(1:nQuantDOF,grainID) 
               if (computePlasticStrain) write(101,*) grainID, plasticStrainOverGrain(grainID) 
               if (computeStressBasal)   write(101,*) grainID, stressBasalOverGrain(grainID) 
            enddo
            
            
            if (exportForCertainGrains) then
               do i= 1, nGrainsToExport
                  grainID= grainIDsToExport(i)
                  write(1000+i,*) times(istep),fieldAverageOverGrain(1:nQuantDOF,grainID) 
                  if (computePlasticStrain) write(4000+i,*) times(istep), plasticStrainOverGrain(grainID) 
                  if (computeStressBasal)   write(4000+i,*) times(istep), stressBasalOverGrain(grainID) 
               enddo
            endif
            
         endif
                  
         
      enddo
      
      if (iproc_mpi ==0) then
         write(*,*)
         write(*,*) 'The following output files were generated:'
         write(*,*) '1. AverageGrainField.out'
         if (computePlasticStrain) write(*,*) '2. grainAvgPlasticStrain.out'
         if (computeStressBasal)   write(*,*) '2. AverageGrainBasalStress.out'
         if (exportForCertainGrains) then
            write(*,*)
            write(*,*) 'and the following grain-specific file(s)'
            do i = 1, nGrainsToExport
               write(*,*)
               write(fileName,'(I4)') grainIDsToExport(i) 
               fileName='AverageGrainField_'//trim(adjustl(fileName))//'.out' 
               write(*,*) fileName
               if (computePlasticStrain) then
                  write(fileName,'(I4)') grainIDsToExport(i); 
                  fileName='grainAvgPlasticStrain_'//trim(adjustl(fileName))//'.out'
                  write(*,*) fileName
               endif
                  
               if (computeStressBasal) then
                  write(fileName,'(I4)') grainIDsToExport(i)
                  fileName='AverageGrainBasalStress_'//trim(adjustl(fileName))//'.out'
                  write(*,*) fileName               
               endif
               
            enddo
            
         endif
      endif
      
      if (iproc_mpi == 0) close(100)
      
      if (computePlasticStrain .or. computeStressBasal) close(101)
      if (computeStressBasal .and. iproc_mpi==0) close(105)
      
      if (exportForCertainGrains) then
         do i= 1, nGrainsToExport
            close(1000+i)
            close(4000+i)
         enddo
      endif
      
      if (multiProcessorInput) then
         do iproc= 0, nprocs_simulation-1
            close(500+iproc)
         enddo           
      else
         close(500)
      endif
      
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR) 
      call mpi_finalize(ierror)
      stop
      end
      
      
      
!###################################################################################################      
      subroutine openFile(unitNumber,fileName,iproc_mpi)
      implicit none
      include 'mpif.h'
      integer, intent (in):: unitNumber, iproc_mpi
      character(len=30), intent(in) :: fileName
      integer :: ierror
      character (len=30) :: AlreadyExistingFileName
      logical:: fileAlreadyOpened
      
      fileAlreadyOpened = .false.
      inquire(unit=unitNumber,NAME=AlreadyExistingFileName , OPENED= fileAlreadyOpened)
      if (fileAlreadyOpened .and. iproc_mpi ==0) then
         write(*,*) 'File Number', unitNumber, 'already exists! Try another file unit number'
         write(*,*) 'File name is ',AlreadyExistingFileName
         call mpi_finalize(ierror)
         stop
      endif

      open (unitNumber,file=fileName)      
      
      end subroutine
      
!###################################################################################################
      SUBROUTINE detectNumProcessors(NPROCS_SIMULATION,solnFileName)
      implicit none
      ! arguments
      integer, intent(out) :: NPROCS_SIMULATION
      ! locals
      logical :: fileExists
      character(len=10) :: strProcNumber
      character(len=30) :: fileName, solnFileName
      
      NPROCS_SIMULATION = 0
      fileExists =.true.
      do while (fileExists)
         CALL GET_SUBSCRIPT(NPROCS_SIMULATION,strProcNumber)
         fileName=trim(adjustl(solnFileName))//strProcNumber
         inquire(FILE=fileName,exist=fileExists)
         if (fileExists) NPROCS_SIMULATION = NPROCS_SIMULATION + 1
      enddo
      
      END SUBROUTINE
      
   !###################################################################################################   
     SUBROUTINE GET_SUBSCRIPT(IPROC,strProcNumber)
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*10 strProcNumber
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
         strProcNumber='_'//ACHAR(48+IN(1))//'.OUT'
      ELSEIF(I.EQ.2)THEN
         strProcNumber='_'//ACHAR(48+IN(1))//ACHAR(48+IN(2))//'.OUT'
      ELSE
         strProcNumber='_'//ACHAR(48+IN(1))//ACHAR(48+IN(2))//ACHAR(48+IN(3)) &
         //'.OUT'
      ENDIF
      
      
      RETURN
      END    
      
      
!###################################################################################################
       SUBROUTINE PARTITION_CALCULATOR(NELST,NELEND,NELNOX,N_TO_BE_DISTRIBUTED,NPROCS,IPROC)
       IMPLICIT NONE
       INTEGER::N_TO_BE_DISTRIBUTED,NPROCS,IPROC,I1,I2,IA,NELST,NELNOX,NELEND



       I1=N_TO_BE_DISTRIBUTED/NPROCS+1
       I2=N_TO_BE_DISTRIBUTED/NPROCS

       IA=N_TO_BE_DISTRIBUTED-I2*NPROCS


       IF((IPROC+1).LE.IA) THEN
              NELST=IPROC*I1+1
              NELNOX=I1
       ELSE
              NELST=IA*I1+(IPROC-IA)*I2+1
              NELNOX=I2
       ENDIF

       NELEND=NELST+NELNOX-1

       RETURN
       END   



!###################################################################################################
      subroutine READGM(G0XYZ,IJK,NX,NELX,iproc_mpi)
      implicit none
      real(8), allocatable :: G0XYZ(:)
      integer, allocatable :: IJK(:)
      integer, intent(out) :: NX, NELX
      !locals
      integer :: IJKNEL(4)
      real(8) :: xc(3)
      integer :: MDIM,NDFADD,NDF,NODE,NGAUSS,LR,iproc_mpi
      character(len=50) :: fileName
      !dummies
      integer :: icyc,len,i,ii,j,jj,nnum,iElem,iNode,NEL,NEQ
      
      fileName = 'loadtime.inp'
      call openFile(201,fileName,iproc_mpi)
      call READSTR(201)
      read(201,*) icyc
      call READSTR(201)
      read(201,*) fileName
      close(201)
      
      len=len_trim(fileName)
      fileName=fileName(1:len)//'.inp'

      LR= 201
      call openFile(201,fileName,iproc_mpi)     

      call READSTR(LR)
      read(LR,*)MDIM,NDFADD
      NDF=MDIM+NDFADD     
      
      call READSTR(LR)
      read(LR,*) NX
      NEQ=NDF*NX

      allocate( G0XYZ(NX*MDIM) )
      G0XYZ = 0.d0
      ! READ INITIAL POSITIONS
      do iNode=1,NX
         read(LR,*) nnum,(xc(ii),ii=1,MDIM)
         do ii=1,MDIM
            G0XYZ(MDIM*(nnum-1)+ii)=xc(ii)
         enddo
      enddo
      
      call READSTR(LR)
      NGAUSS=1
      NODE=4
 

      read(LR,*) NELX

      allocate (IJK(NELX*NODE))
      IJK = 0
      ! READ IN THE IJK/CONNECTIVITY MATRIX
      jj=0
      do iElem=1,NELX
         read(LR,*) NEL,(IJKNEL(j),j=1,NODE)
         do j=1,NODE
            jj=jj+1
            IJK(jj)=IJKNEL(j)
         enddo
      enddo
      
      close(LR)

      return
      end
             

!###################################################################################################
      SUBROUTINE READSTR(LR)
      IMPLICIT REAL*8(A-H,O-Z)
! **** modified by deniz
      ! this code skips 1 line of text in the input file.
      character(150)::dummyString
      read(LR,'(A)') dummyString

      RETURN
      END  

!###################################################################################################
! output, rot,:
!   an active transformation matrix corresponding to the actual rotation due to euler angles {v}_i = TLGT_ij {v0}_j
!   as a (passive) coordinate transformation matrix, it transforms coordinates from crystal to global spatial frame. {v}_i = TLGT_ij {v_cry}_j
! input - euler angles
!     phi   - euler(1)
!     theta - euler(2)
!     omega - euler(3)
!---------------------------------------------------
      subroutine eulerRot(euler,rot)   ! restored to original
      implicit double precision (a-h,o-z)  

!      include 'aba_param.inc'

      dimension euler(3),tlg(3,3),rot(3,3)

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

      rot = TRANSPOSE(tlg)

      return
      end   
      