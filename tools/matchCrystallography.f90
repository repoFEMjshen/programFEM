      PROGRAM matchCrystallography
      
      use crystal_Titanium
      use grains, only : nGrains, &
         readGrainSizes_direct, readGrainPhases_direct, readTexture_direct, &
         grainSize,grainTexture, &
         grainSizesImported,textureImported
            
      implicit none

      ! element groupings
      character :: eulerDirection
      integer, parameter :: methodMisorientation = 2 ! 1=c-axis misorientation angle, 2=crystallographic misorientation acct for symmetries      
      real(8), allocatable :: grainTexture_Matched(:,:)
      real(8), allocatable :: grainTexture_AfterSwap(:,:)
      real(8), allocatable :: grainRot_Matched(:,:,:)
      real(8), allocatable :: grainRot_AfterSwap(:,:,:)
      real(8), allocatable :: listGrainBndAreas(:,:)
      real(8), allocatable :: histogramMisorientation_Matched(:)
      real(8), allocatable :: histogramMisorientation_AfterSwap(:)
      real(8), allocatable :: listGrainMisorientations_Matched(:,:)
      real(8), allocatable :: listGrainMisorientations_AfterSwap(:,:)
      
      real(8), allocatable :: histogramMisorientation_target(:)
      real(8), allocatable :: histogramMisorientation_RecBox(:)
      real(8) :: meanMis_target,meanMis_Matched,meanMis_AfterSwap
      real(8) :: T15_target,T15_Matched,T15_AfterSwap
      real(8) :: T15_RecBox,meanMis_RecBox
      
      real(8) :: errorHistogram_AfterSwap,errorHistogram_Matched
      real(8) :: errorMeanMis_AfterSwap,errorMeanMis_Matched
      real(8) :: errorT15_AfterSwap,errorT15_Matched
      
      integer :: quantityToMatch
      logical :: acceptSwap
      
      real(8), allocatable :: xMisorientation(:)
      real(8), parameter :: PI = 3.14159265359d0
      character(len=128) :: strMisorientationStatsFile
      character(len=128) :: strCommandLine
      character(len=150) :: strLine
      integer :: iColon
      real(8) :: rDummy
      real(8) :: misorientation
      real(8) :: errorNormalization,errorHistogram_Normalized
      real(8) :: error_target
      integer :: iNeighbor,jNeighbor
      integer :: nIterations, iIteration, nIntervals
      integer :: nSlipSystems
      real(8) :: SchmidTensors(3,3,nslip_max)
      
      logical :: success
      integer :: iError
            
      integer, allocatable :: nGrainNeighbors(:)
      integer, allocatable :: grainNeighborList(:,:)
      integer :: nGrainsInFile, nMaxGrainNeighbors
      integer :: staIndex, endIndex, elemIndex
      integer :: grainIdx, iGrainID, jGrainID, neighborGrainID 
      
      integer :: i,j
      
      ! initialize seed
      CALL init_random_seed()
      
      strMisorientationStatsFile = 'MisorientationStatistics.out'
      
      if (iargc() >= 1) then        ! User specified misorientation statistics file
         CALL getarg(1, strCommandLine)
         strMisorientationStatsFile = trim(strCommandLine)
      endif      
      
      !read grain sizes: grainSize(:)
      CALL readGrainSizes_direct('grainSizes.inp',nGrainsInFile)
      if(grainSizesImported) write(*,*) 'grain sizes imported for # grains:',nGrainsInFile
      if(.NOT.grainSizesImported) then
         write(*,*) 'WARNING ! no grain size information (grainSizes.inp) found. assuming default grain size of:',grainSize(1)
      endif

      !read texture: grainTexture(:)
      CALL readTexture_direct('grainTexture.inp',nGrainsInFile)
      if(textureImported) write(*,*) 'texture imported for # grains:',nGrainsInFile
      if(.NOT.textureImported) then
         write(*,*) 'WARNING ! no texture information (grainTexture.inp) found. assuming single crystal, euler:',grainTexture(:,1)
      endif
      
      ! reads grain neighbor and gr boundary areas of the CPFE model
      allocate(nGrainNeighbors(nGrains))
      call readGrainNeighborListsCPFE( &
            'grainNeighborList.dat','grainNeighborAreas.dat', &
            grainNeighborList,listGrainBndAreas, &
            nGrainNeighbors,nMaxGrainNeighbors,nGrains,.true.) ! read header only: nMaxGrainNeighbors
      allocate(listGrainBndAreas(nMaxGrainNeighbors,nGrains))
      allocate(grainNeighborList(nMaxGrainNeighbors,nGrains))
      call readGrainNeighborListsCPFE( &
            'grainNeighborList.dat','grainNeighborAreas.dat', &
            grainNeighborList,listGrainBndAreas, &
            nGrainNeighbors,nMaxGrainNeighbors,nGrains,.false.)
            
      ! initialize HCP slip system definitions
      call init_crys_HCP()      
      nSlipSystems = nslip_HCP

      ! calculate Schmid factors
      call create_schmidt(SchmidTensors,cst_HCP,cmt_HCP,nSlipSystems)

      write(*,*) 'reading misorientation distribution'
      
      open(102,file=trim(strMisorientationStatsFile))
      call readstr(102) ! '      fracGB_t15    meanMis'
      call readstr(102) ! 'CPFE: fracGB_t15    meanMis'
      call readStrReturn(102,strLine,iError) ! 'reconst.BOX: fracGB_t15    meanMis'
      iColon = index(strLine,':')
      read(strLine(iColon+1:len(strLine)),*) T15_RecBox, meanMis_RecBox
      call readStrReturn(102,strLine,iError) ! 'orig.BOX: fracGB_t15    meanMis'
      iColon = index(strLine,':')
      read(strLine(iColon+1:len(strLine)),*) T15_target, meanMis_target
      call readstr(102) ! '    misorientation histograms'
      call readstr(102) ! 'mis (deg)     CPFE    reconst.BOX    orig.BOX'
      read(102,*) nIntervals
      allocate(histogramMisorientation_RecBox(nIntervals))
      allocate(histogramMisorientation_target(nIntervals))
      do i = 1,nIntervals
         read(102,*) misorientation, histogramMisorientation_RecBox(i), histogramMisorientation_target(i)
      enddo
      close(102)
      ! normalize
      histogramMisorientation_RecBox = histogramMisorientation_RecBox / sum(histogramMisorientation_target)
      histogramMisorientation_target = histogramMisorientation_target / sum(histogramMisorientation_target)
   
      allocate(xMisorientation(nIntervals))
      do i=1,nIntervals
         xMisorientation(i) = (90.0/nIntervals)*(i-0.5)
      enddo
      
      !allocate the orientation list
      allocate(grainTexture_Matched(3,nGrains))
      
      grainTexture_Matched = grainTexture
      
      allocate(grainRot_Matched(3,3,nGrains))
      allocate(grainRot_AfterSwap(3,3,nGrains))
      allocate(histogramMisorientation_Matched(nIntervals))
      allocate(histogramMisorientation_AfterSwap(nIntervals))
      allocate(listGrainMisorientations_Matched(nMaxGrainNeighbors,nGrains))
      allocate(listGrainMisorientations_AfterSwap(nMaxGrainNeighbors,nGrains))
      
      do iGrainID = 1, nGrains
         call getRotationMatrix(grainRot_Matched(:,:,iGrainID), &
                              grainTexture_Matched(:,iGrainID))
      enddo
      
      write(*,*) 'calculating grain misorientation histograms...' 
      do iGrainID=1,nGrains
         do iNeighbor = 1,nGrainNeighbors(iGrainID)
            jGrainID = grainNeighborList(iNeighbor,iGrainID)
            if (methodMisorientation == 1) then ! c-axis misorientation
               call calcMisorientationAngle0001( &
                     grainRot_Matched(:,:,iGrainID), &
                     grainRot_Matched(:,:,jGrainID), &
                     misorientation)
            elseif (methodMisorientation == 2) then ! full crystallogprahic misorientation acct for symmetries
               call calcMisorientationAngleHCP( &
                     grainRot_Matched(:,:,iGrainID), &
                     grainRot_Matched(:,:,jGrainID), &
                     misorientation)
            endif
            misorientation = misorientation*180.d0/PI
            listGrainMisorientations_Matched(iNeighbor,iGrainID) = misorientation
         enddo
      enddo
      
      call calculateBoxMisorientationHistogram( &
         (/(iGrainID,iGrainID=1,nGrains)/), nGrains, &
         nGrainNeighbors(:), &
         listGrainMisorientations_Matched, &
         listGrainBndAreas, &
         nMaxGrainNeighbors,nGrains, &
         histogramMisorientation_Matched,nIntervals,meanMis_Matched,T15_Matched)
      histogramMisorientation_Matched = histogramMisorientation_Matched / sum(histogramMisorientation_Matched)
      
      write(*,*) 'Misorientation Distributions'
      write(*,*) 'theta   CPFE (initial)  target'
      do i=1,nIntervals
         write(*,'(F10.1,2F10.3)') xMisorientation(i),histogramMisorientation_Matched(i),histogramMisorientation_target(i)
      enddo
      write(*,*) ''
      write(*,'(A,F10.3)') 'CPFE mean misorientation (initial):', meanMis_Matched
      write(*,'(A,F10.3)') 'Target mean misorientation:', meanMis_target
      write(*,'(A,F10.3)') 'CPFE T15 (initial):', T15_Matched
      write(*,'(A,F10.3)') 'Target T15:', T15_target
      
      write(*,*) 'at each iteration, two grains will be randomly picked, &
                 &and orientations will be swapped, if this increases (or decreases)&
                 &the average misorientation'
      write(*,*) 'how many iterations to do?'
      write(*,'(A,I0,A)') '(better if in the order of: nGrains*nGrains=',nGrains*nGrains,')'
      write(*,*) '(or, enter -1 to specify target relative error)'
      read(*,*) nIterations
      
      error_target = 0.d0
      if (nIterations == -1) then
         write(*,*) 'enter target relative error to do as many iterations as necessary:'
         read(*,*) error_target
         nIterations = nGrains*nGrains
      endif
      
      write(*,*) 'Which misorientation measure to match?'
      write(*,*) '1: Full misorientation histogram'
      write(*,*) '2: Match mean misorientation'
      write(*,*) '3: T15 = Fraction GB < 15 deg misorientation'
      read(*,*) quantityToMatch

      write(*,*) 'starting iterations...'
            
      do iIteration = 1,nIterations
         
         CALL getRandomGrainID(iGrainID,nGrains)
         CALL getRandomGrainID(jGrainID,nGrains)
         
         if (jGrainID == iGrainID) cycle
         
         ! determine whether we should swap iGrain with jGrain
         grainRot_AfterSwap = grainRot_Matched
         grainRot_AfterSwap(:,:,iGrainID) = grainRot_Matched(:,:,jGrainID)
         grainRot_AfterSwap(:,:,jGrainID) = grainRot_Matched(:,:,iGrainID)
         grainTexture_AfterSwap = grainTexture_Matched
         grainTexture_AfterSwap(:,iGrainID) = grainTexture_Matched(:,jGrainID)
         grainTexture_AfterSwap(:,jGrainID) = grainTexture_Matched(:,iGrainID)
         listGrainMisorientations_AfterSwap = listGrainMisorientations_Matched
         ! update the misorientation angles around those grains that were swapped
         ! .. around grain-I
         do iNeighbor = 1,nGrainNeighbors(iGrainID)
            neighborGrainID = grainNeighborList(iNeighbor,iGrainID)

            if (methodMisorientation == 1) then ! c-axis misorientation
               call calcMisorientationAngle0001( &
                  grainRot_AfterSwap(:,:,iGrainID), &
                  grainRot_AfterSwap(:,:,neighborGrainID), &
                  misorientation)         
            elseif (methodMisorientation == 2) then ! full crystallogprahic misorientation acct for symmetries
               call calcMisorientationAngleHCP( &
                  grainRot_AfterSwap(:,:,iGrainID), &
                  grainRot_AfterSwap(:,:,neighborGrainID), &
                  misorientation)       
            else
               write(*,*) 'unrecognised methodMisorientation'
               stop
            endif
            misorientation = misorientation*180.d0/PI

            listGrainMisorientations_AfterSwap(iNeighbor,iGrainID) = misorientation
            do jNeighbor = 1,nGrainNeighbors(neighborGrainID)
               if (grainNeighborList(jNeighbor,neighborGrainID)==iGrainID) then
                  listGrainMisorientations_AfterSwap(jNeighbor,neighborGrainID) = misorientation
                  exit
               endif
            enddo
         enddo
         ! .. around grain-J
         do iNeighbor = 1,nGrainNeighbors(jGrainID)
            neighborGrainID = grainNeighborList(iNeighbor,jGrainID)

            if (methodMisorientation == 1) then ! c-axis misorientation

               call calcMisorientationAngle0001( &
                     grainRot_AfterSwap(:,:,jGrainID), &
                     grainRot_AfterSwap(:,:,neighborGrainID), &
                     misorientation)

            elseif (methodMisorientation == 2) then ! full crystallogprahic misorientation acct for symmetries
               call calcMisorientationAngleHCP( &
                     grainRot_AfterSwap(:,:,jGrainID), &
                     grainRot_AfterSwap(:,:,neighborGrainID), &
                     misorientation)       
            else
               write(*,*) 'unrecognised methodMisorientation'
               stop
            endif
            misorientation = misorientation*180.d0/PI

            listGrainMisorientations_AfterSwap(iNeighbor,jGrainID) = misorientation
            do jNeighbor = 1,nGrainNeighbors(neighborGrainID)
               if (grainNeighborList(jNeighbor,neighborGrainID)==jGrainID) then
                  listGrainMisorientations_AfterSwap(jNeighbor,neighborGrainID) = misorientation
                  exit
               endif
            enddo
         enddo
         
         ! update histogram
         call calculateBoxMisorientationHistogram( &
         (/(iGrainID,iGrainID=1,nGrains)/), nGrains, &
         nGrainNeighbors(:), &
         listGrainMisorientations_AfterSwap, &
         listGrainBndAreas, &
         nMaxGrainNeighbors,nGrains, &
         histogramMisorientation_AfterSwap,nIntervals,meanMis_AfterSwap,T15_AfterSwap)
         histogramMisorientation_AfterSwap = histogramMisorientation_AfterSwap / sum(histogramMisorientation_AfterSwap)

         errorHistogram_AfterSwap = dot_product( &
               histogramMisorientation_AfterSwap - histogramMisorientation_target, &
               histogramMisorientation_AfterSwap - histogramMisorientation_target)
         errorHistogram_Matched = dot_product( &
               histogramMisorientation_Matched - histogramMisorientation_target, &
               histogramMisorientation_Matched - histogramMisorientation_target)
         
         errorMeanMis_AfterSwap = abs(meanMis_AfterSwap - meanMis_target)
         errorMeanMis_Matched = abs(meanMis_Matched - meanMis_target)
         
         errorT15_AfterSwap = abs(T15_AfterSwap - T15_target)
         errorT15_Matched = abs(T15_Matched - T15_target)

         acceptSwap = .false.
         if (quantityToMatch == 1) then ! error measure is on the misorientation histogram
            if (errorHistogram_AfterSwap < errorHistogram_Matched) acceptSwap = .true.
         elseif (quantityToMatch == 2) then ! error measure is on the mean misorientation
            if (errorMeanMis_AfterSwap < errorMeanMis_Matched) acceptSwap = .true.
         elseif (quantityToMatch == 3) then ! error measure is on t15
            if (errorT15_AfterSwap < errorT15_Matched) acceptSwap = .true.
         endif
         
         if (acceptSwap) then
            grainRot_Matched = grainRot_AfterSwap 
            grainTexture_Matched = grainTexture_AfterSwap
            listGrainMisorientations_Matched = listGrainMisorientations_AfterSwap
            histogramMisorientation_Matched = histogramMisorientation_AfterSwap
            meanMis_Matched = meanMis_AfterSwap
            T15_Matched = T15_AfterSwap
            
            errorHistogram_Matched = errorHistogram_AfterSwap
            errorMeanMis_Matched = errorMeanMis_AfterSwap
            errorT15_Matched = errorT15_AfterSwap
         endif

         ! evaluate stopping criterion
         if (error_target /= 0.d0) then
            errorNormalization = dot_product(histogramMisorientation_target,histogramMisorientation_target)
            errorHistogram_Normalized = sqrt(errorHistogram_Matched / errorNormalization)
            if (quantityToMatch == 1) then   ! histogram
               if (errorHistogram_Normalized < error_target) then
                  exit ! reached target error on misorientation distribution
               endif
            elseif (quantityToMatch == 2) then ! mean misorientation
               if (errorMeanMis_Matched/(PI/2) < error_target) then
                  exit ! reached target error on mean misorientation
               endif
            elseif (quantityToMatch == 3) then ! T15
               if (errorT15_Matched < error_target) then
                  exit ! reached target error on mean misorientation
               endif

            endif
         endif
            
      enddo
      
      errorNormalization = dot_product(histogramMisorientation_target,histogramMisorientation_target)
      errorHistogram_Normalized  = sqrt(errorHistogram_Matched / errorNormalization)
      write(*,*) 'final histogram error:',errorHistogram_Normalized
      write(*,'(A,F10.3)') 'final relative error in mean misorientation:',errorMeanMis_Matched/xMisorientation(nIntervals)
      write(*,*) 'number of iterations:',iIteration
      
      write(*,*) 'done.'
      write(*,*) 'theta   CPFE (matched)  target'
      open(131,file='MisorientationStatistics_Matched.out')
      
      
      write(131,*) '      fracGB_t15    meanMis'
      write(131,*) 'CPFE:',T15_Matched,meanMis_Matched
      write(131,*) 'reconst.BOX:',T15_RecBox, meanMis_RecBox
      write(131,*) 'orig.BOX (target):',T15_target, meanMis_target
      write(131,*) '    misorientation histograms'
      write(131,*) 'mis (deg)     CPFE    reconst.BOX    orig.BOX (target)'
      write(*,*) 'mis (deg)     CPFE    reconst.BOX    orig.BOX (target)'
      write(131,*) nIntervals
      do i=1,nIntervals
         write(131,'(F10.1,3F10.4)') xMisorientation(i), &
                     histogramMisorientation_Matched(i), &
                     histogramMisorientation_RecBox(i), &
                     histogramMisorientation_target(i)
         write(*,'(F10.1,3F10.4)') xMisorientation(i), &
                     histogramMisorientation_Matched(i), &
                     histogramMisorientation_RecBox(i), &
                     histogramMisorientation_target(i)
      enddo
      close(131)
      write(*,*) ''
      write(*,'(A,F10.3)') 'target mean misorientation: ',meanMis_target
      write(*,'(A,F10.3)') 'CPFE mean misorientation: ',meanMis_Matched
      write(*,'(A,F10.3)') 'target T15: ',T15_target
      write(*,'(A,F10.3)') 'CPFE T15: ',T15_Matched
      
      open(132,file='grainTexture_Matched.inp')
      write(132,*) nGrains
      do iGrainID=1,nGrains
         write(132,*) iGrainID,grainTexture_Matched(:,iGrainID)
      enddo
      close(132)
      
      END PROGRAM

      subroutine getRandomGrainID(grainID,nGrains)
      
      implicit none
      
      integer, intent(out):: grainID
      integer, intent(in) :: nGrains
      real(8) :: rndNum
      
      CALL RANDOM_NUMBER(rndNum)
         
      grainID = floor(rndNum*nGrains) + 1
      
      ! fix the (very unlikely) limit cases
      if (grainID > nGrains) grainID = nGrains  
      if (grainID < 1) grainID = 1
      
      end subroutine
      
      SUBROUTINE cross(vc, v1, v2)
      implicit none
      real(8), intent(in) :: v1(3),v2(3)
      real(8), intent(out):: vc(3)
      
      vc(1)=+v1(2)*v2(3)-v1(3)*v2(2)
      vc(2)=-v1(1)*v2(3)+v1(3)*v2(1)
      vc(3)=+v1(1)*v2(2)-v1(2)*v2(1)
      
      END SUBROUTINE

      SUBROUTINE readStr(LR)
      IMPLICIT REAL*8(A-H,O-Z)
      !old version:
      !DIMENSION FLAG(20)
      !READ(LR,505) (FLAG(I),I=1,20)
!505   FORMAT(20A4)

! **** modified by deniz
      ! this code skips 1 line of text in the input file.
      character(1024)::dummyString
      read(LR,'(A)') dummyString

      RETURN
      END SUBROUTINE
      
      SUBROUTINE readStrReturn(fileNum,lineStr,error)
      ! this code skips 1 line of text in the input file.
      integer, intent(in) :: fileNum
      character(150), intent(out) :: lineStr
      integer, intent(out) :: error
      read(fileNum,'(A)',iostat=error) lineStr
      
      END
      
      subroutine create_schmidt(sch,cst,cmt,nslip)
      implicit none
      
      real(8),intent(out):: sch(3,3,nslip)
      real(8),intent(in) :: cmt(3,nslip),cst(3,nslip)
      integer,intent(in) :: nslip
      
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