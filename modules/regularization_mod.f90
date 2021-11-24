

      module regularization
      
      use meshtools
      
      implicit none
      
      
      integer, allocatable :: elemNonLocalList(:)  
      integer, allocatable :: elemIndexToElemID(:)  
      integer, allocatable :: elemNonLocalListIndex(:)
      integer, allocatable :: elemNonLocalCount(:)   
      real(8), allocatable :: elemNonLocalWeights(:)
      integer :: nTotNonLocalElementsProc
      integer, private :: reg_nElems, reg_nElemsProc
      real(8), private :: reg_CharacteristicLength
      logical, private :: reg_successList, reg_successWrite
      
      integer, parameter :: reg_GRAINAVG = 1
      integer, parameter :: reg_GAUSSIAN = 2
      integer, parameter :: reg_RVEAVG = 3
      
      contains
      
      SUBROUTINE reg_destruct()
      implicit none
            if (allocated(elemNonLocalList)) deallocate(elemNonLocalList)
            if (allocated(elemNonLocalListIndex)) deallocate(elemNonLocalListIndex)
            if (allocated(elemNonLocalCount)) deallocate(elemNonLocalCount)
            if (allocated(elemIndexToElemID)) deallocate(elemIndexToElemID)
            if (allocated(elemNonLocalWeights)) deallocate(elemNonLocalWeights)
      END SUBROUTINE
      
      SUBROUTINE regularizeFieldGaussian(fieldLocal, nDOF, fieldRegularized, success)
      implicit none
      
      integer, intent(in) :: nDOF
      real(8), intent(in) :: fieldLocal(nDOF,reg_nElems)
      real(8), intent(out):: fieldRegularized(nDOF,reg_nElems)
      logical, intent(out):: success
      
      integer :: iElemIdxProc,jElemIdxProc,iElemID,jElemID,elemIdx,nElemCount
      integer :: staElemIdx,endElemIdx
      real(8) :: weightTimesVol
      
      success = .false.
      if (.not.reg_successList) return
      
      ! reset field
      fieldRegularized = 0.d0

      do iElemIdxProc=1,reg_nElemsProc
      
         nElemCount = elemNonLocalCount(iElemIdxProc)
         iElemID = elemIndexToElemID(iElemIdxProc)
         staElemIdx = 1
         if (iElemIdxProc /= 1) staElemIdx = elemNonLocalListIndex(iElemIdxProc-1)
         endElemIdx = elemNonLocalListIndex(iElemIdxProc)-1
         do elemIdx=staElemIdx,endElemIdx
            jElemID = elemNonLocalList(elemIdx)
            weightTimesVol = elemNonLocalWeights(elemIdx)
            fieldRegularized(:,iElemID) = fieldRegularized(:,iElemID) + fieldLocal(:,jElemID)*weightTimesVol
         enddo

      enddo
      
      success = .true.
         
      END SUBROUTINE
      
      SUBROUTINE calculateRveAverage(fieldLocal, nDOF, fieldRegularized, &
                                     rveAverage, &
                                     nElems,NX,IJK,G0XYZ,&
                                     success)
      implicit none
      
      integer, intent(in) :: nDOF
      real(8), intent(in) :: fieldLocal(nDOF,nElems)
      real(8), intent(out):: fieldRegularized(nDOF,nElems)
      real(8), intent(out):: rveAverage(nDOF)
      integer, intent(in) :: nElems,NX
      integer, intent(in) :: IJK(4*nElems)
      real(8), intent(in) :: G0XYZ(3*NX)
      logical, intent(out):: success
      
      ! locals
      real(8) :: elemVol(nElems)
      real(8) :: xNodes(3,4), elemJacMat(3,3), elemJac
      real(8) :: rveVol
      integer :: iElem
      
      success = .false.

      ! calculate element volumes
      do iElem=1,nElems
         CALL getElemXNodes(iElem,IJK,G0XYZ,4,xNodes)
         CALL calcJacTET4(xNodes, elemJacMat, elemJac, elemVol(iElem))   
      enddo
      
      ! first calculate domain average for this time step
      fieldRegularized = 0.d0
      rveAverage = 0.d0
      rveVol = 0.d0
      do iElem = 1, nElems      
         rveAverage(:) = rveAverage(:) + fieldLocal(:,iElem)*elemVol(iElem)
         rveVol = rveVol + elemVol(iElem)
      enddo
      rveAverage(:) = rveAverage(:) / rveVol
      
      do iElem = 1, nElems
         fieldRegularized(:,iElem) = rveAverage(:)
      enddo
      
      success = .true.
      
      END SUBROUTINE

      SUBROUTINE regularizeFieldGrainAVG(fieldLocal, nDOF, fieldRegularized, &
                                         fieldGrains, &
                                         nElems,NX,IJK,G0XYZ,&
                                         grainIDelem,grainElementIndx,grainElements,nGrains,&
                                         staGrainIdx,endGrainIdx, &
                                         success)
      implicit none
      
      integer, intent(in) :: nDOF
      real(8), intent(in) :: fieldLocal(nDOF,nElems)
      real(8), intent(out):: fieldGrains(nDOF,nGrains)
      real(8), intent(out):: fieldRegularized(nDOF,nElems)
      integer, intent(in) :: nElems,NX,nGrains
      integer, intent(in) :: IJK(4*nElems)
      real(8), intent(in) :: G0XYZ(3*NX)
      integer, intent(in) :: grainIDelem(nElems),grainElements(nElems)
      integer, intent(in) :: grainElementIndx(nGrains)
      integer, intent(in) :: staGrainIdx,endGrainIdx
      logical, intent(out):: success
      real(8) :: elemVol(nElems)
      real(8) :: xNodes(3,4), elemJacMat(3,3), elemJac
      real(8) :: grainVolume
      
      integer :: iElem,jElem,iGrain
      integer :: elemIdx,staIndex,endIndex
      
      success = .false.
      
      ! calculate element volumes
      do iElem=1,nElems
         CALL getElemXNodes(iElem,IJK,G0XYZ,4,xNodes)
         CALL calcJacTET4(xNodes, elemJacMat, elemJac, elemVol(iElem))   
      enddo
      
      ! reset field
      fieldRegularized = 0.d0
      fieldGrains = 0.d0
      
      if (staGrainIdx==0 .and. endGrainIdx==0) then
         success = .true.         
         return
      endif

      do iGrain = staGrainIdx,endGrainIdx
      
         fieldGrains(:,iGrain) = 0.d0  ! reset grain-averaging variable
         grainVolume = 0.d0      ! reset grain volume - used for normalization
         staIndex = 1
         if (iGrain.NE.1) staIndex = grainElementIndx(iGrain-1)
         endIndex = grainElementIndx(iGrain)-1
         do elemIdx=staIndex,endIndex
            iElem = grainElements(elemIdx)
            fieldGrains(:,iGrain) = fieldGrains(:,iGrain) + fieldLocal(:,iElem)*elemVol(iElem)
            grainVolume = grainVolume + elemVol(iElem)
         enddo
         fieldGrains(:,iGrain) = fieldGrains(:,iGrain) / grainVolume
         
         ! now assign the grain-averaged field to the grain elements
         do elemIdx=staIndex,endIndex
            iElem = grainElements(elemIdx)
            fieldRegularized(:,iElem) = fieldGrains(:,iGrain)
         enddo
      enddo
      
      success = .true.
         
      END SUBROUTINE
      SUBROUTINE generateElementNonLocalList(characteristicLength,grainSize, &
                                                         nElems,NX,IJK,G0XYZ, &
                          grainIDelem,grainElementIndx,grainElements,nGrains, &
                                                     staGrainIdx,endGrainIdx, &
                                                          reportTo_l,success)
      implicit none

      real(8), intent(in) :: characteristicLength
      integer, intent(in) :: nElems,NX,nGrains
      integer, intent(in) :: IJK(4*nElems)
      real(8), intent(in) :: G0XYZ(3*NX)
      integer, intent(in) :: grainIDelem(nElems),grainElements(nElems)
      integer, intent(in) :: grainElementIndx(nGrains)
      real(8), intent(in) :: grainSize(nGrains)
      integer, intent(in) :: staGrainIdx,endGrainIdx
      integer, intent(in) :: reportTo_l
      logical, intent(out):: success
      
      integer :: reportTo
      integer :: iElem,jElem,iElemIdx,jElemIdx,iGrain,thisIndex
      real(8) :: elemVol(nElems)
      real(8) :: elemCenters(3,nElems)
      integer :: staIndex,endIndex
      real(8) :: weight, normalizationConstant, dist, distVec(3)
      real(8) :: xNodes(3,4), elemJacMat(3,3), elemJac

      integer :: staElemIDProc,endElemIDProc,staNonLocalElemIdx,endNonLocalElemIdx
      integer :: nElemsProc,iElemsProc
      
      reportTo = reportTo_l
      if (reportTo == 0) reportTo=6 ! print to screen
      
      reg_successList = .false.
      success = .false.
      
      if (staGrainIdx==0 .and. endGrainIdx==0) then
         reg_CharacteristicLength = characteristicLength
         reg_nElems = nElems
         reg_nElemsProc = 0
         reg_successList = .true.
         success = .true.         
         return
      endif
      
      do iElem=1,nElems
         CALL getElemXNodes(iElem,IJK,G0XYZ,4,xNodes)
         CALL getElemXCenter(iElem,IJK,G0XYZ,4,elemCenters(:,iElem))
         CALL calcJacTET4(xNodes, elemJacMat, elemJac, elemVol(iElem))   
      enddo
      
      ! scanning to determine the number of non-local elements per element
      if (reportTo > 0) write(reportTo,'(A)',advance='no') 'scanning...'

      staElemIDProc = 1
      if (staGrainIdx.NE.1) staElemIDProc = grainElementIndx(staGrainIdx-1)
      endElemIDProc = grainElementIndx(endGrainIdx)-1
      nElemsProc = endElemIDProc-staElemIDProc+1  ! # of elements handled by this Processor
      ! allocate the array for # of nonlocal elements and list indices for the elements of this processor
      allocate(elemNonLocalCount(nElemsProc))
      allocate(elemIndexToElemID(nElemsProc))
      allocate(elemNonLocalListIndex(nElemsProc))
      elemNonLocalCount = 0
      elemIndexToElemID = 0
      elemNonLocalListIndex = 0
      iElemsProc = 0           ! reset the index for the elements handled by this Processor
      nTotNonLocalElementsProc = 0 ! reset the index for the non-local couplings handled by this Processor
      do iGrain=staGrainIdx,endGrainIdx
         staIndex = 1
         if (iGrain.NE.1) staIndex = grainElementIndx(iGrain-1)
         endIndex = grainElementIndx(iGrain)-1
         do iElemIdx=staIndex,endIndex
            iElem = grainElements(iElemIdx)
            iElemsProc = iElemsProc + 1
            elemIndexToElemID(iElemsProc)=iElem
            do jElemIdx=staIndex,endIndex
               jElem = grainElements(jElemIdx)
               distVec(:) = elemCenters(:,iElem) - elemCenters(:,jElem)
               dist = dsqrt(DOT_PRODUCT(distVec,distVec))
               if (characteristicLength > 0.d0) then
                  weight = dexp(-(dist / characteristicLength)**2.d0)
               elseif (characteristicLength == 0.d0) then
                  weight = dexp(-(dist / (grainSize(iGrain) / 2.d0))**2.d0)
               elseif (characteristicLength < 0.d0) then
                  weight = dexp(-(dist / (-grainSize(iGrain) * characteristicLength))**2.d0)
               endif
               if (weight > 0.05d0) then
                  nTotNonLocalElementsProc = nTotNonLocalElementsProc + 1
                  elemNonLocalCount(iElemsProc) = elemNonLocalCount(iElemsProc) + 1
                  elemNonLocalListIndex(iElemsProc) = nTotNonLocalElementsProc + 1
               endif
            enddo
         enddo
      enddo
      
      if (reportTo > 0) write(reportTo,*) 'done. total number of non-local couplings:',nTotNonLocalElementsProc
      if (reportTo > 0) write(reportTo,*) 'generating non-local element list...'
      
      ! allocate the array for list indices of the elements of this processor
      allocate(elemNonLocalList(nTotNonLocalElementsProc))
      allocate(elemNonLocalWeights(nTotNonLocalElementsProc))
      elemNonLocalList = 0
      elemNonLocalWeights = 0.d0
      
      iElemsProc = 0           ! reset the index for the elements handled by this Processor
      nTotNonLocalElementsProc = 0 ! reset the index for the non-local couplings handled by this Processor
      do iGrain=staGrainIdx,endGrainIdx
         staIndex = 1
         if (iGrain.NE.1) staIndex = grainElementIndx(iGrain-1)
         endIndex = grainElementIndx(iGrain)-1
         do iElemIdx=staIndex,endIndex
            iElem = grainElements(iElemIdx)
            iElemsProc = iElemsProc + 1
            normalizationConstant = 0.d0
            do jElemIdx=staIndex,endIndex
               jElem = grainElements(jElemIdx)
               distVec(:) = elemCenters(:,iElem) - elemCenters(:,jElem)
               dist = dsqrt(DOT_PRODUCT(distVec,distVec))
               if (characteristicLength > 0.d0) then
                  weight = dexp(-(dist / characteristicLength)**2.d0)
               elseif (characteristicLength == 0.d0) then
                  weight = dexp(-(dist / (grainSize(iGrain) / 2.d0))**2.d0)
               elseif (characteristicLength < 0.d0) then
                  weight = dexp(-(dist / (-grainSize(iGrain) * characteristicLength))**2.d0)
               endif
               if (weight > 0.05d0) then
                  nTotNonLocalElementsProc = nTotNonLocalElementsProc + 1
                  elemNonLocalList(nTotNonLocalElementsProc) = jElem
                  elemNonLocalWeights(nTotNonLocalElementsProc) = weight*elemVol(jElem)
                  normalizationConstant = normalizationConstant + weight*elemVol(jElem)
               endif
            enddo
            
            staNonLocalElemIdx = 1
            if (iElemsProc /= 1) staNonLocalElemIdx = elemNonLocalListIndex(iElemsProc-1)
            endNonLocalElemIdx = elemNonLocalListIndex(iElemsProc) - 1
            
            if (endNonLocalElemIdx-staNonLocalElemIdx+1 /= elemNonLocalCount(iElemsProc)) then
               write(*,*) 'assertion error: nonLocalElem indexing'
               stop
            endif
            
            if (reportTo > 0) write(reportTo,*) 'found ',elemNonLocalCount(iElemsProc),'elements around',iElem
            
            elemNonLocalWeights(staNonLocalElemIdx:endNonLocalElemIdx) = &
               elemNonLocalWeights(staNonLocalElemIdx:endNonLocalElemIdx) / normalizationConstant
         enddo
      enddo
         
      if (reportTo > 0) write(reportTo,*) 'done'
      
      reg_CharacteristicLength = characteristicLength
      reg_nElems = nElems
      reg_nElemsProc = nElemsProc
      reg_successList = .true.
      success = .true.
      
      END SUBROUTINE
      
      SUBROUTINE readElementNonLocalList(strNonlocalFile,NELX,success)
      implicit none

      character(len=*), intent(in) :: strNonlocalFile
      integer, intent(in) :: NELX
      logical, intent(out) :: success

      !Local
      logical :: fileExists
      integer :: error, iElemIdxProc, elemID, nElemCount, iNode, nElems, nElemsProc, I, lastIndex

      reg_successList = .false.
      success = .false.
      
      inquire(file=trim(strNonlocalFile),exist=fileExists)
      if(.NOT.fileExists) then
         return
      endif
      
      open(716,file=trim(strNonlocalFile))
      read(716,*,iostat=error) reg_CharacteristicLength, nElems
      read(716,*,iostat=error) nElemsProc, nTotNonLocalElementsProc
      if (error.ne.0 .or. nElemsProc.NE.NELX) then
         close(716)
         return
      endif
      
      ! allocate data structures: list and weights of element-element couplings for the averaging operation
      allocate(elemNonLocalList(nTotNonLocalElementsProc))
      allocate(elemIndexToElemID(nElemsProc))
      allocate(elemNonLocalCount(nElemsProc))
      allocate(elemNonLocalWeights(nTotNonLocalElementsProc))
      elemNonLocalList = 0
      elemNonLocalCount = 0
      elemNonLocalWeights = 0.d0
            
      lastIndex = 1
      do iElemIdxProc=1,nElemsProc
      
         read(716,*,iostat=error) elemID, nElemCount
         elemNonLocalCount(iElemIdxProc) = nElemCount
         elemIndexToElemID(iElemIdxProc) = elemID
         read(716,*,iostat=error) elemNonLocalList(lastIndex:lastIndex+nElemCount-1)
         read(716,*,iostat=error) elemNonLocalWeights(lastIndex:lastIndex+nElemCount-1)
         lastIndex = lastIndex+nElemCount
                                  
      enddo
      
      if (error.ne.0) then
         close(716)
         return
      endif
      
      reg_nElemsProc = nElemsProc
      reg_nElems = nElems
      success = .true.
      reg_successList = .true.
      close(716)
      
      END SUBROUTINE
      
      SUBROUTINE writeElementNonLocalList(strNonlocalFile,NELX,success)
      implicit none

      character(len=*), intent(in) :: strNonlocalFile      
      integer, intent(in) :: NELX
      logical, intent(out) :: success

      !Local
      !!!--- logical :: fileExists
      integer :: error, iElemIdxProc, jElem, nElemCount, iNode, I, lastIndex

      if (.not.reg_successList) then
         success = .false.
         reg_successWrite = .false.
         return
      endif
                  
      open(716,file=trim(strNonlocalFile),iostat=error)
      if (error.ne.0) then
         success = .false.
         reg_successWrite = .false.         
         return
      endif
      write(716,*) reg_CharacteristicLength, reg_nElems
      write(716,*) reg_nElemsProc, nTotNonLocalElementsProc
      
      lastIndex = 1
      do iElemIdxProc=1,reg_nElemsProc
      
         nElemCount = elemNonLocalCount(iElemIdxProc)
         write(716,*) elemIndexToElemID(iElemIdxProc),nElemCount
         write(716,*) elemNonLocalList(lastIndex:lastIndex+nElemCount-1)
         write(716,*) elemNonLocalWeights(lastIndex:lastIndex+nElemCount-1)
         lastIndex = lastIndex+nElemCount
      enddo
      
      success = .true.
      reg_successWrite = .true.
      close(716)
      
      END SUBROUTINE

      
      END MODULE